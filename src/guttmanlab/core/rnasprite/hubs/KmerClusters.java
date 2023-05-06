package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.BinomialDistribution;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.rnasprite.Kmer;
import guttmanlab.core.rnasprite.RNAInterval;

public class KmerClusters {
	int totalSize;
	int numPerm=5;
	IntervalTree<String> probabilityTree;
	
	
	public KmerClusters(BarcodingDataStreaming data, Collection<Kmer> hubs, Map<String, String> allRNA, String save, int k) throws IOException {
		Collection<Kmer> kmers=getKmers(hubs, k-1);
		writeAnnotation(save+".row.annotation", kmers);
		writeColumnAnnotation(save+".column.annotation", hubs);
		
		MatrixWithHeaders observed=makeMatrix(kmers, allRNA);
		
		this.probabilityTree=data.getProbabilityTree(allRNA);
		this.totalSize=probabilityTree.max().getEnd();
		
		Map<Kmer, Collection<Cluster>> counts=new TreeMap<Kmer, Collection<Cluster>>();
		
		
		Collection<String> names=new TreeSet<String>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			names.addAll(c.getRNANames());
			score(observed, c, allRNA, k, counts);
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		write(save+".names", names);
		
		data.close();
		
		observed.write(save+".observed");
		
		//TODO Score conditional
		MatrixWithHeaders expectedConditional=scoreConditional(observed, counts);
		
		expectedConditional=minEnrichment(expectedConditional);
		//expectedConditional=subset(expectedConditional, 2.0);
		mask(expectedConditional);
		
		expectedConditional.write(save);
		
	}
	
	
	private void write(String string, Collection<String> names) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(String name:names) {writer.write(name+"\n");}
		
		writer.close();
	}


	private void writeColumnAnnotation(String string, Collection<Kmer> kmers) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		writer.write("Hub\tRNA\n");
		for(Kmer k: kmers) {
			for(String rna:k.getRegions()) {
				writer.write(k.getName()+"\t"+rna+"\n");
			}
		}
		
		writer.close();
		
	}


	private void writeAnnotation(String string, Collection<Kmer> kmers) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		writer.write("Hub\tKmer\n");
		for(Kmer k: kmers) {
			writer.write(k.getName()+"\t"+k.toString()+"\n");
		}
		
		writer.close();
	}


	private void mask(MatrixWithHeaders expectedConditional) {
		for(String row: expectedConditional.getRowNames()) {
			for(String col: expectedConditional.getColumnNames()) {
				Kmer k=new Kmer(row);
				if(k.containsRegion(col)) {expectedConditional.set(row, col, -1);}
			}
		}
	}


	private MatrixWithHeaders minEnrichment(MatrixWithHeaders expectedConditional) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(expectedConditional.getRowNames(), expectedConditional.getColumnNames());
		
		for(String row: expectedConditional.getRowNames()) {
			for(String column: expectedConditional.getColumnNames()) {
				double min=getMinEnrichment(row, column, expectedConditional);
				rtrn.set(row, column, min);
			}
		}
		
		return rtrn;
	}

	private double getMinEnrichment(String row, String column, MatrixWithHeaders expectedConditional) {
		Kmer kmer=new Kmer(row);
		kmer.addRegion(column);
		
		double[] vals=new double[kmer.getSize()];
		int count=0;
		for(String region: kmer.getRegionList()) {
			Kmer k=new Kmer();
			k.addRegions(kmer.getRegionList());
			String newRow=k.remove(region).toString();
			double val=0;
			if(expectedConditional.containsRow(newRow) && expectedConditional.containsColumn(region)) {val=expectedConditional.get(newRow, region);}
			vals[count]=val;
			count++;
		}
		return Statistics.min(vals);
	}

	private MatrixWithHeaders subset(MatrixWithHeaders observed) {
		List<String> rows=new ArrayList<String>();
		
		for(String row: observed.getRowNames()) {
			double max=Statistics.max(observed.getRow(row));
			if(max>0) {rows.add(row);}
		}
		
		System.err.println(observed.getRowNames().size()+" "+rows.size());
		
		return observed.submatrixByRowNames(rows);
	}
	
	private MatrixWithHeaders subset(MatrixWithHeaders observed, double val) {
		List<String> rows=new ArrayList<String>();
		
		for(String row: observed.getRowNames()) {
			double max=Statistics.max(observed.getRow(row));
			if(max>val) {rows.add(row);}
		}
		
		System.err.println(observed.getRowNames().size()+" "+rows.size());
		
		return observed.submatrixByRowNames(rows);
	}

	private MatrixWithHeaders scoreConditional(MatrixWithHeaders observed, Map<Kmer, Collection<Cluster>> counts) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(observed.getRowNames(), observed.getColumnNames());
		
		int counter=0;
		for(String kmer: observed.getRowNames()) {
			System.err.println(counter+" "+observed.getRowNames().size());
			Kmer k=new Kmer(kmer);
			Collection<Cluster> clusters=counts.get(k);
			if(counts.containsKey(k) && !clusters.isEmpty()) {
				Map<String, Double> permScoresByRegion=perm(clusters, k, numPerm);
				for(String col: permScoresByRegion.keySet()) {
					if(observed.containsColumn(col)) {
						double score=observed.get(kmer, col);
						double ratio=0;
						if(score>0) {ratio=score/permScoresByRegion.get(col);}
						rtrn.set(kmer, col, ratio);
					}
				}
			}
			counter++;
		}
		
		return rtrn;
	}

	private Map<String, Double> perm(Collection<Cluster> clusters, Kmer k, int numPerm2) {
		
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(int i=0; i<numPerm2; i++) {
			Map<String, Integer> regions=perm(clusters, k);
			update(regions, rtrn);
		}
		
		rtrn=divide(rtrn, numPerm2);
		return rtrn;
	}

	

	private Map<String, Integer> perm(Collection<Cluster> clusters, Kmer k) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		for(Cluster c: clusters) {
			int size=c.size()-k.getRegions().size();
			Cluster rand=this.getPerm(size, probabilityTree, totalSize);
			update(rand, rtrn);
		}
		return rtrn;
	}

	private void update(Cluster rand, Map<String, Integer> rtrn) {
		for(String region: rand.getRNANames()) {
			int count=0;
			if(rtrn.containsKey(region)) {count=rtrn.get(region);}
			count++;
			rtrn.put(region, count);
		}
		
	}

	private void update(Map<String, Integer> regions, Map<String, Double> rtrn) {
		for(String region: regions.keySet()) {
			double count=0;
			if(rtrn.containsKey(region)) {count=rtrn.get(region);}
			count+=regions.get(region);
			rtrn.put(region, count);
		}
		
	}

	private Map<String, Double> divide(Map<String, Double> map, int numPerm2) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		for(String key: map.keySet()) {
			double val=map.get(key);
			double score=val/numPerm2;
			rtrn.put(key, score);
		}
		return rtrn;
	}

	
	
	
	
	
	
	private Cluster[] getPerms(Cluster c, int numPerm, IntervalTree<String> probabilityTree, int totalSize) {
		Cluster[] rtrn=new Cluster[numPerm];
		for(int i=0; i<numPerm; i++) {
			rtrn[i]=getPerm(c, probabilityTree, totalSize);
		}
		return rtrn;
	}

	private Cluster getPerm(Cluster c, IntervalTree<String> probabilityTree, int totalSize) {
		int size=c.getAllRNARegions().size();
		Cluster newCluster=new Cluster(c.getBarcode());
		Collection<RNAInterval> rnas=getRandom(probabilityTree, size, totalSize);
		newCluster.addRNAReads(rnas);
		return newCluster;
	}
	
	private Cluster getPerm(int size, IntervalTree<String> probabilityTree, int totalSize) {
		Cluster newCluster=new Cluster("rand");
		Collection<RNAInterval> rnas=getRandom(probabilityTree, size, totalSize);
		newCluster.addRNAReads(rnas);
		return newCluster;
	}
	
	private Collection<RNAInterval> getRandom(IntervalTree<String> tree, int size, int totalSize) {
		Collection<RNAInterval> rtrn=new TreeSet<RNAInterval>();
		for(int i=0; i<size; i++) {
			double random=Math.random();
			int pos=(int)(totalSize*random);
			String rna=tree.overlappingValueIterator(pos, pos+1).next();
			RNAInterval temp=new RNAInterval(new SingleInterval(rna, 0, 1));
			temp.setName(rna);
			rtrn.add(temp);
		}
		return rtrn;
	}
	
	
	private MatrixWithHeaders norm(MatrixWithHeaders observed, MatrixWithHeaders[] expected) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(observed.getRowNames(), observed.getColumnNames());
		
		List<String> rowList=new ArrayList<String>();
		for(String row: observed.getRowNames()) {
			double maxScore=Statistics.max(observed.getRow(row));
			if(maxScore>0) {rowList.add(row);}
			for(String column: observed.getColumnNames()) {
				double score=observed.get(row, column);
				double[] permScores=get(expected, row, column);
				double enrichment=0;
				if(score>0) {
					enrichment=score/Statistics.mean(permScores, score);
				}
				rtrn.set(row, column, enrichment);
			}
		}
		
		rtrn=rtrn.submatrixByRowNames(rowList);
		
		return rtrn;
	}

	private double[] get(MatrixWithHeaders[] expected, String row, String column) {
		double[] rtrn=new double[expected.length];
		
		for(int i=0; i<expected.length; i++) {
			rtrn[i]=expected[i].get(row, column);
		}
		
		return rtrn;
	}

	/*private void score(MatrixWithHeaders[] expected, Cluster[] permutedClusters, List<String> allRNA, int k) {
		for(int i=0; i<expected.length; i++) {
			score(expected[i], permutedClusters[i], allRNA, k);
		}
		
	}*/

	private void score(MatrixWithHeaders observed, Cluster c, Map<String, String> allRNA, int k, Map<Kmer, Collection<Cluster>> counts) {
		Kmer full=getHits(c, allRNA);
		Collection<Kmer> kHits=getKmers(full, k-1);
		populate(observed, kHits, full.getRegionList());
		
		Cluster renamed=c.renameRNA(allRNA);
		
		for(Kmer k1: kHits) {
			Collection<Cluster> list=new ArrayList<Cluster>();
			if(counts.containsKey(k1)) {list=counts.get(k1);}
			list.add(renamed);
			counts.put(k1, list);
		}
	}

	private MatrixWithHeaders[] makeMatrices(MatrixWithHeaders observed, int numPerm2) {
		MatrixWithHeaders[] rtrn=new MatrixWithHeaders[numPerm2];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=new MatrixWithHeaders(observed.getRowNames(), observed.getColumnNames());
		}
		
		return rtrn;
	}

	private Kmer getHits(Cluster c, Map<String, String> allRNA) {
		Kmer rtrn=new Kmer();
		
		for(String name: c.getRNANames()) {
			if(allRNA.containsKey(name)) {rtrn.addRegion(allRNA.get(name));}
		}
		
		return rtrn;
	}

	private Collection<Kmer> getKmers(Kmer full, int i) {
		return full.enumerateSubK(i);
	}

	private Collection<Kmer> getKmers(Map<String, String> allRNA, int i) {
		Kmer k=new Kmer();
		k.addRegions(allRNA.values());
		return k.enumerateSubK(i);
	}
	
	
	private Collection<Kmer> getKmers(Collection<Kmer> hubs, int i) {
		Collection<Kmer> rtrn=new ArrayList<Kmer>();
		for(Kmer hub: hubs) {
			System.err.println(hub.getName());
			rtrn.addAll(hub.enumerateSubK(i));
		}
		
		return rtrn;
	}

	private MatrixWithHeaders makeMatrix(Collection<Kmer> kmers, Map<String, String> allRNA) {
		List<String> rows=new ArrayList<String>();
		for(Kmer k: kmers) {
			rows.add(k.toString());
		}
		
		List<String> classes=getValues(allRNA);
		
		return new MatrixWithHeaders(rows, classes);
	}
	
	private List<String> getValues(Map<String, String> allRNA) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String rna: allRNA.keySet()) {
			String rnaClass=allRNA.get(rna);
			if(!rtrn.contains(rnaClass)) {rtrn.add(rnaClass);}
		}
		
		return rtrn;
	}


	private void populate(MatrixWithHeaders mwh, Collection<Kmer> kHits, List<String> regionList) {
		for(Kmer k: kHits) {
			String row=k.toString();
			for(String col: regionList) {
				if(!k.containsRegion(col) && mwh.containsRow(row)) {
					mwh.incrementCount(row, col);
				}
			}
		}
		
	}

	private MatrixWithHeaders makeMatrix(Collection<Cluster> list, Kmer k, List<String> allRNA) {
		List<String> columns=allRNA;
		List<String> rows=getRows(list);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(Cluster c: list) {
			String row=c.getBarcode();
			for(String rna: c.getRNANames()) {
				if(rtrn.containsColumn(rna)) {
					rtrn.set(row, rna, 1);
				}
			}
		}
		
		return rtrn;	
	}
	
	private List<String> getRows(Collection<Cluster> clustersInHub) {
		List<String> rtrn=new ArrayList<String>();
		
		for(Cluster c: clustersInHub) {rtrn.add(c.getBarcode());}
		
		return rtrn;
	}
	
	private Map<Kmer, Collection<Cluster>> getClusters(BarcodingDataStreaming data, Kmer hub) {
		Map<Kmer, Collection<Cluster>> rtrn=new TreeMap<Kmer, Collection<Cluster>>();
		int count=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Kmer k=getKmers(c, hub);
			if(k.getSize()>0) {add(rtrn, k, c);}
			count++;
			if(count%1000000==0) {System.err.println(count);}
		}
		data.close();
		
		this.totalSize=count;
		
		return rtrn;
	}
	
	
	
	private void add(Map<Kmer, Collection<Cluster>> rtrn, Kmer k, Cluster c) {
		if(!rtrn.containsKey(k)) {rtrn.put(k, new TreeSet<Cluster>());}
		rtrn.get(k).add(c);
	}
	
	private Kmer getKmers(Cluster c, Kmer hub) {
		Kmer k=numHits(c, hub);
		k.setName(hub.getName());
		return k;
	}
	
	private Kmer numHits(Cluster c, Kmer hub) {
		Kmer rtrn=new Kmer();
		for(String region: hub.getRegions()){
			if(c.getRNANames().contains(region)){rtrn.addRegion(region);}
		}
		return rtrn;
	}


	private MatrixWithHeaders reorderByHub(MatrixWithHeaders kmerCounts, Collection<Kmer> hubs) {
		List<String> newRows=new ArrayList<String>();
		
		Map<String, List<String>> rowsByHub=new TreeMap<String, List<String>>();
		
		for(String row: kmerCounts.getRowNames()) {
			String bestHub=maxHub(hubs, kmerCounts, row);
			if(bestHub!=null) {
				if(!rowsByHub.containsKey(bestHub)) {rowsByHub.put(bestHub, new ArrayList<String>());}
				rowsByHub.get(bestHub).add(row);
			}
		}
		
		
		for(String hub: rowsByHub.keySet()) {
			newRows.addAll(rowsByHub.get(hub));
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(newRows, kmerCounts.getColumnNames());
		for(String row: newRows) {
			for(String col: kmerCounts.getColumnNames()) {
				rtrn.set(row, col, kmerCounts.get(row, col));
			}
		}
		
		return rtrn;
	}


	private String maxHub(Collection<Kmer> hubs, MatrixWithHeaders kmerCounts, String row) {
		String bestHub=null;
		double maxscore=0;
		
		for(Kmer hub: hubs) {
			double score=score(kmerCounts, row, hub);
			if(score>maxscore) {
				maxscore=score;
				bestHub=hub.getName();
			}
		}
		
		if(maxscore<3) {bestHub=null;}
		return bestHub;
		
	}


	private double score(MatrixWithHeaders kmerCounts, String row, Kmer hub) {
		double sum=0;
		for(String col: hub.getRegions()) {
			sum+=kmerCounts.get(row, col);
		}
		return sum;
	}


	private MatrixWithHeaders filter(MatrixWithHeaders kmerCounts) {
		List<String> subCol=new ArrayList<String>();
		
		for(String col: kmerCounts.getRowNames()) {
			double sum=sum(kmerCounts.getRow(col));
			if(sum>3) {subCol.add(col);}
		}
		
		return kmerCounts.submatrixByRowNames(subCol);
	}


	private double sum(double[] column) {
		double sum=0;
		
		for(int i=0; i<column.length; i++) {sum+=column[i];}
		
		return sum;
	}


	private List<String> getColumns(Collection<Cluster> clustersInHub) {
		List<String> rtrn=new ArrayList<String>();
		
		for(Cluster c: clustersInHub) {rtrn.add(c.getBarcode());}
		
		return rtrn;
	}


	private Collection<Cluster> getClustersInHub(BarcodingDataStreaming data, List<String> rows) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		Kmer k=new Kmer();
		k.addRegions(rows);
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			if(c.containsRNA(k, true)) {rtrn.add(c);}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		data.close();
		
		return rtrn;
	}

	
	
	private static Map<String, String> parseRNAs(String string) throws IOException {
		List<String> lines= BEDFileIO.loadLines(string);
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String line: lines) {
			String rna=line.split("\t")[1].replaceAll("\"", "");
			System.err.println(line.split("\t")[1]+" "+rna);
			String className=line.split("\t")[2];
			rtrn.put(rna, className);
		}
		
		return rtrn;
	}
	
	private static Collection<Kmer> parseHubs(String string) throws IOException {
		List<String> lines= BEDFileIO.loadLines(string);
		
		Map<String, Collection<String>> temp=new TreeMap<String, Collection<String>>();
		for(String line: lines) {
			String hub=line.split("\t")[0];
			String className=line.split("\t")[2];
			Collection<String> set=new TreeSet<String>();
			if(temp.containsKey(hub)) {set=temp.get(hub);}
			set.add(className);
			temp.put(hub, set);
		}
		
		
		Collection<Kmer> rtrn=new ArrayList<Kmer>();
		for(String hub: temp.keySet()) {
			Kmer kmer=new Kmer();
			kmer.setName(hub);
			kmer.addRegions(temp.get(hub));
			rtrn.add(kmer);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Map<String, String> allRNAs=parseRNAs(args[1]);
			Collection<Kmer> hubs=parseHubs(args[1]);
			String save=args[2];
			int kmer=Integer.parseInt(args[3]);
			new KmerClusters(data, hubs, allRNAs, save, kmer);
			
		}
		else {System.err.println(usage);}
	}
	

	private static String usage=" args[0]=clusters \n args[1]=all rnas \n args[2]=save \n args[3]=kmer size";
	
}

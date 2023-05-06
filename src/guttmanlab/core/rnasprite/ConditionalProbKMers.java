package guttmanlab.core.rnasprite;

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

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class ConditionalProbKMers {
	//IntervalTree<String> probabilityTree;
	int totalSize;
	Map<String, Integer> probabilities;
	private double threshold=2;
	private int min=2;
	
	public ConditionalProbKMers(BarcodingDataStreaming data, Collection<Kmer> hubs, String save, int numPerm) throws IOException {
		/*this.probabilityTree=data.getProbabilityTree();
		this.totalSize=probabilityTree.max().getEnd();
		this.probabilities=getProbabilities();*/
		
		this.probabilities=new TreeMap<String, Integer>();
		
		FileWriter writer=new FileWriter(save);
		
		for(Kmer hub: hubs) {
			//TODO: Go through all and match any keeping track of sizes for each kmer
			Map<Kmer, Map<Integer, Integer>> sizePerKmer=getClusterSizes(data, hub);
			
			//Conditional probabilities
			Map<Kmer, Double> enrichment=computeConditionalProbability(sizePerKmer, writer, numPerm);
		}
		
		writer.close();
	}
	
	
	public ConditionalProbKMers(BarcodingDataStreaming data, Kmer hub, String save) throws IOException {
		this.probabilities=new TreeMap<String, Integer>();
		
		Map<Kmer, Map<Integer, Integer>> sizePerKmer=getClusterSizes(data, hub);
			
		//Conditional probabilities
		Map<Kmer, Double> enrichment=computeConditionalProbability(sizePerKmer, 100);
		Map<Kmer, Collection<Cluster>> clusters=getClusters(data, hub);
		write(enrichment, clusters, save);
	}

	private void write(Map<Kmer, Double> enrichment, Map<Kmer, Collection<Cluster>> clusters, String save) throws IOException {
		MatrixWithHeaders matrix=null;
		for(Kmer k: enrichment.keySet()) {
			double score=enrichment.get(k);
			if(score>threshold) {
				Collection<Cluster> list=clusters.get(k);
				if(list!=null && list.size()>min) {
					System.err.println(k+" "+score+" "+list.size());
					MatrixWithHeaders mwh=makeMatrix(list, k);
					matrix=add(mwh, matrix);
				}
			}
		}
		matrix.write(save);
	}
	
	private MatrixWithHeaders makeMatrix(Collection<Cluster> list, Kmer k) {
		List<String> columns=k.getRegionList();
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


	private MatrixWithHeaders add(MatrixWithHeaders mwh, MatrixWithHeaders matrix) {
		if(matrix==null) {return mwh;}
		
		List<String> columns=merge(mwh.getColumnNames(), matrix.getColumnNames());
		List<String> rows=merge(mwh.getRowNames(), matrix.getRowNames());
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(String row: mwh.getRowNames()) {
			for(String column: mwh.getColumnNames()) {
				rtrn.set(row, column, mwh.get(row, column));
			}
		}
		
		for(String row: matrix.getRowNames()) {
			for(String column: matrix.getColumnNames()) {
				rtrn.set(row, column, matrix.get(row, column));
			}
		}
		
		return rtrn;
	}


	private List<String> merge(List<String> l1, List<String> l2) {
		List<String> rtrn=new ArrayList<String>();
		rtrn.addAll(l1);
		for(String name: l2) {
			if(!rtrn.contains(name)) {rtrn.add(name);}
		}
		return rtrn;
	}


	private Map<Kmer, Collection<Cluster>> getClusters(BarcodingDataStreaming data, Kmer hub) {
		Map<Kmer, Collection<Cluster>> rtrn=new TreeMap<Kmer, Collection<Cluster>>();
		int count=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Kmer k=getKmers(c, hub);
			if(k.getSize()>0) {add2(rtrn, k, c);}
			count++;
			if(count%1000000==0) {System.err.println(count);}
		}
		data.close();
		
		this.totalSize=count;
		
		return rtrn;
	}
	
	
	
	private void add2(Map<Kmer, Collection<Cluster>> rtrn, Kmer k, Cluster c) {
		if(!rtrn.containsKey(k)) {rtrn.put(k, new ArrayList<Cluster>());}
		rtrn.get(k).add(c);
	}
	
	


	private Map<Kmer, Map<Integer, Integer>> getClusterSizes(BarcodingDataStreaming data, Kmer hub) {
		Map<Kmer, Map<Integer, Integer>> rtrn=new TreeMap<Kmer, Map<Integer, Integer>>();
		int count=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Kmer k=getKmers(c, hub);
			if(k.getSize()>0) {add(rtrn, k, c);}
			count++;
			addProb(c.getRNANames());
			if(count%1000000==0) {System.err.println(count);}
		}
		data.close();
		
		this.totalSize=count;
		//Merge into cumulative
		return makeCumulative(rtrn);
	}
	
	private void addProb(Collection<String> rnaNames) {
		for(String rna: rnaNames) {
			int count=0;
			if(this.probabilities.containsKey(rna)) {count=this.probabilities.get(rna);}
			count++;
			this.probabilities.put(rna,  count);
		}
	}


	private Map<Kmer, Map<Integer, Integer>> makeCumulative(Map<Kmer, Map<Integer, Integer>> observed) {
		Map<Kmer, Map<Integer, Integer>> rtrn=new TreeMap<Kmer, Map<Integer, Integer>>();
		
		for(Kmer k: observed.keySet()) {
			Collection<Kmer> subK=getSubK(k);
			add(rtrn, subK, observed.get(k));
		}
		return rtrn;
	}
	
	private void add(Map<Kmer, Map<Integer, Integer>> counts, Collection<Kmer> kmers, Map<Integer, Integer> observed) {
		for(Kmer kmer: kmers) {
			add(counts, kmer, observed);
		}
	}
	
	
	private void add(Map<Kmer, Map<Integer, Integer>> counts, Kmer overlap, Map<Integer, Integer> add) {
		Map<Integer, Integer> temp;
		if(counts.containsKey(overlap)) {temp=counts.get(overlap);}
		else {temp=new TreeMap<Integer, Integer>();}
		
		merge(temp, add);
		counts.put(overlap, temp);
		
	}
	
	
	private void merge(Map<Integer, Integer> temp, Map<Integer, Integer> add) {
		for(Integer size: add.keySet()) {
			int a1=0;
			if(temp.containsKey(size)) {a1=temp.get(size);}
			int sum=add.get(size)+a1;
			temp.put(size, sum);
		}
	}


	private void add(Map<Kmer, Map<Integer, Integer>> rtrn, Kmer k, Cluster c) {
		Map<Integer, Integer> temp;
		if(rtrn.containsKey(k)) {
			temp=rtrn.get(k);
		}
		else {temp=new TreeMap<Integer, Integer>();}
		
		int size=c.getClusterSize();
		int num=0;
		if(temp.containsKey(size)) {num=temp.get(size);}
		num++;
		
		temp.put(size, num);
		rtrn.put(k, temp);
		
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


	private Collection<Cluster> matchAny(BarcodingDataStreaming data, Kmer hub) {
		return data.getClusters(hub, false);
	}


	private Collection<Kmer> getSubK(Kmer k) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		//rtrn.add(k);
		
		for (int i=1; i<=k.getSize(); i++) {
			Collection<Kmer> subKs=k.enumerateSubK(i);
			//System.err.println(i+" "+k.getSize()+" "+subKs.size());
			rtrn.addAll(subKs);
		}
		return rtrn;
	}
	
	
	private Map<Kmer, Double> computeConditionalProbability(Map<Kmer, Map<Integer, Integer>> clusterSizesPerKmer, FileWriter writer, int numPerm) throws IOException {
		Map<Kmer, Double> rtrn=new TreeMap<Kmer, Double>();
		int counter=0;
		for(Kmer k: clusterSizesPerKmer.keySet()) {
			if(k.getSize()>1) {
				double observed=size(clusterSizesPerKmer.get(k));
				Map<Kmer, double[]> randomVals=computeConditionalProbability(k, clusterSizesPerKmer, numPerm);
				writer.write(k.toString()+"\t"+k.getName()+"\t"+k.getSize()+"\t"+observed);
				double expected=0;
				for(Kmer sub: randomVals.keySet()) {
					expected=Math.max(expected, Statistics.mean(randomVals.get(sub)));
					//writer.write("\t"+sub.toString()+"\t"+Statistics.mean(randomVals.get(sub)));
				}
				double enrichment=observed/expected;
				rtrn.put(k, enrichment);
				writer.write("\t"+expected+"\t"+enrichment);
				writer.write("\n");
			}
			counter++;
			//System.err.println(counter+" "+clusterSizesPerKmer.size()+" "+k);
		}
		return rtrn;
	}
	
	private Map<Kmer, Double> computeConditionalProbability(Map<Kmer, Map<Integer, Integer>> clusterSizesPerKmer, int numPerm) throws IOException {
		Map<Kmer, Double> rtrn=new TreeMap<Kmer, Double>();
		for(Kmer k: clusterSizesPerKmer.keySet()) {
			if(k.getSize()>1) {
				double observed=size(clusterSizesPerKmer.get(k));
				Map<Kmer, double[]> randomVals=computeConditionalProbability(k, clusterSizesPerKmer, numPerm);
				double expected=0;
				for(Kmer sub: randomVals.keySet()) {
					expected=Math.max(expected, Statistics.mean(randomVals.get(sub)));
				}
				double enrichment=observed/expected;
				rtrn.put(k, enrichment);
			}
		}
		return rtrn;
	}

	private int size(Map<Integer, Integer> map) {
		int sum=0; 
		
		for(Integer size: map.keySet()) {sum+=map.get(size);}
		
		return sum;
	}


	private Map<Kmer, double[]> computeConditionalProbability(Kmer k,  Map<Kmer, Map<Integer, Integer>> clusterSizesPerKmer, int numPerm) {
		Map<Kmer, double[]> rtrn=new TreeMap<Kmer, double[]>();
		
		for(String region: k.getRegions()) {
			Kmer sub=subK(k, region);
			if(clusterSizesPerKmer.containsKey(sub)) {
				Map<Integer, Integer> clusterSizes=clusterSizesPerKmer.get(sub);
				//fix subK and randomize the rest
				double[] random=scoreRandom(clusterSizes, sub, region, numPerm);
				rtrn.put(sub, random);
			}
			else {System.err.println("skipping "+sub);}
		}
		
		return rtrn;
	}
	
	private double[] scoreRandom(Map<Integer, Integer> clusterSizes, Kmer sub, String region, int numPerm) {
		double[] rtrn=new double[numPerm];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=scoreRandom(clusterSizes, sub, region);
		}
		
		return rtrn;
	}
	
	private double scoreRandom(Map<Integer, Integer> clusterSizes, Kmer sub, String region) {
		double sum=0;
		for(Integer clusterSize: clusterSizes.keySet()) {
			int num=clusterSizes.get(clusterSize);
			double prob=probability(clusterSize, sub, region);
			sum+=(num*prob);
			/*for(int i=0; i<num; i++) {
				sum+=randomHasHit(clusterSize, sub, region); //TODO Bernoulli dist?
			}*/
		}
		return sum;
	}


	private double probability(int clusterSize, Kmer sub, String region) {
		double p=getP(region);
		int numTrials=clusterSize-sub.getSize();
		BinomialDistribution dist=new BinomialDistribution(numTrials, p);
		
		double rtrn=1-dist.probability(0);
		//System.err.println(region+" "+p+" "+numTrials+" "+rtrn);
		return rtrn;
	}
	
	private double getP(String region) {
		return (double)this.probabilities.get(region)/(double)totalSize;
	}
	
	/*private int randomHasHit(int clusterSize, Kmer sub, String region) {
		int size=clusterSize-sub.getSize();
		boolean hasRNA=getRandom(probabilityTree, size, totalSize, region);
		if(hasRNA) {return 1;}
		return 0;
	}*/


	private boolean getRandom(IntervalTree<String> tree, int size, int totalSize, String region) {
		for(int i=0; i<size; i++) {
			double random=Math.random();
			int pos=(int)(totalSize*random);
			String rna=tree.overlappingValueIterator(pos, pos+1).next();
			if(rna.equals(region)) {return true;}
		}
		return false;
	}
	
	
	private Collection<Cluster> getClusters(Collection<Cluster> data, Kmer k) {
		Collection<Cluster> rtrn=new HashSet<Cluster>();
		
		for(Cluster c: data) {
			if(c.containsRNA(k, true)) {rtrn.add(c);}
		}
		
		return rtrn;
	}

	private Kmer subK(Kmer k, String region) {
		Kmer rtrn=new Kmer();
		rtrn.setName(k.getName());
		for(String r: k.getRegions()) {
			if(!r.equals(region)) {rtrn.addRegion(r);}
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
		if(args.length>2) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Collection<Kmer> hubs=parseHubs(args[1]);
			String save=args[2];
			new ConditionalProbKMers(data, hubs.iterator().next(), save);
			
		}
		else {System.err.println(usage);}
	}
	
	



	private static String usage=" args[0]=clusters \n args[1]=RNA hub file \n args[2]=save";
	
}

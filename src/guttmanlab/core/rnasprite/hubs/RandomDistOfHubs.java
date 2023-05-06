package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
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
import guttmanlab.core.math.Statistics;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.rnasprite.Kmer;
import guttmanlab.core.rnasprite.RNAInterval;
import guttmanlab.core.sars.batch.ProbabilityTree;

public class RandomDistOfHubs {
	int totalSize;
	int numPerm=100;
	IntervalTree<String> probabilityTree;
	
	
	public RandomDistOfHubs(BarcodingDataStreaming data, IntervalTree<String> probabilityTree, Kmer hub, Map<String, String> allRNA, String save) throws IOException {
		this.probabilityTree=probabilityTree;
		this.totalSize=probabilityTree.max().getEnd();
		
		Map<String, Double> probabilites=makeProbabilities(probabilityTree);
		
		Collection<Cluster> clusters=getClusters(data, hub, allRNA);
		
		System.err.println(clusters.size());
		
		double observed=getClusters(clusters, hub, allRNA).size();
		Map<String, double[]> perms=scoreConditional(clusters, hub, allRNA, probabilites);
		write(save, perms, observed, hub);
	}
	

	private Collection<Cluster> getClusters(BarcodingDataStreaming data, Kmer hub, Map<String, String> allRNA) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			int count=count(c, hub, allRNA);
			if(count>=hub.getSize()) {rtrn.add(c);}
			counter++;
			if(counter%10000==0) {System.err.println(counter);}
		}
		
		data.close();
		
		return rtrn;
	}


	private int count(Cluster c, Kmer hub, Map<String, String> allRNA) {
		Collection<String> set=new TreeSet<String>();
		
		for(String r: c.getRNANames()) {
			if(allRNA.containsKey(r)) {
				String name=allRNA.get(r);
				if(hub.containsRegion(name)) {set.add(name);}
			}
		}
		return set.size();
	}


	private void write(String save, Map<String, double[]> perms, double observed, Kmer hub) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String r: perms.keySet()) {
			double[] vals=perms.get(r);
			double avg=Statistics.mean(vals);
			double enrichment=observed/avg;
			double z=Statistics.zScore(observed, vals);
			double p=Statistics.percentLessThan(observed, vals);
			writer.write(hub.toString()+"\t"+r+"\t"+observed+"\t"+avg+"\t"+enrichment+"\t"+z+"\t"+p);
			for(int i=0; i<vals.length; i++) {
				writer.write("\t"+vals[i]);
			}
			writer.write("\n");
		}
		
		writer.close();
	}


	private Map<String, Double> makeProbabilities(IntervalTree<String> probabilityTree2) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		Iterator<Node<String>> iter=probabilityTree2.iterator();
		
		while(iter.hasNext()) {
			Node<String> next=iter.next();
			double length=next.getEnd()-next.getStart();
			double ratio=length/this.totalSize;
			rtrn.put(next.getValue(), ratio);
		}
		
		return rtrn;
	}


	private Map<String, double[]> scoreConditional(BarcodingDataStreaming data, Kmer hub, Map<String, String> allRNA, Map<String, Double> probabilites) {
		Map<String, double[]> rtrn=new TreeMap<String, double[]>();
		
		
		
		
		for(String region: hub.getRegions()) {
			System.err.println(region);
			Kmer sub=hub.remove(region);
			Collection<Cluster> clusters=getClusters(data, sub, allRNA);
			double[] permVals=new double[numPerm];
			//int i=0;
			for(Cluster c: clusters) {
				double[] vals=perm(c, numPerm, region, sub.getSize(), probabilites);
				permVals=add(permVals, vals);
				//i++;
				//System.err.println(i+" "+clusters.size());
			}
			rtrn.put(region, permVals);
		}
		return rtrn;
	}
	
	
	private Map<String, double[]> scoreConditional(Collection<Cluster> data, Kmer hub, Map<String, String> allRNA, Map<String, Double> probabilites) {
		Map<String, double[]> rtrn=new TreeMap<String, double[]>();
		
		for(String region: hub.getRegions()) {
			System.err.println(region);
			Kmer sub=hub.remove(region);
			Collection<Cluster> clusters=getClusters(data, sub, allRNA);
			double[] permVals=new double[numPerm];
			for(Cluster c: clusters) {
				double[] vals=perm(c, numPerm, region, sub.getSize(), probabilites);
				permVals=add(permVals, vals);
			}
			rtrn.put(region, permVals);
		}
		return rtrn;
	}
	
	private Collection<Cluster> getClusters(Collection<Cluster> data, Kmer sub, Map<String, String> allRNA) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		//while(data.hasNext()) {
		for(Cluster c: data) {	
			//Cluster c=data.next();
			if(has(c, sub, allRNA)) {rtrn.add(c);}
			counter++;
			if(counter%10000==0) {System.err.println(counter);}
		}
		
		//data.close();
		
		return rtrn;
	}


	private boolean has(Cluster c, Kmer sub, Map<String, String> allRNA) {
		Kmer hits=getHits(c, sub, allRNA);
		
		for(String region: sub.getRegions()) {
			if(!hits.containsRegion(region)) {return false;}
		}
		return true;
	}


	private double[] add(double[] permVals, double[] vals) {
		double[] rtrn=new double[permVals.length];
		for(int i=0; i<permVals.length; i++) {rtrn[i]=permVals[i]+vals[i];}
		return rtrn;
	}


	private double[] perm(Cluster c, int numPerm2, String region, int subSize, Map<String, Double> probabilites) {
		double[] rtrn=new double[numPerm2];
		
		double p=getP(region, probabilites);
		int numTrials=c.getClusterSize()-subSize;
		BinomialDistribution dist=new BinomialDistribution(numTrials, p);
		
		double cutoff=1-dist.probability(0);
		
		for(int i=0; i<numPerm2; i++) {
			if(Math.random()<cutoff) {rtrn[i]=1;}
		}
		
		return rtrn;
	}

	
	private double getP(String region, Map<String, Double> probabilites) {
		return probabilites.get(region);
	}

	


	/*private void scoreConditional(Kmer hub, Map<Kmer, Collection<Cluster>> counts, String save) throws IOException {
		
		FileWriter writer=new FileWriter(save);
		
		double observed=0;
		if(counts.containsKey(hub)) {
			observed=counts.get(hub).size();
		}
		Map<Kmer, double[]> subKScores=new TreeMap<Kmer, double[]>();
		int count=0;
		for(String remaining: hub.getRegions()) {
			System.err.println(count+" "+hub.getSize());	
			Kmer subK=hub.remove(remaining);
			Collection<Cluster> clusters=counts.get(subK);
			Map<String, double[]> perm=perm(clusters, subK, numPerm);
			double[] permScoresByRegion=new double[numPerm];
			if(perm.containsKey(remaining)) {
				permScoresByRegion=perm.get(remaining);
			}
			subKScores.put(subK, permScoresByRegion);
			write(writer, subK, remaining, observed, permScoresByRegion);
			count++;
		}
		
		writer.close();
	}*/

	private void write(FileWriter writer, Kmer subK, String remaining, double observed, double[] permScoresByRegion) throws IOException {
		double expected=Statistics.mean(permScoresByRegion);
		double enrichment=observed/expected;
		double p=Statistics.percentLessThan(observed, permScoresByRegion);
		double z=Statistics.zScore(observed, permScoresByRegion);
		writer.write(subK+"\t"+remaining+"\t"+observed+"\t"+expected+"\t"+enrichment+"\t"+p+"\t"+z);
		
		for(int i=0; i<permScoresByRegion.length; i++) {writer.write("\t"+permScoresByRegion[i]);}
		
		writer.write("\n");
	}


	private Map<String, double[]> perm(Collection<Cluster> clusters, Kmer k, int numPerm2) {
		Map<String, double[]> rtrn=new TreeMap<String, double[]>();
		
		if(clusters==null || clusters.isEmpty()) {return rtrn;}
		
		
		
		for(int i=0; i<numPerm2; i++) {
			Map<String, Double> regions=perm(clusters, k);
			update(regions, rtrn, i, numPerm2);
		}
		
		
		return rtrn;
	}

	

	private Map<String, Double> perm(Collection<Cluster> clusters, Kmer k) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		for(Cluster c: clusters) {
			int size=c.size()-k.getRegions().size();
			Cluster rand=this.getPerm(size, probabilityTree, totalSize);
			update(rand, rtrn);
		}
		return rtrn;
	}

	private void update(Cluster rand, Map<String, Double> rtrn) {
		for(String region: rand.getRNANames()) {
			double count=0;
			if(rtrn.containsKey(region)) {count=rtrn.get(region);}
			count++;
			rtrn.put(region, count);
		}
		
	}

	private void update(Map<String, Double> regions, Map<String, double[]> rtrn, int index, int numPerm) {
		for(String region: regions.keySet()) {
			if(!rtrn.containsKey(region)) {rtrn.put(region, new double[numPerm]);}
			double[] array=rtrn.get(region);
			array[index]=regions.get(region);
			rtrn.put(region, array);
		}
		
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
	

	



	private boolean match(Kmer full, Kmer k) {
		for(String region: k.getRegions()) {
			if(!full.containsRegion(region)) {return false;}
		}
		return true;
	}


	

	private Kmer getHits(Cluster c, Kmer k, Map<String, String> allRNA) {
		Kmer rtrn=new Kmer();
		
		for(String name: c.getRNANames()) {
			if(allRNA.containsKey(name)) {
				if(k.containsRegion(allRNA.get(name))) {
					rtrn.addRegion(allRNA.get(name));
				}
			}
		}
		
		return rtrn;
	}

	private Collection<Kmer> getKmers(Kmer full, int i) {
		return full.enumerateSubK(i);
	}

	
	
	private static Map<String, String> parseRNAs(String string) throws IOException {
		List<String> lines= BEDFileIO.loadLines(string);
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String line: lines) {
			String rna=line.split("\t")[1].replaceAll("\"", "");
			String className=line.split("\t")[2];
			rtrn.put(rna, className);
		}
		
		return rtrn;
	}
	
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>4) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			IntervalTree<String> probTree=ProbabilityTree.parseProbabilityTree(args[1]);
			
			
			Map<String, String> allRNAs=parseRNAs(args[2]);
			Kmer hub=new Kmer(args[3]);
			String save=args[4];
			
			new RandomDistOfHubs(data, probTree, hub, allRNAs, save);
			System.err.println("done");
			
		}
		else {System.err.println(usage);}
	}
	

	private static String usage=" args[0]=clusters \n args[1]=prob tree \n args[2]=all rnas \n args[3]=kmer \n args[4]=save";
	
}

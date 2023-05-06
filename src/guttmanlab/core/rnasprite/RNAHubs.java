package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.Statistics;

public class RNAHubs {

	
	public RNAHubs(BarcodingDataStreaming data, Collection<Kmer> hubs, String save, int numPerm) throws IOException {
		Map<Kmer, Integer> observed=new TreeMap<Kmer, Integer>();
		Map<Kmer, Integer>[] perm=initialize(numPerm);
		
		IntervalTree<String> probabilityTree=data.getProbabilityTree();
		int totalSize=probabilityTree.max().getEnd();
		System.err.println(totalSize);
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			getKmers(c, hubs, observed);
			
			Cluster[] permutedClusters=getPerms(c, numPerm, probabilityTree, totalSize);
			getKmers(permutedClusters, hubs, perm);
			
			counter++;
			if(counter%10000 ==0) {System.err.println("perms "+counter);}
		}
		System.err.println(counter);
		data.close();
		write(save, observed, perm);
		
		
		//make cumulative kmers
		Map<Kmer, Integer> observedCDF=makeCumulative(observed);
		Map<Kmer, Integer>[] permCDF=makeCumulative(perm);
		
		write(save+".cumulative", observedCDF, permCDF);
		
		
		
	}
	
	

	/*public RNAHubs(BarcodingDataStreaming data, Collection<Kmer> hubs, String save, int numPerm) throws IOException {
		Kmer allRNAs=getRNAs(hubs);
		
		Map<Kmer, Integer> observed=new TreeMap<Kmer, Integer>();
		Map<Kmer, Integer>[] perm=initialize(numPerm);
		
		IntervalTree<String> probabilityTree=data.getProbabilityTree();
		int totalSize=probabilityTree.max().getEnd();
		System.err.println(totalSize);
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			getKmers(c, allRNAs, observed);
			
			Cluster[] permutedClusters=getPerms(c, numPerm, probabilityTree, totalSize);
			getKmers(permutedClusters, allRNAs, perm);
			
			counter++;
			if(counter%10000 ==0) {System.err.println("perms "+counter);}
		}
		System.err.println(counter);
		data.close();
		write(save, observed, perm);
		
		
		//make cumulative kmers
		Map<Kmer, Integer> observedCDF=makeCumulative(observed);
		Map<Kmer, Integer>[] permCDF=makeCumulative(perm);
		
		write(save+".cumulative", observedCDF, permCDF);
		
	}*/
	
	
	private Kmer getRNAs(Collection<Kmer> hubs) {
		Kmer rtrn=new Kmer();
		rtrn.setName("all");
		
		for(Kmer k: hubs) {rtrn.addRegions(k.getRegions());}
		
		return rtrn;
	}


	private Map<Kmer, Integer>[] makeCumulative(Map<Kmer, Integer>[] perm) {
		Map<Kmer, Integer>[] rtrn=new Map[perm.length];
		
		for(int i=0; i<perm.length; i++) {
			rtrn[i]=makeCumulative(perm[i]);
		}
		
		return rtrn;
	}

	private Map<Kmer, Integer> makeCumulative(Map<Kmer, Integer> observed) {
		Map<Kmer, Integer> rtrn=new TreeMap<Kmer, Integer>();
		
		for(Kmer k: observed.keySet()) {
			Collection<Kmer> subK=getSubK(k);
			add(rtrn, subK, observed.get(k));
		}
		return rtrn;
	}

	private Collection<Kmer> getSubK(Kmer k) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		//rtrn.add(k);
		
		for (int i=1; i<=k.getSize(); i++) {
			Collection<Kmer> subKs=k.enumerateSubK(i);
			rtrn.addAll(subKs);
		}
		return rtrn;
	}

	private Map<Kmer, Integer>[] initialize(int numPerm) {
		Map<Kmer, Integer>[] rtrn=new Map[numPerm];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=new TreeMap<Kmer, Integer>();
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


	private void getKmers(Cluster[] clusters, Collection<Kmer> hubs, Map<Kmer, Integer>[] expected) {
		for(int i=0; i<clusters.length; i++) {
			getKmers(clusters[i], hubs, expected[i]);
		}
	}
	
	private void getKmers(Cluster[] clusters, Kmer hub, Map<Kmer, Integer>[] expected) {
		for(int i=0; i<clusters.length; i++) {
			getKmers(clusters[i], hub, expected[i]);
		}
	}
	
	
	private void getKmers(Cluster c, Kmer hub, Map<Kmer, Integer> observed) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		Kmer k=numHits(c, hub);
		if(k.getSize()>0) {rtrn.add(k);}
		add(observed, rtrn);
	}
	
	private void getKmers(Cluster c, Collection<Kmer> hubs, Map<Kmer, Integer> observed) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		for(Kmer hub: hubs) {
			Kmer k=numHits(c, hub);
			k.setName(hub.getName());
			if(k.getSize()>0) {rtrn.add(k);}
		}
		
		add(observed, rtrn);
	}
	
	private Collection<Kmer> getKmers(Cluster c, Collection<Kmer> hubs) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		for(Kmer hub: hubs) {
			Kmer k=numHits(c, hub);
			if(k.getSize()>0) {rtrn.add(k);}
		}
		
		return rtrn;
	}
	
	private void write(String save, Map<Kmer, Integer> observed, Map<Kmer, Integer>[] expected) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Kmer> allKmers=new TreeSet<Kmer>();
		allKmers.addAll(observed.keySet());
		for(int i=0; i<expected.length; i++) {allKmers.addAll(expected[i].keySet());}
		
		for(Kmer k: allKmers) {
			int observedScore=get(observed, k);
			int[] expectedScores=get(expected, k);
			writer.write(k.toString()+"\t"+k.getName()+"\t"+k.getSize()+"\t"+observedScore);
			for(int i=0; i<expectedScores.length; i++) {writer.write("\t"+expectedScores[i]);}
			writer.write("\n");
		}
		
		writer.close();
	}
	

	private int[] get(Map<Kmer, Integer>[] expected, Kmer k) {
		int[] scores=new int[expected.length];
		
		for(int i=0; i<expected.length; i++) {
			scores[i]=get(expected[i],k);
		}
		
		return scores;
	}


	private int get(Map<Kmer, Integer> observed, Kmer k) {
		int score=0;
		if(observed.containsKey(k)) {score=observed.get(k);}
		return score;
	}


	private void write(String save, Map<Kmer, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(Kmer k: counts.keySet()) {
			writer.write(k.toString()+"\t"+k.getSize()+"\t"+counts.get(k)+"\n");
		}
		writer.close();
	}
	
	private void add(Map<Kmer, Integer> counts, Collection<Kmer> kmers, int observed) {
		for(Kmer kmer: kmers) {
			add(counts, kmer, observed);
		}
	}
	
	private void add(Map<Kmer, Integer> counts, Collection<Kmer> kmers) {
		for(Kmer kmer: kmers) {
			add(counts, kmer);
		}
	}
	
	private void add(Map<Kmer, Integer> counts, Kmer overlap, int observed) {
		int add=observed;
		
		int count=0;
		if(counts.containsKey(overlap)) {count=counts.get(overlap);}
		count+=add;
		counts.put(overlap, count);
		
	}

	private void add(Map<Kmer, Integer> counts, Kmer overlap) {
		int count=0;
		if(counts.containsKey(overlap)) {count=counts.get(overlap);}
		count++;
		counts.put(overlap, count);
		
	}

	private Kmer numHits(Cluster c, Kmer hub) {
		Kmer rtrn=new Kmer();
		for(String region: hub.getRegions()){
			if(c.getRNANames().contains(region)){rtrn.addRegion(region);}
		}
		return rtrn;
	}

	private Kmer getAllRNAs(Collection<Kmer> kmers) {
		Kmer rtrn=new Kmer();
		for(Kmer kmer: kmers) {rtrn.addRegions(kmer.getRegions());}
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
			new RNAHubs(data, hubs, save, 100);
			
		}
		else {System.err.println(usage);}
	}
	
	



	private static String usage=" args[0]=clusters \n args[1]=RNA hub file \n args[2]=save";
	
}

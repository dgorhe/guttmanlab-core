package guttmanlab.core.nanobody;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import guttmanlab.core.sequence.Sequence;

public class ClusterProteinsBySequence {

	static int length=15;
	static int minDistance=2;
	static int kmerSize=5;
	static String badKmer="QVTVSS";
	
	public ClusterProteinsBySequence(File file, String save) throws IOException {
		/*Collection<String> uniqueSequences=parse(files);
		System.err.println(uniqueSequences.size());*/
		
		Map<String, Integer> kmerCounts=kmerCount(file, kmerSize);
		
		write(save, kmerCounts);
		
		//TODO cluster
		/*Collection<Cluster> clusters=enumerateKmers(uniqueSequences, kmerSize, minDistance);
		
		write(clusters, save);
		
		Collection<Cluster> filtered=collapse(clusters);
		write(filtered, save+".filtered");*/
		
		/*clusters=collapse(clusters, input); //TODO score
		write(clusters, save+".collapsed");*/
		
	}
	
	

	

	

	private Map<String, Integer> kmerCount(File file, int kmerSize2) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		Map<String, Integer> list=parse(file);
		System.err.println(list.size());
		
		int counter=0;
		
		for(String seq: list.keySet()) {
			
			Collection<String> kmers=enumerateKmers(seq, kmerSize2);
			for(String kmer: kmers) {
				int count=0;
				if(rtrn.containsKey(kmer)) {count=rtrn.get(kmer);}
				count+=list.get(seq);
				rtrn.put(kmer, count);
			}
				
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		
		return rtrn;
	}

	private Collection<String> enumerateKmers(String seq, int kmerSize2) {
		Collection<String> rtrn=new TreeSet<String>();
		Map<Integer, String> kmers=Sequence.enumerateKmer(seq, kmerSize2);
		for(Integer pos: kmers.keySet()) {
			String sub=kmers.get(pos);
			rtrn.add(sub);
			//rtrn.add(pos+"_"+sub);
			//System.err.println(pos+" "+sub+" "+seq);
		}
		return rtrn;
	}

	private Collection<Cluster> collapse(Collection<Cluster> clusters) {
		Collection<String> subsets=new TreeSet<String>();
		
		for(Cluster c: clusters) {
			if(!subsets.contains(c.getName())) {
				Collection<String> subs=getSubset(c, clusters);
				subsets.addAll(subs);
				System.err.println(c.getName()+" "+subs.size()+" "+subs);
			}
		}
		
		
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Cluster c: clusters) {
			if(!subsets.contains(c.getName())) {rtrn.add(c);}
		}
		
		return rtrn;
	}

	private Collection<String> getSubset(Cluster c, Collection<Cluster> clusters) {
		Collection<String> rtrn=new TreeSet<String>();
		for(Cluster other: clusters) {
			if(other.isSubset(c)) {rtrn.add(other.getName());}
		}
		return rtrn;
	}

	/*private Collection<String> parse(File[] files) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn.addAll(parse(files[i]));
		}
		return rtrn;
	}*/

	private Map<String, Integer> parse(File file) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String seq;
			if(tokens.length>1) {seq=tokens[3];}
			else {seq=tokens[0];}
			if(!seq.contains(badKmer) && !seq.contains("*")) {
				int count=0;
				if(rtrn.containsKey(seq)) {count=rtrn.get(seq);}
				count++;
				rtrn.put(seq, count);
			}
			if(counter%1000000==0) {System.err.println(counter);}
			counter++;
		}
		reader.close();
		return rtrn;
	}

	private Collection<Cluster> collapse(Collection<Cluster> clusters, String input) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		Map<String, Cluster> map=makeMap(clusters);
		Map<String, Integer> counts=getCounts(input, map.keySet());
		
		Collection<String> visitedProteins=new TreeSet<String>();
		for(String protein: map.keySet()) {
			if(!visitedProteins.contains(protein)) {
				Cluster mergedClster=map.get(protein);
				int sum=sum(mergedClster, counts);
				mergedClster.setScore(sum);
				rtrn.add(mergedClster);
				visitedProteins.addAll(mergedClster.proteins);
			}
		}
		
		/*for(Cluster c: clusters) {
			Cluster mergedClster=getAllClusters(c, map);
			int sum=sum(mergedClster, counts);
			mergedClster.setScore(sum);
			rtrn.add(mergedClster);
		}*/
		
		/*
		Collection<Cluster> multiples=getSingles(input, map.keySet()); //TODO Retain singles with multiple counts
		for(Cluster c: multiples) {
			//System.err.println(c);
			rtrn.add(c);
		}*/
		
		
		return rtrn;
	}
	
	
	private Collection<Cluster> getSingles(String input, Set<String> keySet) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String line;
		int counter=0;
		while ((line = reader.readLine()) != null) {
			String seq=line.split("\t")[0];
			if(!seq.contains("*") && !keySet.contains(seq)) {
				int count=0;
				if(rtrn.containsKey(seq)) {
					count=rtrn.get(seq);
				}
				count++;
				rtrn.put(seq, count);
				if(counter%1000000==0) {System.err.println(counter);}
			}
			counter++;
		}
		
		reader.close();
		
		Collection<Cluster> clusters=new TreeSet<Cluster>();
		clusters=filterMultiple(rtrn);
		
		return clusters;
	}
	
	private Collection<Cluster> filterMultiple(Map<String, Integer> map) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(String s: map.keySet()) {
			int count=map.get(s);
			if(count>1) {
				Collection<String> list=new TreeSet<String>();
				list.add(s);
				Cluster c=new Cluster(list);
				c.setScore(count);
				rtrn.add(c);
			}
		}
		return rtrn;
	}

	private Map<String, Integer> getCounts(String input, Set<String> keySet) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String line;
		int counter=0;
		while ((line = reader.readLine()) != null) {
			String seq=line.split(" ")[1];
			if(!seq.contains("X") && keySet.contains(seq)) {
				int count=Integer.parseInt(line.split(" ")[0]);
				//System.err.println(seq+" "+count);
				if(rtrn.containsKey(seq)) {
					count+=rtrn.get(seq);
				}
				//count++;
				rtrn.put(seq, count);
				if(counter%1000000==0) {System.err.println(counter);}
			}
			counter++;
		}
		
		reader.close();
		return rtrn;
	}

	private Map<String, Cluster> makeMap(Collection<Cluster> clusters) {
		Map<String, Cluster> rtrn=new TreeMap<String, Cluster>();
		
		for(Cluster c: clusters) {
			for(String protein: c.proteins) {
				Collection<String> list=new TreeSet<String>();
				list.addAll(c.proteins);
				if(rtrn.containsKey(protein)) {
					list.addAll(rtrn.get(protein).proteins);
				}
				rtrn.put(protein, new Cluster(list));	
			}
		}
		
		Map<String, Cluster> finalRtrn=new TreeMap<String, Cluster>();
		for(String protein: rtrn.keySet()) {
			Collection<String> clusterP=rtrn.get(protein).proteins;
			for(String newPro: clusterP) {
				clusterP.addAll(rtrn.get(newPro).proteins);
			}
			Cluster newCluster=new Cluster(clusterP);
			finalRtrn.put(protein, newCluster);
		}
		
		return finalRtrn;
	}

	private Cluster getAllClusters(Cluster c, Map<String, Cluster> map) {
		Collection<String> list=new TreeSet<String>();
		
		for(String p: c.proteins) {
			Cluster c1=map.get(p);
			list.addAll(c1.proteins);
		}
		
		
		return new Cluster(list);
	}

	/*private void getMultiples(Collection<String> kmerFiles) throws IOException {
		Map<String, Integer> allKmers=new TreeMap<String, Integer>();
		
		
		for(String file: kmerFiles) {
			
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
			String line;
			int counter=0;
			while ((line = reader.readLine()) != null) {
				String[] tokens=line.split("\t");
				String seq=tokens[0];
				int size=tokens.length-1;
				int count=size;
				if(allKmers.containsKey(seq)) {
					count=allKmers.get(seq)+size;
				}
				allKmers.put(seq, count);
				
				counter++;
			}
			
			reader.close();
			
		}
		
		
	}*/


	private static void write(String save, Map<String, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String k: counts.keySet()) {
			int count=counts.get(k);
			writer.write(k+"\t"+count+"\n");
		}
		
		writer.close();
	}

	
	private static Collection<Cluster> enumerateKmers(Collection<String> sequences, int k, int minDistance) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(int i=0; i<=length-k; i++) {
			System.err.println("position "+i);
			Map<String, Collection<String>> kmers=enumerate(sequences, k, i);
			
			Collection<Cluster> related=new TreeSet<Cluster>();
			for(String kmer: kmers.keySet()) {
				Collection<String> list=kmers.get(kmer);
				if(list.size()>1) {
					if(list.size()>100) {System.err.println(kmer+" "+list.size());}
					Collection<Cluster> map=getClose(list, minDistance);
					related.addAll(map);
				}
			}
			rtrn.addAll(related);
			//write(related, save+"."+i+".clusters");
		}
		return rtrn;
	}
	
	private static Map<String, Collection<String>> enumerate(Collection<String> sequences, int k, int pos) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		int counter=0;
		for(String seq: sequences) {
			String kmer=seq.substring(pos, pos+k);
			Collection<String> list=new TreeSet<String>();
			if(rtrn.containsKey(kmer)) {
				list=rtrn.get(kmer);
			}
			list.add(seq);
			rtrn.put(kmer, list);
			if(counter%1000000==0) {System.err.println(counter+" "+seq+" "+kmer+" "+kmer.length());}
			counter++;
		}
		
		return rtrn;
	}

	private static Collection<Cluster> enumerateKmers(String fileName, int k, String save, int minDistance) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(int i=0; i<=length-k; i++) {
			System.err.println("position "+i);
			Map<String, Collection<String>> kmers=enumerateKmers(fileName, k, i);
			
			Collection<Cluster> related=new TreeSet<Cluster>();
			for(String kmer: kmers.keySet()) {
				Collection<String> list=kmers.get(kmer);
				if(list.size()>1) {
					Collection<Cluster> map=getClose(list, minDistance);
					related.addAll(map);
				}
			}
			rtrn.addAll(related);
			//write(related, save+"."+i+".clusters");
		}
		return rtrn;
	}
	
	private static Map<String, Collection<String>> enumerateKmers(String fileName, int k, int pos) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String line;
		int counter=0;
		while ((line = reader.readLine()) != null) {
			String seq=line.split(" ")[1];
			if(!seq.contains("X") && !seq.contains(badKmer)) {
				String kmer=seq.substring(pos, pos+k);
				Collection<String> list=new TreeSet<String>();
				//int count=0;
				if(rtrn.containsKey(kmer)) {
					list=rtrn.get(kmer);
				}
				list.add(seq);
				rtrn.put(kmer, list);
				if(counter%1000000==0) {System.err.println(counter+" "+seq+" "+kmer+" "+kmer.length());}
			}
			counter++;
		}
		
		reader.close();
		return rtrn;
	}


	private static void write(Collection<Cluster> clusters, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		for(Cluster cluster: clusters) {
			String name="Cluster_"+counter;
			cluster.setName(name);
			writer.write(name+"\t"+cluster+"\n");
			counter++;
		}
		
		writer.close();
	}


	/*private static Collection<String> enumerateKmers(String fileName, int k, int maxFileSize, String save) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		Map<String, Collection<String>> kmersToFull=new TreeMap<String, Collection<String>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String line;
		int counter=0;
		while ((line = reader.readLine()) != null) {
			String seq=line.split("\t")[0];
			if(!seq.contains("*")) {
				Map<Integer, String> kmers=Sequence.enumerateKmer(seq, k);
				for(Integer pos: kmers.keySet()) {
					String kmer=kmers.get(pos);
					Collection<String> list=new TreeSet<String>();
					if(kmersToFull.containsKey(kmer)) {
						list=kmersToFull.get(kmer);
					}
					list.add(seq);
					kmersToFull.put(kmer, list);
				}
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
			if(counter%maxFileSize==0) {
				System.err.println(counter+" reset");
				String out=save+"."+(counter/maxFileSize)+".kmers";
				write(kmersToFull, out);
				rtrn.add(out);
				kmersToFull=new TreeMap<String, Collection<String>>();
			}
		}
		
		reader.close();
		
		return rtrn;
	}*/

	private int getMinScore(Map<String, Integer> distanceMap, int max) {
		int min=max;
		
		for(String s: distanceMap.keySet()) {
			int score=distanceMap.get(s);
			min=Math.min(min, score);
		}
		
		return min;
	}

	private void write(String save, Sequence seq1, Map<String, Integer> distanceMap) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String s: distanceMap.keySet()) {
			writer.write(seq1.getSequenceBases()+"\t"+s+"\t"+distanceMap.get(s)+"\n");
		}
		
		writer.close();
	}

	private static int distance(Sequence seq1, Sequence seq2) {
		char[] char1=seq1.getSequenceBases().toCharArray();
		char[] char2=seq2.getSequenceBases().toCharArray();
		
		int distance=0;
		for(int i=0; i<char1.length; i++) {
			if(char1[i]!=char2[i]) {distance++;}
		}
		return distance;
	}
	
	private static int distance(String seq1, String seq2) {
		char[] char1=seq1.toCharArray();
		char[] char2=seq2.toCharArray();
		
		int distance=0;
		for(int i=0; i<char1.length; i++) {
			if(char1[i]!=char2[i]) {distance++;}
		}
		return distance;
	}

	private List<String> getNames(List<Sequence> list) {
		List<String> rtrn=new ArrayList<String>();
		
		for(Sequence seq: list) {rtrn.add(seq.getName());}
		
		
		return rtrn;
	}
	
	private static Map<String, Integer> getDistance(Map<String, Integer> list1, int distance) {
		//make kmer sets
		
		
		//match within sets
		
		
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		int counter=0;
		for(String seq: list1.keySet()) {
			Collection<String> close= getSequenceWithinDistance(seq, list1.keySet(), distance);
			
			int sum=0;
			for(String c: close) {
				sum+=list1.get(c);
			}
			
			if(sum>1) {System.out.println(seq+"\t"+sum+"\t"+close.size());}
			
			rtrn.put(seq,  sum);
			counter++;
			if(counter%1000 ==0) {System.err.println(counter+" "+list1.size());}
		}
		return rtrn;
	}
	
	
	private static Collection<Cluster> getClose(Collection<String> list1, int distance) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		int counter=0;
		for(String seq: list1) {
			Collection<String> close= getSequenceWithinDistance(seq, list1, distance);
			if(close.size()>1) {
				Cluster c=new Cluster(close);
				rtrn.add(c);
			}
			counter++;
			if(counter%1000 ==0) {System.err.println(counter+" "+list1.size());}
		}
		return rtrn;
	}
	
	
	
	
	
	private static Collection<String> getSequenceWithinDistance(String seq, Collection<String> list1, int distance) {
		Collection<String> rtrn=new TreeSet<String>();
		for(String seq2: list1) {
			if(isClose(seq, seq2, distance)) {rtrn.add(seq2);}
		}
		return rtrn;
	}

	private static boolean isClose(String seq, String seq2, int maxDistance) {
		if(seq.length()!=seq2.length()) {return false;}
		
		int distance=distance(seq, seq2);
		//System.err.println(seq+" "+seq2+" "+distance);
		
		
		//compare in windows
		/*for(int i=0; i<(seq.length()-distance); i+=distance) {
			String sub1=seq.substring(i, i+distance);
			String sub2=seq2.substring(i, i+distance);
			if(distance(sub1, sub2)>=distance) {return false;}
			//if(!sub1.equals(sub2)) {return false;}
		}*/
		
		return distance<=maxDistance;
	}

	
	
	

	
	private static int distanceToClosest(String seq, Map<String, Integer> list1) {
		int minDistance=15;
		for(String seq2: list1.keySet()) {
			if(seq.length()==seq2.length()) {
				int dist=distance(seq, seq2);
				minDistance=Math.min(dist, minDistance);
			}
		}
		return minDistance;
	}

	private static void write(Map<String, Collection<String>> related, Map<String, Integer> list1) {
		for(String seq: related.keySet()) {
			int sum=sum(related.get(seq), list1);
			if(sum>1) {
				System.out.println(seq+"\t"+sum+"\t"+related.get(seq).size());
			}
		}
		
	}
	
	private static int sum(Cluster cluster, Map<String, Integer> list1) {
		int sum=0;
		for(String c: cluster.proteins) {
			sum+=list1.get(c);
		}
		return sum;
	}

	private static int sum(Collection<String> close, Map<String, Integer> list1) {
		int sum=0;
		for(String c: close) {
			sum+=list1.get(c);
		}
		return sum;
	}

	


	
	private Collection<Cluster> collapse(Map<String, Collection<String>> relatedSequences) {
		//go through each sequence and get all related as well as relateds
		
		
		Map<String, Cluster> updatedClusters=new TreeMap<String, Cluster>();
		for(String seq: relatedSequences.keySet()) {
			Cluster c=new Cluster(relatedSequences.get(seq));
			updatedClusters.put(seq, c);
		}
		
		
		Collection<Cluster> clusters=new TreeSet<Cluster>();
		
		for(String seq: relatedSequences.keySet()) {
			Collection<String> list=relatedSequences.get(seq);
			Cluster c=updatedClusters.get(seq);
			for(String related: list) {
				Cluster c1=updatedClusters.get(related);
				Cluster merged=merge(c, c1);
				updatedClusters.put(seq, merged);
			}
		}
		
		for(String seq: updatedClusters.keySet()) {
			Cluster c=updatedClusters.get(seq);
			clusters.add(c);
		}
		
		System.err.println(relatedSequences.size()+" "+clusters.size());
		return clusters;
	}

	
	private Cluster merge(Cluster c, Cluster c1) {
		Collection<String> list=new TreeSet<String>();
		list.addAll(c.proteins);
		list.addAll(c1.proteins);
		Cluster newCluster=new Cluster(list);
		return newCluster;
	}

	private static Map<String, Integer> parse(String fileName) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String line;
		int counter=0;
		while ((line = reader.readLine()) != null) {
			String seq=line.split("\t")[0];
			if(!seq.contains("*")) {
				int count=0; 
				if(rtrn.containsKey(seq)) {
					count=rtrn.get(seq);
				}
				count++;
				rtrn.put(seq, count);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		reader.close();
		
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			System.err.println("V3");
			File file=new File(args[0]);
			String save=args[1];
			ClusterProteinsBySequence.kmerSize=Integer.parseInt(args[2]);
			new ClusterProteinsBySequence(file, save);
		}
		else {System.err.println(usage);}
		
	}

	static String usage=" args[0]=files \n args[1]=save \n args[2]=kmer size (to search) ";
	
	//KFFEFWVDGTRCRARRGAMVKT

	
	
}

package guttmanlab.core.nanobody;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.sequence.Sequence;

public class SequenceDistance {

	/*public SequenceDistance(String input, int minDistance, String save) throws IOException {
		Map<String, Integer> list1=parse(input);
		Map<String, Collection<String>> relatedSequences=getSequencesWithinDistance(list1, minDistance);
		Collection<Cluster> clusters=collapse(relatedSequences, list1);
		write(save, clusters);
	}*/
	
	private void write(String save, Collection<Cluster> clusters) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		//TODO sort and name by rank
		int counter=0;
		for(Cluster c: clusters) {
			String name="cluster"+counter;
			c.setName(name);
			writer.write(c.toString()+"\n");
			counter++;
		}
		
		writer.close();
	}

	private Map<String, Collection<String>> getSequencesWithinDistance(Map<String, Integer> list1, int minDistance) throws IOException {
		//for each sequence get all within distance
		Map<String, Collection<String>> kmers=enumerateKmers(list1.keySet(), 9);
		
		Map<String, Collection<String>> related=new TreeMap<String, Collection<String>>();
		
		for(String k: kmers.keySet()) {
			Collection<String> list=kmers.get(k);
			Map<String, Collection<String>> map=getClose(list, minDistance);
			update(map, related);
		}
		
		return related;
		
	}

	private static Map<String, Collection<String>> enumerateKmers(Collection<String> list1, int k) {
		Map<String, Collection<String>> kmersToFull=new TreeMap<String, Collection<String>>();
		
		int counter=0;
		for(String seq: list1) {
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
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+list1.size());}
		}
		
		return kmersToFull;
	}

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
	
	
	private static Map<String, Collection<String>> getClose(Collection<String> list1, int distance) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		int counter=0;
		for(String seq: list1) {
			Collection<String> close= getSequenceWithinDistance(seq, list1, distance);
			rtrn.put(seq,  close);
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
		
		return distance<maxDistance;
	}

	public static Map<String, Collection<String>> getSequencesWithinDistance(String[] args) throws IOException {
		Map<String, Integer> list1=parse(args[0]);
		//for each sequence get all within distance
		//Map<String, Integer> counts=getDistance(list1, 3);
		Map<String, Collection<String>> kmers=enumerateKmers(list1.keySet(), 10);
		
		
		Map<String, Collection<String>> related=new TreeMap<String, Collection<String>>();
		
		for(String k: kmers.keySet()) {
			if(kmers.get(k).size()>1) {
				//System.err.println(k+" "+kmers.get(k).size());
				Collection<String> list=kmers.get(k);
				Map<String, Collection<String>> map=getClose(list, Integer.parseInt(args[1]));
				update(map, related);
			}
		}
		
		//write(related, list1);
		
		return related;
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

	private static void update(Map<String, Collection<String>> map, Map<String, Collection<String>> related) {
		for(String seq: map.keySet()) {
			Collection<String> list=new TreeSet<String>();
			if(related.containsKey(seq)) {
				list=related.get(seq);
			}
			list.addAll(map.get(seq));
			related.put(seq, list);
		}
	}


	
	/*private Collection<Cluster> collapse(Map<String, Collection<String>> relatedSequences, Map<String, Integer> list1) {
		//go through each sequence and get all related as well as relateds
		
		
		Map<String, Cluster> updatedClusters=new TreeMap<String, Cluster>();
		for(String seq: relatedSequences.keySet()) {
			Cluster c=new Cluster(seq, relatedSequences.get(seq));
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
			int sum=sum(c, list1);
			c.setScore(sum);
			clusters.add(c);
		}
		
		System.err.println(relatedSequences.size()+" "+clusters.size());
		return clusters;
	}*/

	
	/*private Cluster merge(Cluster c, Cluster c1) {
		Collection<String> list=new TreeSet<String>();
		list.addAll(c.proteins);
		list.addAll(c1.proteins);
		Cluster newCluster=new Cluster(c.getName()+"_"+c1.getName(), list);
		return newCluster;
	}*/

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
	
	
	/*public static void main(String[] args) throws IOException {
		if(args.length>2) {
			new SequenceDistance(args[0], Integer.parseInt(args[1]), args[2]);
		}
		else {System.err.println(usage);}
		
	}*/

	static String usage=" args[0]=list1 \n args[1]=max distance \n args[2]=save";

	
	
}

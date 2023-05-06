package guttmanlab.core.nanobody;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;

public class CompareSets {

	private static int distanceToClosest(String seq, Map<String, Integer> list1) {
		int minDistance=seq.length();
		String closest=seq;
		
		for(String seq2: list1.keySet()) {
			if(seq.length()==seq2.length()) {
				int dist=distance(seq, seq2);
				if(dist<minDistance) {
					minDistance=dist;
					closest=seq2;
				}
			}
		}
		return minDistance;
	}
	
	
	private static Map<Cluster, Integer> distanceToClosest(Collection<Cluster> clusters, Map<String, Integer> list1) {
		Map<Cluster, Integer> rtrn=new TreeMap<Cluster, Integer>();
		int counter=0;
		for(Cluster cluster: clusters) {
			int minDistance=100;
			for(String protein: cluster.proteins) {
				int distance=distanceToClosest(protein, list1);
				minDistance=Math.min(minDistance, distance);
			}
			cluster.setScore(minDistance);
			rtrn.put(cluster, minDistance);
			counter++;
			if(counter%100==0) {System.err.println(counter+" "+clusters.size());}
		}
		
		
		return rtrn;
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
	
	
	private static Map<String, Integer> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(String line: lines) {
			String seq=line.split("\t")[0];
			//String sub=seq.substring(seq.length()-10, seq.length());
			//System.err.println(sub+" "+sub.length());
			
			if(!seq.contains("*")) {
			
				int count=0; 
				if(rtrn.containsKey(seq)) {
					count=rtrn.get(seq);
				}
				count++;
				rtrn.put(seq, count);
			}
		}
		
		return rtrn;
	}
	
	private static void write(String string, Map<Cluster, Integer> distance) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Cluster c: distance.keySet()) {
			writer.write(c.toString()+"\n");
		}
		
		writer.close();
	}
	
	/*private static Collection<Cluster> parseClusters(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			int count=Integer.parseInt(tokens[2]);
			if(count>=2) {
				Collection<String> proteins=new TreeSet<String>();
				for(int i=3; i<tokens.length; i++) {
					proteins.add(tokens[i]);
				}
				Cluster newCluster=new Cluster(tokens[0], proteins);
				rtrn.add(newCluster);
			}
		}
		
		return rtrn;
	}*/
	
	
	/*public static void main(String[] args) throws IOException {
		Map<String, Integer> list1=parse(args[0]);
		Collection<Cluster> clusters=parseClusters(args[1]);
		Map<Cluster, Integer> distance=distanceToClosest(clusters, list1);
		
		write(args[2], distance);
	}*/


	


	
	
}

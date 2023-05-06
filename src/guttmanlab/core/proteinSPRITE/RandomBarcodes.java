package guttmanlab.core.proteinSPRITE;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Kmer;

public class RandomBarcodes {

	
	static int numPossible=24;
	static int numRounds=5;
	
	private static void enumeratePossibleBarcodes(int totalNumberOfClusters, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		Map<String, Integer> roundCount=new TreeMap<String, Integer>();
		Map<Kmer, Integer> uniqueKmers=new TreeMap<Kmer, Integer>();
		for(int i=0; i<totalNumberOfClusters; i++) {
			List<String> regions=random();
			add(regions, roundCount);
			Kmer k=makeKmer(regions);
			add(uniqueKmers, k);
			writer.write(k.toString(".")+"\n");
			double ratio=(double)uniqueKmers.size()/((double)i+1);
			int diff=(i+1)-uniqueKmers.size();
			if(i%10000==0) {System.err.println((i+1) +" "+uniqueKmers.size()+" "+ratio+" "+diff);}
		}
		writer.close();
		
		double ratio=(double)uniqueKmers.size()/(double)totalNumberOfClusters;
		int diff=totalNumberOfClusters-uniqueKmers.size();
		System.err.println(uniqueKmers.size()+" "+totalNumberOfClusters+" "+ratio+" "+diff);
		
		print(uniqueKmers);
		
		MatrixWithHeaders distanceMatrix=computeDistanceMatrix(uniqueKmers);
		distanceMatrix=filter(distanceMatrix);
		distanceMatrix.write(save+".distance");
		
	}
	
	private static MatrixWithHeaders filter(MatrixWithHeaders distanceMatrix) {
		Collection<String> rows=new TreeSet<String>();
		
		for(String r1: distanceMatrix.getRowNames()) {
			for(String r2: distanceMatrix.getColumnNames()) {
				double distance=distanceMatrix.get(r1, r2);
				if(!r1.equals(r2)) {
					if(distance>0) {
						rows.add(r1);
						rows.add(r2);
					}
				}
			}
		}
		
		System.err.println(distanceMatrix.getRowNames().size()+" "+rows.size());
		distanceMatrix=distanceMatrix.submatrixByRowNames(rows);
		distanceMatrix=distanceMatrix.submatrixByColumnNames(rows);
		return distanceMatrix;
	}

	private static MatrixWithHeaders computeDistanceMatrix(Map<Kmer, Integer> uniqueKmers) {
		List<String> list=toList(uniqueKmers);
		MatrixWithHeaders rtrn=new MatrixWithHeaders(list, list);
		
		int counter=0;
		for(Kmer k1: uniqueKmers.keySet()) {
			for(Kmer k2: uniqueKmers.keySet()) {
				int distance=distance(k1, k2);
				if(distance==0 || distance==1) {
					rtrn.set(k1.toString("."), k2.toString("."), k1.getSize()-distance);
				}
				
				if(distance==1) {System.out.println(k1.toString(".")+" "+k2.toString(".")+" "+distance);}
				counter++;
				if(counter%1000000==0) {System.err.println(counter);}
			}
			
			
		}
		
		
		return rtrn;
	}

	private static int distance(Kmer k1, Kmer k2) {
		return k1.distance(k2);
	}

	private static List<String> toList(Map<Kmer, Integer> uniqueKmers) {
		List<String> rtrn=new ArrayList<String>();
		
		for(Kmer k: uniqueKmers.keySet()) {
			rtrn.add(k.toString("."));
		}
		
		return rtrn;
	}

	private static void add(Map<Kmer, Integer> uniqueKmers, Kmer k) {
		int count=0;
		if(uniqueKmers.containsKey(k)) {
			count=uniqueKmers.get(k);
		}
		count++;
		uniqueKmers.put(k, count);
	}

	private static void print(Map<Kmer, Integer> roundCount) {
		for(Kmer r: roundCount.keySet()) {
			System.out.println(r.toString(";")+"\t"+roundCount.get(r));
		}
		
	}

	private static void add(List<String> regions, Map<String, Integer> roundCount) {
		for(String region: regions) {
			int count=0; 
			if(roundCount.containsKey(region)) {
				count=roundCount.get(region);
			}
			count++;
			roundCount.put(region, count);
		}
	}

	private static Kmer makeKmer(List<String> regions) {
		Kmer rtrn=new Kmer();
		
		for(String r: regions) {
			rtrn.addRegion(r);
		}
		
		return rtrn;
	}

	private static List<String> random() {
		List<String> rtrn=new ArrayList<String>();
		for(int i=0; i<numRounds; i++) {
			String round="ROUND"+i;
			int num=(int)(Math.random()*numPossible);
			String val=round+"_"+num;
			rtrn.add(val);
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		enumeratePossibleBarcodes(Integer.parseInt(args[0]), args[1]);
	}
	
}

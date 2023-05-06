package guttmanlab.core.proteinSPRITE;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.rnasprite.Kmer;

public class BarcodeDistances {

	int numPossible=24;
	int numRounds=6;
	int numPerm=100;
	
	public BarcodeDistances(int totalNumberOfClusters, String save) throws IOException {
		String[] barcodeFiles=makeSave(numRounds, save+".temp");
		int numberOfUniqueBarcodes=enumeratePossibleBarcodes(totalNumberOfClusters, barcodeFiles);
		
		double ratio=(double)numberOfUniqueBarcodes/(double)totalNumberOfClusters;
		System.err.println(totalNumberOfClusters+" "+numberOfUniqueBarcodes+" "+ratio);
		
		String finalFile=save+".final";
		int counts=collisionsK1(barcodeFiles, finalFile, numberOfUniqueBarcodes);
		
		ratio=(double)counts/(double)totalNumberOfClusters;
		
		System.err.println(counts+" "+ratio);
	}
	
	
	
	public BarcodeDistances(File input, String save) throws IOException {
		String[] barcodeFiles=makeSave(numRounds, save+".temp");
		
		int numberOfUniqueBarcodes=parseBarcodes(input, barcodeFiles);
		
		String finalFile=save+".final";
		int counts=collisionsK1(barcodeFiles, finalFile, numberOfUniqueBarcodes);
		
		double ratio=(double)counts/(double)numberOfUniqueBarcodes;
		
		System.err.println(counts+" "+ratio);
	}
	
	

	private int parseBarcodes(File input, String[] barcodeFiles) throws IOException {
		FileWriter[] writers=makeWriters(barcodeFiles);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		
		int count=0;
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			List<String> regions=parse(nextLine, "\\.");
			Kmer k=makeKmer(regions);
			for(int pos=0; pos<regions.size(); pos++) {
				Kmer sub=makeKmer(regions, regions.get(pos));
				writers[pos].write(sub.toString(";")+"\t"+k.toString()+"\n");
			}
			count++;
			if(count%100000==0) {System.err.println(count);}
		}
		reader.close();
		close(writers);
		return count;
	}



	private List<String> parse(String nextLine, String delim) {
		List<String> rtrn=new ArrayList<String>();
		String[] tokens=nextLine.split(delim);
		
		for(int i=0; i<tokens.length-1; i++) {rtrn.add(tokens[i]);}
		
		return rtrn;
	}



	private String[] makeSave(int numRounds2, String save) {
		String[] rtrn=new String[numRounds2];
		
		for(int i=0; i<numRounds2; i++) {
			rtrn[i]=save+"."+i;
		}
		
		return rtrn;
	}

	private int enumeratePossibleBarcodes(int totalNumberOfClusters, String[] barcodeFiles) throws IOException {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		FileWriter[] writers=makeWriters(barcodeFiles);
		for(int i=0; i<totalNumberOfClusters; i++) {
			List<String> regions=random();
			Kmer k=makeKmer(regions);
			for(int pos=0; pos<regions.size(); pos++) {
				Kmer sub=makeKmer(regions, regions.get(pos));
				writers[pos].write(sub.toString(";")+"\t"+k.toString()+"\n");
			}
			
			rtrn.add(k);
			if(i%100000==0) {System.err.println(i);}
		}
		close(writers);
		return rtrn.size();
	}

	private FileWriter[] makeWriters(String[] save) throws IOException {
		FileWriter[] rtrn=new FileWriter[save.length];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=new FileWriter(save[i]);
		}
		
		return rtrn;
	}

	private void close(FileWriter[] writers) throws IOException {
		for(int i=0; i<writers.length; i++) {writers[i].close();}
	}

	private Kmer makeKmer(List<String> regions) {
		Kmer rtrn=new Kmer();
		
		for(String r: regions) {
			rtrn.addRegion(r);
		}
		
		return rtrn;
	}
	
	private Kmer makeKmer(List<String> regions, String string) {
		Kmer rtrn=new Kmer();
		
		for(String r: regions) {
			if(!r.equals(string)) {rtrn.addRegion(r);}
		}
		
		return rtrn;
	}

	private List<String> random() {
		//Kmer k=new Kmer();
		List<String> rtrn=new ArrayList<String>();
		for(int i=0; i<numRounds; i++) {
			String round="ROUND"+i;
			long num=Math.round(Math.random()*numPossible);
			String val=round+"_"+num;
			rtrn.add(val);
		}
		return rtrn;
	}

	private int collisionsK1(String[] barcodeFiles, String save, int allBarcodes) throws IOException {
		Collection<Kmer> kmersToRemove=new TreeSet<Kmer>();
		
		for(int i=0; i<barcodeFiles.length; i++) {
			Map<Kmer, Collection<Kmer>> rtrn=new TreeMap<Kmer, Collection<Kmer>>();
			
			System.err.println(barcodeFiles[i]);
			Collection<Kmer> subset=getCounts(barcodeFiles[i]);
			
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFiles[i])));
		
			int count=0;
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				parse(nextLine, rtrn, subset);
				count++;
				if(count%100000==0) {System.err.println(count);}
			}
			kmersToRemove.addAll(filter(rtrn));
			reader.close();
		}
		
		
		for(Kmer k: kmersToRemove) {System.out.println(k.toString(";"));}
		
	
		//System.out.println(allBarcodes+" "+kmersToRemove.size());
		
		
		return kmersToRemove.size();
		
	}
	
	private Collection<Kmer> filter(Map<Kmer, Collection<Kmer>> temp) {
		//Map<Kmer, Collection<Kmer>> rtrn=new TreeMap<Kmer, Collection<Kmer>>();
		
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(Kmer k: temp.keySet()) {
			int size=temp.get(k).size();
			if(size>1) {rtrn.addAll(temp.get(k));}
		}
		
		return rtrn;
	}

	private Collection<Kmer> getCounts(String barcodeFile) throws IOException {
		Map<Kmer, Integer> counts=new TreeMap<Kmer, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		
		int counter=0;
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			parse(nextLine, counts);
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		reader.close();
		
		
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(Kmer k: counts.keySet()) {
			int count=counts.get(k);
			if(count>1) {rtrn.add(k);}
		}
		
		//System.err.println(rtrn.size());		
		return rtrn;
	}

	private Integer[] getReverse(Set<Integer> set) {
		Integer[] rtrn=new Integer[set.size()];
		
		int position=set.size()-1;
		for(Integer val: set) {
			rtrn[position]= val;
			position--;
		}
		
		return rtrn;
	}

	private int remove(Collection<Kmer> original, Collection<Kmer> kmersToRemove) {
		int count=0;
		
		for(Kmer k: original) {
			if(kmersToRemove.contains(k)) {
				kmersToRemove.remove(k);
				count++;
			}
		}
		
		return count;
	}

	private void add(int size, Kmer sub, Map<Integer, Collection<Kmer>> mergedSizes) {
		if(!mergedSizes.containsKey(size)) {mergedSizes.put(size, new TreeSet<Kmer>());}
		mergedSizes.get(size).add(sub);
	}

	private Collection<Kmer> hasSameSubs(Kmer kmer1, Map<Kmer, Collection<Kmer>> mergedKmers) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		Collection<Kmer> list1=mergedKmers.get(kmer1);
		for(Kmer kmer2: mergedKmers.keySet()) {
			if(!kmer1.equals(kmer2)) {
				Collection<Kmer> list2=mergedKmers.get(kmer2);
				if(equals(list1, list2)) {rtrn.add(kmer2);}
			}
		}
		return rtrn;
	}

	private boolean equals(Collection<Kmer> list1, Collection<Kmer> list2) {
		if(list1.size()!=list2.size()) {return false;}
		
		Iterator<Kmer> iter1=list1.iterator();
		Iterator<Kmer> iter2=list2.iterator();
		
		while(iter1.hasNext()) {
			Kmer k1=iter1.next();
			Kmer k2=iter2.next();
			if(!k1.equals(k2)) {return false;}
		}
		return true;
	}

	private void parse(String nextLine, Map<Kmer, Collection<Kmer>> rtrn, Collection<Kmer> subset) {
		String[] tokens=nextLine.split("\t");
		
		Kmer sub=parse(tokens[0]);
		
		if(subset.contains(sub)) {
		
			Kmer full=parse(tokens[1]);
			
			if(!rtrn.containsKey(sub)) {rtrn.put(sub, new TreeSet<Kmer>());}
			rtrn.get(sub).add(full);
			
			//if(rtrn.get(sub).size()>1) {listWithMultiple.add(sub);}
		}
	}
	
	private void parse(String nextLine, Map<Kmer, Integer> rtrn) {
		String[] tokens=nextLine.split("\t");
		
		Kmer sub=parse(tokens[0]);
		
		int count=0;
		
		if(rtrn.containsKey(sub)) {count=rtrn.get(sub);}
		count++;
		rtrn.put(sub, count);
	}
	
	/*private void parse(String nextLine, Map<Kmer, Collection<Kmer>> rtrn) {
		String[] tokens=nextLine.split("\t");
		
		Kmer k=parse(tokens[0]);
		Collection<Kmer> list=new ArrayList<Kmer>();
		for(int i=1; i<tokens.length; i++) {
			list.add(parse(tokens[i]));
		}
		
		for(Kmer sub: list) {
			if(!rtrn.containsKey(sub)) {rtrn.put(sub, new TreeSet<Kmer>());}
			rtrn.get(sub).add(k);
		}
		
	}*/

	private Kmer parse(String string) {
		String[] barcodes=string.split(";");
		Kmer rtrn=new Kmer();
		rtrn.addRegions(barcodes);
		return rtrn;
	}

	public static void main(String[] args) throws NumberFormatException, IOException {
		new BarcodeDistances(new File(args[0]), args[1]);
	}
	
}



package guttmanlab.core.nanobody;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.math.Statistics;

public class Merge {

	
	private static void write(String save, Map<String, Double>[] kmerCounts, File[] files) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<String> allKmers=new TreeSet<String>();
		for(int i=0; i<kmerCounts.length; i++) {allKmers.addAll(kmerCounts[i].keySet());}
		
		writer.write("kmer");
		for(int i=0; i<files.length; i++) {
			writer.write("\t"+files[i].getName());
		}
		writer.write("\n");
		
		for(String kmer: allKmers) {
			double[] vals=getVals(kmer, kmerCounts);
			writer.write(kmer);
			for(int i=0; i<vals.length; i++) {writer.write("\t"+vals[i]);}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	private static double[] getVals(String kmer, Map<String, Double>[] kmerCounts) {
		double[] rtrn=new double[kmerCounts.length];
		
		for(int i=0; i<kmerCounts.length; i++) {
			rtrn[i]=get(kmerCounts[i], kmer);
		}
		return rtrn;
	}
	
	private static double get(Map<String, Double> map, String kmer) {
		if(map.containsKey(kmer)) {return map.get(kmer);}
		return 0;
	}
	
	
	
	private static Map<String, Double>[] parse(File[] files, int minThreshold) throws NumberFormatException, IOException {
		Map<String, Double>[] rtrn=new Map[files.length];
		double[] total=new double[files.length];
		
		Collection<String> kmers=new TreeSet<String>();
		for(int i=0; i<files.length; i++) {
			total[i]=getKmers(files[i], minThreshold, kmers); //TODO get total counts and normalize the values
		}
		
		for(int i=0; i<files.length; i++) {rtrn[i]=parse(files[i], kmers, total[i]);}
		
		return rtrn;
	}
	
	
	private static double getKmers(File file, int minThreshold, Collection<String> kmers) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		double sum=0;
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String kmer=tokens[0];
			Integer val=Integer.parseInt(tokens[1]);
			sum+=val;
			if(val>=minThreshold) {
				rtrn.add(kmer);
			}
			if(counter%1000000==0) {System.err.println(counter);}
			counter++;
		}
		reader.close();
		kmers.addAll(rtrn);
		
		return sum;
	}

	private static Map<String, Double> parse(File file, Collection<String> allKmers, double total) throws NumberFormatException, IOException {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String kmer=tokens[0];
			Integer val=Integer.parseInt(tokens[1]);
			double norm=((double)val/total)*Math.pow(10, 7);
			if(allKmers.contains(kmer)) {
				rtrn.put(kmer, norm);
			}
			if(counter%1000000==0) {System.err.println(counter);}
			counter++;
		}
		reader.close();
		return rtrn;
	}
	
	private static Map<String, Integer> parse(File file) throws NumberFormatException, IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String kmer=tokens[0];
			Integer val=Integer.parseInt(tokens[1]);
			rtrn.put(kmer, val);
			if(counter%1000000==0) {System.err.println(counter);}
			counter++;
		}
		reader.close();
		return rtrn;
	}
	
	private static void filter(String file, int max, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			if(counter==0) {writer.write(nextLine+"\n");}
			else {
				String[] tokens=nextLine.split("\t");
				double[] vals=parse(tokens);
				if(Statistics.max(vals)>max) {writer.write(nextLine+"\n");}
			}
			if(counter%1000000==0) {System.err.println(counter);}
			counter++;
		}
		reader.close();
		writer.close();
	}

	private static double[] parse(String[] tokens) {
		double[] rtrn=new double[tokens.length-1];
		
		for(int i=1; i<tokens.length; i++) {
			rtrn[i-1]=Double.parseDouble(tokens[i]);
		}
		
		return rtrn;
	}
	
	private static double sum(Map<String, Integer> map) {
		double sum=0;
		
		for(String k: map.keySet()) {
			int val=map.get(k);
			sum+=val;
		}
		return sum;
	}


	public static void main(String[] args) throws NumberFormatException, IOException {
		if(args.length>2) {
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			int threshold=Integer.parseInt(args[2]);
			Map<String, Double>[] maps=parse(files, threshold);
			
			/*for(int i=0; i<maps.length; i++) {
				System.out.println(files[i].getName()+"\t"+sum(maps[i]));
			}*/
			
			write(save, maps, files);
		}
		else {System.err.println(usage);}
		
		
		//filter(args[0], Integer.parseInt(args[1]), args[2]);
		
	}
	
	static String usage=" args[0]=files \n args[1]=save \n args[2]=min threshold";

	
	

	
	
}

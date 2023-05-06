package guttmanlab.core.proteinSPRITE;

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
import java.util.TreeMap;
import java.util.TreeSet;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class MaxPercentByCluster {
	
	int binResolution=1000;

	public MaxPercentByCluster(File clusters, String save) throws IOException, InterruptedException{
		
		//Map<String, Double> percent=new TreeMap<String, Double>();
		Map<Double, Integer> percentCount=new TreeMap<Double, Integer>();
		
		int empty=0;
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(clusters)));
		int counter=0;
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			Map<String, Integer> proteins=new TreeMap<String, Integer>();
			String[] tokens=nextLine.split("\t");
			String barcode=tokens[0];
			for(int i=1; i<tokens.length; i++){
				if(tokens[i].startsWith("BEAD")){
					add(proteins, tokens[i].split(":")[0]);
				}
			}
			
			if(!proteins.isEmpty()){
				double maxPercent=maxPercent(proteins);
				add(percentCount, maxPercent);
				//percent.put(barcode, maxPercent);
			}
			else {empty++;}
			
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		reader.close();
		
		write(percentCount, save);
		System.out.println(empty+"\t"+counter);
	}
	
	private void add(Map<Double, Integer> percentCount, double maxPercent) {
		int count=0;
		if(percentCount.containsKey(maxPercent)){count=percentCount.get(maxPercent);}
		count++;
		percentCount.put(maxPercent, count);
		
	}

	private void write(Map<Double, Integer> percentCount, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Double percent: percentCount.keySet()){
			writer.write(percent+"\t"+percentCount.get(percent)+"\n");
		}
		
		writer.close();
	}

	private void binAndWrite(Map<String, Double> percent, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		double increment=0.01;
		for(double i=0; i<1; i+=increment){
			int count=count(percent, i, i+increment);
			double fraction=(double)count/(double)percent.size();
			writer.write(i+"\t"+count+"\t"+fraction+"\n");
		}
		writer.close();
	}

	private int count(Map<String, Double> percent, double i, double d) {
		int count=0;
		
		for(String b: percent.keySet()){
			double val=percent.get(b);
			if(val>=i && val<=d){count++;}
		}
		
		return count;
	}

	private void add(Map<String, Integer> proteins, String string) {
		int count=0;
		if(proteins.containsKey(string)){
			count=proteins.get(string);
		}
		count++;
		proteins.put(string, count);
		
	}

	private double maxPercent(Map<String, Integer> proteins) {
		int total=0;
		
		for(String p: proteins.keySet()){total+=proteins.get(p);}
		
		double max=0.0;
		for(String p: proteins.keySet()){
			int val=proteins.get(p);
			double percent=(double)val/(double)total;
			max=Math.max(max, percent);
		}
		
		return max;
	}

	private void add(Map<String, Integer> countsWithDNA, Collection<String> proteins) {
		 for(String protein: proteins){
			 int count=0;
			 if(countsWithDNA.containsKey(protein)){count=countsWithDNA.get(protein);}
			 count++;
			 countsWithDNA.put(protein, count);
		 }
	}

	private void write(Map<String, Integer> countsWithDNA, Map<String, Integer> countsWithoutDNA, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		Collection<String> proteins=new TreeSet<String>();
		proteins.addAll(countsWithoutDNA.keySet());
		proteins.addAll(countsWithDNA.keySet());
		
		
		for(String protein: proteins){
			int with=get(countsWithDNA,protein);
			int wo=get(countsWithoutDNA,protein);
			writer.write(protein+"\t"+with+"\t"+wo+"\n");
		}
		
		writer.close();
	}

	

	
	
	private int get(Map<String, Integer> countsWithDNA, String protein) {
		if(countsWithDNA.containsKey(protein)){
			return countsWithDNA.get(protein);
		}
		return 0;
	}

	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>1){
		File clusters=new File(args[0]);
		String save=args[1];
		
		new MaxPercentByCluster(clusters, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=save";
}

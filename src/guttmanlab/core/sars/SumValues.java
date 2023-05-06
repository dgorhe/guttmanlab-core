package guttmanlab.core.sars;

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

import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;

public class SumValues {

	
	private static void write(String save, Map<String, Integer> merged, long sumTotal, int binSize, long totalSize) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String key: merged.keySet()){
			writer.write(key+"\t"+merged.get(key)+"\t"+sumTotal+"\t"+binSize+"\t"+totalSize+"\n");
		}
		
		writer.close();
	}
	
	private static Map<String, Integer> merge(Map<String, Integer> vals, Map<String, Integer> merged) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(String key: vals.keySet()){
			int score1=get(merged, key);
			int score2=get(vals, key);
			int mergedScore=score1+score2;
			rtrn.put(key, mergedScore);
		}
		
		return rtrn;
	}
	
	private static int get(Map<String, Integer> merged, String key) {
		if(merged.containsKey(key)){return merged.get(key);}
		return 0;
	}

	private static Map<String, Integer> parse(File file) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String key=tokens[0];
			int val=new Integer(tokens[1]);
			rtrn.put(key, val);
		}
		reader.close();
		return rtrn;
	}
	
	private static long getTotal(File file) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine=reader.readLine();
		reader.close();
		String[] tokens=nextLine.split("\t");
		long val=new Long(tokens[2]);
		return val;
	}
	
	private static int getBinSize(File file) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine=reader.readLine();
		reader.close();
		String[] tokens=nextLine.split("\t");
		int val=new Integer(tokens[3]);
		return val;
	}
	
	private static long getTotalSize(File file) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine=reader.readLine();
		reader.close();
		String[] tokens=nextLine.split("\t");
		int val=new Integer(tokens[4]);
		return val;
	}
	
	private static void writeQuantile(String save, Map<String, Integer>[] counts, long[] totalCounts, int binSize, long totalSize, double quantile) throws IOException {
		FileWriter writer=new FileWriter(save);
			
		Collection<String> allKeys=getKeys(counts);
		
		
		//double[] totalVals=getCounts(totalCounts);
		for(String key: allKeys){
			int[] countVals=getCounts(counts, key);
			int index=getQuantile(countVals, totalCounts, quantile);
			writer.write(key+"\t"+countVals[index]+"\t"+totalCounts[index]+"\t"+binSize+"\t"+totalSize+"\n");
		}
			
		writer.close();
		
		
	}
	
	
	private static void writeSum(String save, Map<String, Integer>[] counts, long[] totalCounts, int binSize, long totalSize) throws IOException {
		FileWriter writer=new FileWriter(save);
			
		Collection<String> allKeys=getKeys(counts);
		
		for(String key: allKeys){
			int[] countVals=getCounts(counts, key);
			//double[] totalVals=getCounts(totalCounts, key);
			
			
			writer.write(key+"\t"+sum(countVals)+"\t"+sum(totalCounts)+"\t"+binSize+"\t"+totalSize+"\n");
		}
			
		writer.close();
		
		
	}
	
	
	private static void writeAverage(String save, Map<String, Integer>[] counts, long[] totalCounts, int binSize, long totalSize) throws IOException {
		FileWriter writer=new FileWriter(save);
			
		Collection<String> allKeys=getKeys(counts);
		
		for(String key: allKeys){
			int[] countVals=getCounts(counts, key);
			//double[] totalVals=getCounts(totalCounts, key);
			
			
			writer.write(key+"\t"+average(countVals)+"\t"+average(totalCounts)+"\t"+binSize+"\t"+totalSize+"\n");
		}
			
		writer.close();
		
		
	}
	
	private static double average(double[] countVals) {
		return sum(countVals)/(double)countVals.length;
	}
	
	private static int average(int[] countVals) {
		return new Double(sum(countVals)/(double)countVals.length).intValue();
	}
	
	private static long average(long[] countVals) {
		return new Double(sum(countVals)/(double)countVals.length).longValue();
	}

	
	private static double sum(double[] vals){
		double sum=0;
		for(int i=0; i<vals.length; i++){sum+=vals[i];}
		return sum;
	}
	
	private static int sum(int[] vals){
		int sum=0;
		for(int i=0; i<vals.length; i++){sum+=vals[i];}
		return sum;
	}
	
	private static long sum(long[] vals){
		long sum=0;
		for(int i=0; i<vals.length; i++){sum+=vals[i];}
		return sum;
	}
	
	private static Collection<String> getKeys(Map<String, Integer>[] counts) {
		Collection<String> rtrn=new TreeSet<String>();
		for(int i=0; i<counts.length; i++){rtrn.addAll(counts[i].keySet());}
		return rtrn;
	}

	private static int getQuantile(int[] countVals, long[] totalVals, double quantile) {
		double[] ratio=new double[countVals.length];
		for(int i=0; i< ratio.length; i++){
			ratio[i]=getRatio(countVals[i], totalVals[i]);
		}
		int index=getMedianPosition(ratio, quantile);
		//System.err.println(index+" "+countVals[index]+" "+totalVals[index]);
		return index;
	}

	//TODO non-zero??
	private static int getMedianPosition(double[] ratio, double quantile) {
		double val=Statistics.quantile(ratio.clone(), quantile);
		//System.err.println(val);
		//System.err.println(val+" "+Statistics.min(ratio)+" "+Statistics.max(ratio)+" "+Statistics.quantile(ratio, 0.5)+" "+Statistics.quantile(ratio, 0.9));
		for(int i=0; i<ratio.length; i++){
			if(ratio[i]==val){return i;}
		}
		return 0;
	}

	private static double getRatio(int d, long e) {
		return (double)d/(double)e;
	}

	private static int[] getCounts(Map<String, Integer>[] counts, String key) {
		int[] vals=new int[counts.length];
		
		for(int i=0; i<counts.length; i++){
			int score=0;
			if(counts[i].containsKey(key)){score=counts[i].get(key);}
			vals[i]=score;
		}
		
		return vals;
	}
	
	private static double[] getCounts(int [] counts) {
		double[] vals=new double[counts.length];
		
		for(int i=0; i<counts.length; i++){
			vals[i]=counts[i];
		}
		
		return vals;
	}

	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		
		long sumTotal=0;
		int binSize=0;
		long totalSize=0;
		Map<String, Integer> merged=new TreeMap<String, Integer>();
		
		Map<String, Integer>[] counts=new Map[files.length];
		long[] totalCounts=new long[files.length];
		
		for(int i=0; i<files.length; i++){
			counts[i]=parse(files[i]);
			totalCounts[i]=getTotal(files[i]);
			System.err.println(files[i].getName()+" "+totalCounts[i]);
			
			binSize=Math.max(binSize, getBinSize(files[i]));
			totalSize=Math.max(totalSize, getTotalSize(files[i]));
		}
		
		writeQuantile(save+".median", counts, totalCounts, binSize, totalSize, 0.5);
		writeQuantile(save+".90Percent", counts, totalCounts, binSize, totalSize, 0.9);
		writeSum(save+".sum", counts, totalCounts, binSize, totalSize);
		writeAverage(save+".average", counts, totalCounts, binSize, totalSize);
		
		//write(save, merged, sumTotal, binSize, totalSize);
		
	}

	

	

	

	

	

	
	
}

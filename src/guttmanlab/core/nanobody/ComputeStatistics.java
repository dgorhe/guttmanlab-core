package guttmanlab.core.nanobody;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class ComputeStatistics {
	
	List<String> names;
	Map<String, List<Integer>> map;
	List<Integer> totals;
	List<Integer> uniqueSequenceCounts;
	
	public ComputeStatistics(File input, String saveDir, List<String> names) throws IOException {
		this.names=names;
		parse(input);
		
		for(int index=0; index<names.size()-1; index++) {
			System.err.println("running index "+index);
			String save=saveDir+"/"+names.get(index)+"."+names.get(index+1)+".scores";
			significanceByIndex(index, save);
		}
		
	}
	
	private void significanceByIndex(int index, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int newTotal=sumIfIndex1(map, index);
		
		for(String cluster: map.keySet()) {
			List<Integer> vals=map.get(cluster);
			double freq2=computeFreq(vals.get(index), newTotal);
			int n=totals.get(index+1);
			int k=vals.get(index+1);
			
			BinomialDistribution d=new BinomialDistribution(n, freq2);
			double p2=1-d.cumulativeProbability(k);
			
			double expected=(double)n*freq2;
			double enrichment=(double)k/expected;
			
			if(p2<0.05) {
				writer.write(cluster+"\t"+vals.get(index)+"\t"+vals.get(index+1)+"\t"+enrichment+"\t"+p2+"\n");
			}
			
		}
		
		writer.close();
		System.err.println(totals.get(index)+"\t"+newTotal);
		
	}

	private int sumIfIndex1(Map<String, List<Integer>> map2, int index) {
		int rtrn=0;
		
		for(String cluster: map2.keySet()) {
			List<Integer> vals=map2.get(cluster);
			if(vals.get(index+1)>0) {rtrn+=vals.get(index);}
		}
		
		return rtrn;
	}

	private double computeFreq(int val, int total) {
		double newVal=Math.max(1, val);
		return newVal/(double)total;
	}

	private void parse(File input) throws NumberFormatException, IOException {
		this.map=new TreeMap<String, List<Integer>>();
		
		List<Integer> sum=initializeI(names.size(), 0);
		List<Integer> counts=initializeI(names.size(), 0);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens=nextLine.split("\t");
			String cluster=tokens[0];
			List<Integer> vals=new ArrayList<Integer>();
			for(int i=1; i<tokens.length; i++) {
				int val=Integer.parseInt(tokens[i]);
				vals.add(val);
				int newSum=sum.get(i-1)+val;
				sum.set(i-1, newSum);
				int count=0;
				if(val>0) {count=1;}
				int newCount=counts.get(i-1)+count;
				counts.set(i-1, newCount);
			}
			map.put(cluster, vals);
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		reader.close();
		
		
		for(int i=0; i<names.size(); i++) {
			System.err.println(names.get(i)+"\t"+sum.get(i)+"\t"+counts.get(i));
		}
		
		this.totals=sum;
		this.uniqueSequenceCounts=counts;
		
	}

	
	
	private List<Integer> initializeI(int size, int initialVal) {
		List<Integer> rtrn=new ArrayList<Integer>();
		
		for(int i=0; i<size; i++) {
			rtrn.add(initialVal);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		File input=new File(args[0]);
		String saveDir=args[1];
		
		List<String> names=new ArrayList<String>();
		names.add(args[2]);
		names.add(args[3]);
		
		new ComputeStatistics(input, saveDir, names);
	}
	
}

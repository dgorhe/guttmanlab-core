package guttmanlab.core.spidr;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

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
import java.util.TreeSet;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class NormalizeCountMatrix {
	
	static int numPerm=100;
	static double minFold=2;

	public static MatrixWithHeaders normalize(MatrixWithHeaders matrix) {
		MatrixWithHeaders norm=new MatrixWithHeaders(matrix.getRowNames(), matrix.getColumnNames());
		double total=getTotal(matrix);
		Map<String, Double> totals=getTotals(matrix);
		
		//Collection<String> set=new TreeSet<String>();
		
		//go column by column
		int count=1;
		for(String experiment: matrix.getColumnNames()) {
			double experimentTotal=matrix.get("total", experiment);
			double fraction=experimentTotal/(total-experimentTotal);
			System.err.println(experiment+" "+fraction+" "+count+" "+matrix.getColumnNames().size());
			for(String region: matrix.getRowNames()) {
				double observed=matrix.get(region, experiment);
				double sum=totals.get(region)-observed;
				double expected=percentile(sum, fraction, 0.95);
				double enrich=observed/expected;
				//if(enrich>minFold) {set.add(region);}
				norm.set(region, experiment, enrich);
			}
			count++;
		}
		
		//System.err.println(norm.getRowNames().size()+" "+set.size());
		
		//norm=norm.submatrixByRowNames(set);
		
		//System.err.println("filtered "+norm.getRowNames().size()+" "+set.size());
		
		return norm;
	}
	
	
	
	
	public static void normalize(File input, String save, double minCount) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		
		double total=0; // all counts for all samples
		Map<String, Double> totalPerSample=null; //total for each sample
		String[] samples=null;
		
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			if(counter==0){
				samples=tokens;
				writer.write("sample");
				for(int i=1; i<tokens.length; i++) {writer.write("\t"+tokens[i]);}
				writer.write("\n");
			}
			else if(counter==1) {
				//TOTAL
				totalPerSample=getTotals(tokens, samples);
				total=Statistics.sum(totalPerSample.values());
			}
			else{
				String region=tokens[0];
				double sum=sum(tokens);
				if(sum>minCount) {
					writer.write(region);
					for(int i=1; i<tokens.length; i++) {
						double experimentTotal=totalPerSample.get(samples[i]);
						double fraction=experimentTotal/(total-experimentTotal);
						double observed=Double.parseDouble(tokens[i]);
						double expected=percentile((sum-observed), fraction, 0.95); //TODO bootstrap samples
						double enrich=observed/expected;
						writer.write("\t"+enrich);
					}
					writer.write("\n");
				}
			}
			
			
			
			
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		reader.close();
		writer.close();
	}
	

	private static double sum(String[] tokens) {
		double rtrn=0;
		for(int i=1; i<tokens.length; i++) {
			rtrn+=Double.parseDouble(tokens[i]);
		}
		return rtrn;
	}




	private static Map<String, Double> getTotals(String[] tokens, String[] samples) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(int i=1; i<tokens.length; i++) {
			rtrn.put(samples[i], Double.parseDouble(tokens[i]));
		}
		
		return rtrn;
	}




	private static Map<String, Double> getTotals(MatrixWithHeaders matrix) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(String region: matrix.getRowNames()) {
			double[] vals=matrix.getRow(region);
			double sum=Statistics.sum(vals);
			rtrn.put(region, sum);
		}
		
		return rtrn;
	}

	private static double getTotal(MatrixWithHeaders matrix) {
		double[] vals=matrix.getRow("total");
		return Statistics.sum(vals);
	}

	private static double[] getPerm(double sum, double fraction) {
		double[] rtrn=new double[numPerm];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=perm(sum, fraction);
		}
		
		return rtrn;
	}

	private static double percentile(double sum, double fraction, double percentile) {
		double mean=sum*fraction;
		double rtrn=1;
		if(mean>0) {
			PoissonDistribution dist=new PoissonDistribution(mean);
			rtrn=Math.max(rtrn,dist.inverseCumulativeProbability(percentile));
		}
		return rtrn;
	}
	
	private static double perm(double sum, double fraction) {
		double mean=sum*fraction;
		PoissonDistribution dist=new PoissonDistribution(mean);
		return dist.sample();
	}

	private static double getSum(MatrixWithHeaders matrix, String region, String experiment) {
		double[] vals=matrix.getRow(region);
		double val=matrix.get(region, experiment);
		return Statistics.sum(vals)-val;
	}
	
	private static MatrixWithHeaders filterMatrix(File input, int min) throws IOException {
		List<String> columns=new ArrayList<String>();
		List<String> rows=new ArrayList<String>();
		
		List<String> lines=new ArrayList<String>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			if(counter>0) {
				double sum=0;
				for(int i=1; i<tokens.length; i++) {
					sum+=Double.parseDouble(tokens[i]);
				}
				if(sum>min) {
					rows.add(tokens[0]);
					lines.add(nextLine);
					
					//System.out.println(nextLine);
				}
			}
			else {
				for(int i=1; i<tokens.length; i++) {
					columns.add(tokens[i]);
				}
				
				//System.out.println(nextLine);
			}
			counter++;
			if(counter%10000==0) {System.err.println(counter);}
		}
		reader.close();
		
		MatrixWithHeaders rtrn=makeMatrix(rows, columns, lines);
		System.err.println("done making matrix");
		return rtrn;
		
	}
	
	private static MatrixWithHeaders makeMatrix(List<String> rows, List<String> columns, List<String> lines) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			String row=tokens[0];
			for(int i=1; i<tokens.length; i++) {
				String column=columns.get(i-1);
				double val=Double.parseDouble(tokens[i]);
				rtrn.set(row, column, val);
			}
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			File file=new File(args[0]);
			double minCount=Double.parseDouble(args[1]);
			normalize(file, args[2], minCount);
		}
		else {System.err.println(usage);}
		
	}
	
	static String usage=" args[0]=count matrix \n args[1]=min count \n args[2]=save";
	
}

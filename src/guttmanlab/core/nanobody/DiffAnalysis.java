package guttmanlab.core.nanobody;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;

import guttmanlab.core.math.Statistics;

public class DiffAnalysis {
	
	double alpha=0.0001;
	
	public DiffAnalysis(File file, List<Double> sums) throws NumberFormatException, IOException {
		//Map<String, List<Double>> map=new TreeMap<String, List<Double>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		
		double p=getP(sums);
		System.err.println(p);
		
		List<Double> sum=initialize(4,0);
		
		while ((nextLine = reader.readLine()) != null) {
			if(counter!=0) {
				String[] tokens=nextLine.split(",");
				String sequence=tokens[1];
				List<Double> vals=new ArrayList<Double>();
				for(int i=2; i<tokens.length; i++) {
					double val=Double.parseDouble(tokens[i]);
					vals.add(val);
					add(sum, val, i-2);
				}
				//map.put(sequence, vals);	
				double sum1=vals.get(0)+vals.get(1);
				double sum2=vals.get(2)+vals.get(3);
				
				int n=new Double(Statistics.sum(vals)).intValue();
				int k=new Double(sum1).intValue();
				
				
				BinomialDistribution d=new BinomialDistribution(n, p);
				double p2=1-d.cumulativeProbability(k);
				double pval=Math.min(p2, 1-p2);
				
				double expected=(double)n*p;
				double enrichment=sum1/expected;
				
				double fold=(sum1+1)/(sum2+1);
				
				if(pval>0.25)  {
					System.out.println(sequence+"\t"+sum1+"\t"+sum2+"\t"+enrichment+"\t"+fold+"\t"+pval);
				}
			}
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		reader.close();
		
		
		/*for(int i=0; i<sum.size(); i++) {
			System.out.println(i+"\t"+sum.get(i));
		}*/
		
	}
	
	

	private double getP(List<Double> sums) {
		double sum=Statistics.sum(sums);
		double num=sums.get(0)+sums.get(1);
		return num/sum;
	}



	private List<Double> initialize(int size, double initialVal) {
		List<Double> rtrn=new ArrayList<Double>();
		
		for(int i=0; i<size; i++) {
			rtrn.add(initialVal);
		}
		
		return rtrn;
	}
	
	private void add(List<Double> sum, double val, int i) {
		double newSum=sum.get(i)+val;
		sum.set(i, newSum);
	}


	public static void main(String[] args) throws NumberFormatException, IOException {
		File file=new File(args[0]);
		List<Double> sums=new ArrayList<Double>();
		sums.add(5711061.0);
		sums.add(5672723.0);
		sums.add(5705348.0);
		sums.add(5702849.0);
		
		new DiffAnalysis(file, sums);
	}
	
}

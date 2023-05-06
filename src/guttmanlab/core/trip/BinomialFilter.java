package guttmanlab.core.trip;

import java.io.IOException;
import java.util.List;

import org.apache.commons.math3.distribution.BinomialDistribution;

import guttmanlab.core.annotation.io.BEDFileIO;

public class BinomialFilter {

	public BinomialFilter(String input) throws IOException {
		List<String> lines=BEDFileIO.loadLines(input);
		
		for(String line: lines) {
			String[] tokens=line.split(" ");
			String barcode=tokens[0];
			
			int totalCount1=Integer.parseInt(tokens[1])+Integer.parseInt(tokens[2]);
			int totalCount2=Integer.parseInt(tokens[3])+Integer.parseInt(tokens[4]);
			int spliced1=Integer.parseInt(tokens[2]);
			int spliced2=Integer.parseInt(tokens[4]);
			double p=(double)(spliced1+spliced2)/(double)(totalCount1+totalCount2);
			
			double p1=binomialProb(totalCount1, p, spliced1);
			double p2=binomialProb(totalCount2, p, spliced2);
			
			double min=Math.min(p1, p2);
			double max=Math.max(p1, p2);
			min=Math.min(min, 1-max);
			
			
			double r1=(double)spliced1/(double)totalCount1;
			double r2=(double)spliced2/(double)totalCount2;
			
			
			
			
			
			double pComp=binomialProb(totalCount2, r1, spliced2);
			double pComp2=binomialProb(totalCount1, r2, spliced1);
			
			//System.out.println(barcode+"\t"+totalCount1+"\t"+totalCount2+"\t"+r1+"\t"+r2+"\t"+min);
			
			System.out.println(barcode+"\t"+min+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+tokens[4]);
		}
		
	}
	
	
	private static double binomialProb(int totalCount, double p, int observedCount) {
		BinomialDistribution dist=new BinomialDistribution(totalCount, p);
		//return dist.probability(observedCount);
		return dist.cumulativeProbability(observedCount);
	}
	
	public static void main(String[] args) throws IOException {
		new BinomialFilter(args[0]);
	}
	
}

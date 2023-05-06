package guttmanlab.core.util;

import jsc.distributions.Binomial;

public class SplitAmplifyBinomial {

	public static void main(String[] args) {
		double startingP=Double.parseDouble(args[0]);
		int numberOfCells=Integer.parseInt(args[1]);
		int numberOfWells=Integer.parseInt(args[2]);
		int numberOfRounds=10;
		
		/*
		//int startingK=new Double(numberOfCells*startingP).intValue();
		//System.err.println(startingK);
		
		
		double p=startingP;
		System.out.println((0)+"\t"+p);
		//int k=startingK;
		for(int i=0; i<numberOfRounds; i++) {
			if(p<1) {
				Binomial b=new Binomial(numberOfCells, p);
				//int k=startingK*2;
				
				
				int k=getK(numberOfCells, p, numberOfWells, b);
				double newP=((double)k/(double)numberOfCells);
				double enrich=newP/p;
				System.out.println((i+1)+"\t"+newP);
				p=newP;
			}
			else {System.out.println((i+1)+"\t"+p);}
			//startingK=k;
		}*/
		
		
		Binomial b=new Binomial(numberOfCells, startingP);
		int k=2;
		double newP=1-b.cdf(k);
		System.err.println(newP);
		
	}

	private static int getK(int numberOfCells, double p, int numberOfWells, Binomial b) {
		for(int i=0; i<numberOfCells; i++) {
			double newP=1-b.cdf(i);
			double expected=newP*numberOfWells;
			//System.err.println(i+" "+expected);
			if(expected<1) {return i;}
		}
		return numberOfCells;
	}
		
		/*Binomial b=new Binomial(16, 1.0/16.0);
		int k=2;
		double newP=1-b.cdf(k);
		System.err.println(newP);
		}*/
	
}

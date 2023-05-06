package guttmanlab.core.splicing.speckle.kineticmodel;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.math3.distribution.PoissonDistribution;


public class ProductionModel {
	static int maxTime=240;
	static int numIter=4;

	private static Map<Double, Double> fitProductionModel(int geneLength, double initiationRate, double elongationRate, double labelTime, double decayRate) {
			//TODO Compute time step from initiation rate and elongation rate
			double timeStep=1.0/initiationRate;
			
			//Initial conditions
			double numberOfPolymerasesAtInitial=initialNumbers(initiationRate, elongationRate, geneLength);
			List<Double> polymerasePositions=new ArrayList<Double>();
			for(int i=0; i<numberOfPolymerasesAtInitial; i++) {
				//randomly pick position on gene
				double genePosition=randomGenePosition(geneLength);
				polymerasePositions.add(genePosition);
				//System.err.println(genePosition);
			}
			
			Map<Double, Integer> completedAtTime=new TreeMap<Double, Integer>();
			
			//Step through labeling time interval
			for(double time=-labelTime; time<=maxTime; time+=timeStep) {
				double ir=0;
				if(time<=0) {ir=initiationRate;}
				List<Double> positions=iterate(polymerasePositions, ir, elongationRate, timeStep);
				int numberCompleted=getCompleted(positions, geneLength);
				completedAtTime.put(time, numberCompleted);
				polymerasePositions=getRemainingPolymerasePositions(positions, geneLength);
			}
			
		
			//Decay
			//TODO Add decay
			return cumulativeValues(completedAtTime, decayRate, timeStep);
	}
	

	private static Map<Double, Double> cumulativeValues(Map<Double, Integer> completedAtTime, double decayRate, double timeStep) {
		Map<Double, Double> rtrn=new TreeMap<Double, Double>();
		double sum=0;
		double updatedSum=0;
		double cumalitiveDecay=0;
		for(Double time: completedAtTime.keySet()) {
			sum+=completedAtTime.get(time);
			double mean=decayRate*timeStep*updatedSum;
			double remove=0;
			if(mean>0) {remove=sample(mean);}
			//if(mean>0) {remove=mean;}
			updatedSum-=remove;
			updatedSum+=completedAtTime.get(time);
			cumalitiveDecay+=remove;
			//System.out.println(time+"\t"+completedAtTime.get(time)+"\t"+sum+"\t"+updatedSum+"\t"+remove+"\t"+cumalitiveDecay); //TODO this should be variables
			rtrn.put(time, updatedSum);
		}
		return rtrn;
	}

	private static double randomGenePosition(int length) {
		double rand=Math.random();
		return rand*length;
	}

	private static List<Double> iterate(List<Double> polymerasePositions, double initiationRate, double elongationRate, double timeStep) {
		List<Double> rtrn=new ArrayList<Double>();
		double distanceToMove=elongationRate*timeStep;
		//System.err.println(distanceToMove);
		//For each existing polymerase, elongate at elongation rate
		for(Double pos: polymerasePositions) {
			double newPos=pos+distanceToMove;
			rtrn.add(newPos);
		}	
		
		if(initiationRate>0) {
			//Initiate from start at initiation rate
			double numToInitiate=initiationRate*timeStep;
			int sample=sample(numToInitiate);
			//int sample=1;
			
			for(int i=0; i<sample; i++) {
				double newPos=0;
				rtrn.add(newPos);
			}
		}
		
		return rtrn;
	}

	private static int sample(double numToInitiate) {
		PoissonDistribution dist=new PoissonDistribution(numToInitiate);
		return dist.sample();
	}

	private static List<Double> getRemainingPolymerasePositions(List<Double> positions, int genomicLength) {
		List<Double> rtrn=new ArrayList<Double>();
		for(double position: positions) {
			if(position<=genomicLength) {
				rtrn.add(position);
			}
		}
		
		return rtrn;
	}

	private static int getCompleted(List<Double> positions, int genomicLength) {
		int counter=0;
		for(double position: positions) {
			if(position>genomicLength) {counter++;}
		}
		return counter;
	}

	private static double initialNumbers(double initiationRate, double elongationRate, int genomicLength) {
		double timeToComplete=genomicLength/elongationRate;
		double meanPolymerasesPerGene=initiationRate*timeToComplete;
		
		//System.err.println(timeToComplete+" "+meanPolymerasesPerGene);
		return meanPolymerasesPerGene;
	}
	
	private static double compare(double[] x, double[] y, Map<Double, Double> fit) {
		double ss=0;
		for(int i=0; i<x.length; i++) {
			double time=x[i];
			double observed=y[i];
			double model=getClosest(fit, time);
			ss+=(Math.abs(model-observed));
		}
		return ss;
	}
	
	private static double getClosest(Map<Double, Double> fit, double time) {
		double timeAtMin=Double.MAX_VALUE;
		double min=Double.MAX_VALUE;
		
		for(Double t: fit.keySet()) {
			double diff=Math.abs(t-time);
			if(diff<min) {
				min=diff;
				timeAtMin=t;
			}
		}
		
		return fit.get(timeAtMin);
	}

	
	private static double[] optimizeFit(double[] initiationRates, double[] decayRates, int geneLength, double elongationRate, double labelTime, double[] x, double[] y) {
		double[] initialConditions= {initiationRates[0], decayRates[0]};
		for(int i=0; i<numIter; i++) {
			double initiationRate=minimize(initiationRates, initialConditions[1], geneLength, elongationRate, labelTime, x, y);
			initialConditions[0]=initiationRate;
			double decayRate=minimize(initialConditions[0], decayRates, geneLength, elongationRate, labelTime, x, y);
			initialConditions[1]=decayRate;
			//initiationRate=minimize(initiationRates, decayRate, geneLength, elongationRate, labelTime, x, y);
			
			Map<Double, Double> fit=fitProductionModel(geneLength, initiationRate, elongationRate, labelTime, decayRate);
			double ss=compare(x,y,fit);
			System.err.println(i+" "+initiationRate+" "+decayRate+" "+ss);
		}
		return initialConditions;
	}


	private static double minimize(double initiationRate, double[] decayRates, int geneLength, double elongationRate,double labelTime, double[] x, double[] y) {
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		
		for(int i=0; i<decayRates.length; i++) {
			double decayRate=decayRates[i];
			Map<Double, Double> fit=fitProductionModel(geneLength, initiationRate, elongationRate, labelTime, decayRate);
			double ss=compare(x,y,fit);
			if(ss<minSS) {
				minSS=ss;
				rtrn=decayRate;
			}
		}
		return rtrn;
	}


	private static double minimize(double[] initiationRates, double decayRate, int geneLength, double elongationRate, double labelTime, double[] x, double[] y) {
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		
		for(int i=0; i<initiationRates.length; i++) {
			double initiationRate=initiationRates[i];
			Map<Double, Double> fit=fitProductionModel(geneLength, initiationRate, elongationRate, labelTime, decayRate);
			double ss=compare(x,y,fit);
			System.err.println(initiationRate+" "+decayRate+" "+ss);
			if(ss<minSS) {
				minSS=ss;
				rtrn=initiationRate;
			}
		}
		return rtrn;
	}


	public static void main(String[] args) {
		int geneLength=3971;
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		double[] y= {511,687,850,803,632,501,213,131,82,63,31,56};
		double elongationRate=3000;
		double labelTime=10;
		
		double[] initiationRates= {75,80,85,90,95,100,105,110,115,120};
		double[] decayRates= {0.02,0.025, 0.03, 0.035, 0.04, 0.05};
		//double initiationRate=300;
		//double decayRate=0.1;
		
		double[] initialConditions=optimizeFit(initiationRates, decayRates, geneLength, elongationRate, labelTime, x, y);
		
		Map<Double, Double> fit=fitProductionModel(geneLength, initialConditions[0], elongationRate, labelTime, initialConditions[1]);
		for(Double t: fit.keySet()) {
			System.out.println(t+"\t"+fit.get(t));
		}
		
		
		/*for(int i=0; i<initiationRates.length; i++) {
			double initiationRate=initiationRates[i];
			for(int j=0; j<decayRates.length; j++) {
				double decayRate=decayRates[j];
				Map<Double, Integer> fit=fitProductionModel(geneLength, initiationRate, elongationRate, labelTime, decayRate);
				double ss=compare(x,y,fit);
				System.err.println(initiationRate+"\t"+decayRate+"\t"+ss);
			}
		}*/
		
		
	}


	
		
	



	
	
}

package guttmanlab.core.splicing.speckle.kineticmodel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.PoissonDistribution;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;


public class ExonModel {
	static int maxTime=240;
	static int numIter=4;

	private static Map<Double, List<Integer>> fitProductionModel(int geneLength, double initiationRate, double elongationRate, double labelTime, double decayRate) {
		Map<Double, List<Integer>> rtrn=new TreeMap<Double, List<Integer>>();
		
		double timeStep=1.0/initiationRate;
			
		//Initial conditions
		double numberOfPolymerasesAtInitial=initialNumbers(initiationRate, elongationRate, geneLength);
		List<Double> polymerasePositions=new ArrayList<Double>();
		for(int i=0; i<numberOfPolymerasesAtInitial; i++) {
			//randomly pick position on gene
			double genePosition=randomGenePosition(geneLength);
			polymerasePositions.add(genePosition);
		}
			
		//Step through labeling time interval
		int numberCompleted=0;
		for(double time=-labelTime; time<=maxTime; time+=timeStep) {
			/**
			 * Step 1: Move the elongating polymerase by elongation step
			 * Remove all that are now completed
			 * Keep those still elongating for next round
			 */
			double ir=0;
			if(time<=0) {ir=initiationRate;}
			List<Double> positions=iterate(polymerasePositions, ir, elongationRate, timeStep);
			int newlyCompleted=getCompleted(positions, geneLength);
			polymerasePositions=getRemainingPolymerasePositions(positions, geneLength);
			
			
			/**
			 * Step 2: Sample from the decay rate to decide how many completed transcripts to remove
			 */
			//how many to degrade
			double mean=decayRate*timeStep*numberCompleted;
			int remove=0;
			if(mean>0) {remove=sample(mean);}
			numberCompleted=numberCompleted-remove;
			numberCompleted+=newlyCompleted;
			
			List<Integer> list=makeList(polymerasePositions, numberCompleted, geneLength);
			//System.out.println(time+"\t"+numberCompleted);
			rtrn.put(time, list);
		}
		
		
		return rtrn;
	}
	
	
	
	private static void fitProductionModel(Gene transcript, double initiationRate, double elongationRate, double labelTime, double decayRate, double splicingRate) {
		Map<Double, List<Integer>> rtrn=new TreeMap<Double, List<Integer>>();
		
		int geneLength=transcript.getGenomicLength();
		double timeStep=1.0/initiationRate;
			
		//Initial conditions
		double numberOfPolymerasesAtInitial=initialNumbers(initiationRate, elongationRate, geneLength);
		List<Double> polymerasePositions=new ArrayList<Double>();
		List<Gene> stillTranscribing=new ArrayList<Gene>();
		for(int i=0; i<numberOfPolymerasesAtInitial; i++) {
			//randomly pick position on gene
			double genePosition=randomGenePosition(geneLength);
			if(genePosition<transcript.getGenomicLength()) {
				polymerasePositions.add(genePosition);
				
				//Test whether to splice the introns before the current Pol II position by sampling from the rate distribution
				//Use this to generate initial transcript set
				double absolutePosition=genePosition+transcript.getReferenceStartPosition();
				Gene generatedTranscript=getTranscript(absolutePosition, transcript, splicingRate, elongationRate);
				stillTranscribing.add(generatedTranscript);
				
				SingleInterval region=new SingleInterval(transcript.getReferenceName(), transcript.getReferenceStartPosition(), Double.valueOf(absolutePosition).intValue());
				System.out.println(region.toBED());
			}
		}
			
		
		
	}
	

	//TODO This does NOT work yet
	private static Gene getTranscript(double genePosition, Gene transcript, double splicingRate, double elongationRate) {
		//TODO Test whether to splice the introns before the current Pol II position by sampling from the rate distrinution
		//TODO Use this to generate initial transcript set
		
		Collection<Gene> junctions=transcript.getJunctions();
		Collection<SingleInterval> blocks=new TreeSet<SingleInterval>();
		
		//Iterate through all introns, 
		for(Gene junction: junctions) {
			Annotation intron=junction.getIntrons().iterator().next();
			//if intron.end < gene position
			if(intron.getReferenceEndPosition()< genePosition) {
				//test for splicing
				double distance=genePosition-intron.getReferenceEndPosition();
				double timeStep=distance/elongationRate;
				double spliceProb=splicingRate*timeStep;
				if(Math.random()<spliceProb) {
					//splice
					blocks.add(junction.getFirstBlock());
					
					SingleInterval second=new SingleInterval(junction.getReferenceName(), intron.getReferenceEndPosition(),  Double.valueOf(genePosition).intValue());
					
					blocks.add(second);
				}
				else {
					blocks.add(junction.getSingleInterval());
				}
				
			}
		}
		
		BlockedAnnotation annotation=new BlockedAnnotation(blocks);
		Gene rtrn=new Gene(annotation);
		return rtrn;
	}



	private static List<Integer> makeList(List<Double> polymerasePositions, int numberCompleted, int geneLength) {
		List<Integer> rtrn=new ArrayList<Integer>();
		
		for(Double position: polymerasePositions) {
			int pos= Double.valueOf(position).intValue();
			rtrn.add(pos);
		}
		for(int i=0; i< numberCompleted; i++) {rtrn.add(geneLength);}
		
		return rtrn;
	}


	private static void write(double time, List<Double> polymerasePositions, int numberCompleted) {
		for(Double pos: polymerasePositions) {
			System.out.println("18S\t0\t"+ Double.valueOf(pos).intValue());
		}
		
		for(int i=0; i<numberCompleted; i++) {
			System.out.println("18S\t0\t1870");
		}
		
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

	


	private static MatrixWithHeaders makeCountMatrix(Map<Double, List<Integer>> vals, Gene gene) {
		Map<Integer, String> positions=new TreeMap<Integer, String>();
		List<String> columns=new ArrayList<String>();
		Iterator<SingleInterval> blocks=gene.getBlocks();
		int i=1;
		while(blocks.hasNext()) {
			SingleInterval exon=blocks.next();
			String name="exon_"+i;
			columns.add(name);
			int distance=gene.getReferenceEndPosition()-exon.getReferenceEndPosition();
			positions.put(distance, name);
			i++;
		}
		
		Iterator<Annotation> introns=gene.getIntrons().iterator();
		i=1;
		while(introns.hasNext()) {
			SingleInterval intron=introns.next().getSingleInterval();
			String name="intron_"+i;
			columns.add(name);
			int distance=gene.getReferenceEndPosition()-intron.getReferenceEndPosition();
			positions.put(distance, name);
			i++;
		}
		
		List<String> rows=new ArrayList<String>();
		for(Double time: vals.keySet()) {
			rows.add( Double.valueOf(time).toString());
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(Double time: vals.keySet()) {
			String row= Double.valueOf(time).toString();
			List<Integer> list=vals.get(time);
			for(int position: list) {
				for(int position2: positions.keySet()) {
					if(position>=position2) {
						String name=positions.get(position2);
						rtrn.incrementCount(row, name);
					}
				}
			}
		}
		
		
		return rtrn;
		
	}
	

	public static void main(String[] args) throws IOException {
		int geneLength=7069;
		double elongationRate=3000;
		double labelTime=20;
		
		double initiationRate=75;
		double decayRate=0.02;
		double splicingRate=0.2;
		
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		Gene transcript=genes.iterator().next();
		
		Map<Double, List<Integer>> vals=fitProductionModel(geneLength, initiationRate, elongationRate, labelTime, decayRate);
		
		MatrixWithHeaders matrix=makeCountMatrix(vals, transcript);
		
		matrix.write(args[1]);
		
		
		//fitProductionModel(transcript, initiationRate, elongationRate, labelTime, decayRate, splicingRate);
	}



	


	
		
	



	
	
}

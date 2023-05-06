package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.predicate.StrandFilter;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import jsc.distributions.Poisson;
import net.sf.samtools.util.CloseableIterator;

public class SimulateEnrichment {

	private double numberOfMolecules=6600000000.0;
	private double percentSignal=.5;
	private double enrichment=10.0;
	private double yield=0.0001;
	private int windowSize=100;
	private double inputYield=0.01; //TODO Replace
	
	private Random random;
	
	public SimulateEnrichment(File bamFile, String save, double yield, double percentSignal) throws IOException{
		this.random=new Random(12345678);
		this.yield=yield;
		this.percentSignal=percentSignal;
		
		//Compute abundance
		Map<Annotation, Double> abundance=computeEnrichment(bamFile);
		
		//Randomly generate signal positions
		Collection<Annotation> trueSignal=generateSignal(abundance);
		
		//Deplete total number of molecules by yield (poisson sample from dist set to yield)
		Map<Annotation, Double> inputToIP=depleteStartingMaterial(abundance, this.yield);
		
		//Enrich signal
		Map<Annotation, Double> enrichmentSample=enrich(inputToIP, trueSignal);
		
		//Input is sampled at 1/100 total
		Map<Annotation, Double> inputSample=depleteStartingMaterial(abundance, inputYield);
		
		//write enrichment sample and input
		write(save, enrichmentSample, inputSample, trueSignal);
	}
	
	
	
	
	private Map<Annotation, Double> enrich(Map<Annotation, Double> abundance, Collection<Annotation> trueSignal) {
		Map<Annotation, Double> rtrn=new TreeMap<Annotation, Double>();
		Poisson enriched=new Poisson(this.enrichment); //TODO Replace with Normal
		Poisson background=new Poisson(1.0);
		
		for(Annotation window: abundance.keySet()){
			double factor=1.0;
			Poisson dist=background;
			if(trueSignal.contains(window)){
				dist=enriched;
				factor=10.0;
			}
			//double factor=dist.random();
			double newVal=abundance.get(window)*factor;
			System.out.println(dist.mean()+" "+factor+" "+abundance.get(window)+" "+newVal);
			rtrn.put(window, newVal);
		}
		return rtrn;
	}




	private void write(String save, Map<Annotation, Double> enrichmentSample, Map<Annotation, Double> inputSample, Collection<Annotation> trueSignal) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Annotation window: enrichmentSample.keySet()){
			double enrichment=enrichmentSample.get(window);
			double input=inputSample.get(window);
			double ratio=enrichment/input;
			boolean enriched=false;
			if(trueSignal.contains(window)){enriched=true;}
			writer.write(window.toUCSC(window.getOrientation())+"\t"+enrichment+"\t"+input+"\t"+ratio+"\t"+enriched+"\n");
		}
		
		writer.close();
	}




	private Map<Annotation, Double> depleteStartingMaterial(Map<Annotation, Double> abundance, double yield) {
		//Deplete total number of molecules by yield (poisson sample from dist set to yield)
		Map<Annotation, Double> rtrn=new TreeMap<Annotation, Double>();
		
		
		for(Annotation window: abundance.keySet()){
			double numMolecules=(abundance.get(window)*this.numberOfMolecules*yield);
			Poisson p=new Poisson(numMolecules);
			double sampleYield=p.random();
			rtrn.put(window, sampleYield);
			
		}
		return rtrn;
	}




	private Map<Annotation, Double> computeEnrichment(File bamFile) {
		Map<Annotation, Double> negative=computeEnrichment(bamFile, Strand.NEGATIVE);
		Map<Annotation, Double> positive=computeEnrichment(bamFile, Strand.POSITIVE);
		positive.putAll(negative);
		
		Map<Annotation, Double> merged=convertToFractionTotal(positive);
		return merged;
	}


	private Map<Annotation, Double> convertToFractionTotal(Map<Annotation, Double> map) {
		Map<Annotation, Double> rtrn=new TreeMap<Annotation, Double>();
		
		double sum=0;
		for(Annotation w: map.keySet()){
			sum+=map.get(w);
		}
		
		for(Annotation w: map.keySet()){
			double val=map.get(w)/sum;
			rtrn.put(w, val);
		}
		
		return rtrn;
	}




	private Map<Annotation, Double> computeEnrichment(File bamFile1, Strand strand){
		BAMSingleReadCollection bam1=new BAMSingleReadCollection(bamFile1);
		bam1.addFilter(new StrandFilter(strand));
			
		Map<Annotation, Double> rtrn=new TreeMap<Annotation, Double>();
			
		//Iterate over windows
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows1=bam1.getPopulatedWindows(windowSize, windowSize);
			
		int counter=0;
		while(windows1.hasNext()){
			PopulatedWindow<SAMFragment> window=windows1.next();
			double count=window.getNumberOfAnnotationsInWindow(); //TODO Compute TPM
			double norm=(double)count/(double)window.size();
			Annotation newAnnotation=new SingleInterval(window.getReferenceName(), window.getReferenceStartPosition(), window.getReferenceEndPosition(), strand);
			rtrn.put(newAnnotation, norm);
			counter++;
			if(counter%10000 ==0){System.err.println(counter+" "+window.toUCSC(strand));}
		}
			
		windows1.close();
		return rtrn;
	}


	private void write(String save, Collection<Annotation> trueSignal) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Annotation signal: trueSignal){
			writer.write(signal.toBED()+"\n");
		}
		
		writer.close();
	}




	private Collection<Annotation> generateSignal(Map<Annotation, Double> abundance) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		for(Annotation window: abundance.keySet()){
			double random=this.random.nextDouble();
			if(random<percentSignal){
				rtrn.add(window);
			}
		}
		return rtrn;
	}




	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File bamFile=new File(args[0]);
			String save=args[1];
			double yield=new Double(args[2]);
			double percentSignal=new Double(args[3]);
			SimulateEnrichment e=new SimulateEnrichment(bamFile, save, yield, percentSignal);
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage="args[0]=bam file \n args[1]=save \n args[2]=yield \n args[3]=percent signal";
}

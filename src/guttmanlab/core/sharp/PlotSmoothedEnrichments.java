package guttmanlab.core.sharp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.barcodeidentification.PeakCalling;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PlotSmoothedEnrichments {
	
	double alpha=0.001;
	double minEnrichment=3.0;
	int minCount=20;
	Map<Annotation, Integer> maxInputScore;
	double percentileToUse=.9;
	
	public PlotSmoothedEnrichments(File sampleBAM, File inputBAM, int binSize, Map<String, IntervalTree<Gene>> genes, String save, double alpha) throws IOException {
		this.alpha=alpha;
		Map<SingleInterval, Integer> sampleCounts=getCounts(sampleBAM, binSize, genes);
		Map<SingleInterval, Integer> inputCounts=getCounts(inputBAM, binSize, genes);
		
		Map<Gene, Collection<SingleInterval>> windowsByGene=makeGeneWindows(genes, binSize);
		Map<Gene, Integer> inputGeneScores=inputGeneScores(windowsByGene, inputCounts);
		
		
		
		Pair<Integer> total1=getTotalCount(sampleBAM, binSize);
		Pair<Integer> total2=getTotalCount(inputBAM, binSize);
		
		
		Map<SingleInterval, Double> windowEnrichment=getWindowEnrichment(sampleCounts, inputCounts, inputGeneScores, windowsByGene, total1, total2);
		
		writeSmoothedBedgraph(windowEnrichment, save);
		
		System.err.println("done writing");
	}
	
	
	
	
	private Pair<Integer> getTotalCount(File bam, int binSize) {
		SAMFileReader inputReader= new SAMFileReader(bam);
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		Collection<SingleInterval> windows=new TreeSet<SingleInterval>();
		
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				totalCount++;
				windows.addAll(SAMFragment.allBins(read, binSize));
			}
			if(totalCount%1000000 ==0){System.err.println(totalCount+" "+read.getReferenceName());}
		}
					
		reads.close();
		inputReader.close();
		
		Pair<Integer> rtrn=new Pair<Integer>(totalCount, windows.size());
		
		return rtrn;
	}




	private Map<SingleInterval, Double> getWindowEnrichment(Map<SingleInterval, Integer> counts1,Map<SingleInterval, Integer> counts2, Map<Gene, Integer> inputGeneScores, Map<Gene, Collection<SingleInterval>> windowsByGene, Pair<Integer> total1, Pair<Integer> total2) {
		Map<SingleInterval, Double> differentialEnrich=new TreeMap<SingleInterval, Double>();
		
		//this.maxInputScore=getInputScores(inputPeaks, samplePeaks, geneScores);
		
		for(Gene gene: windowsByGene.keySet()) {
			int inputGeneScore=inputGeneScores.get(gene);
			for(SingleInterval window: windowsByGene.get(gene)) {
				int window1=get(counts1,window);
				int window2=get(counts2, window);
				window2=Math.max(1, window2);
				window2=Math.max(window2, inputGeneScore);
				
				double enrichment=getWindowNormEnrichment(window1, window2, total1.getValue1(), total2.getValue1(), total1.getValue2(), total2.getValue2());
				double windowP=getWindowNormPValue(window1, window2, total1.getValue1(), total2.getValue1(), total1.getValue2(), total2.getValue2());
				
				
				//double enrichment=getEnrichment(window1, window2, total1, total2);
				if(windowP<alpha) {enrichment=0;}
				if(differentialEnrich.containsKey(window)) {
					enrichment=Math.min(differentialEnrich.get(window), enrichment);
				}
				differentialEnrich.put(window, enrichment);
			}
		}
		

		return differentialEnrich;
	}


	private double getWindowNormPValue(int sampleReadCount, int inputReadCount, int sampleTotalCount, int inputTotalCount, int sampleWindows, int inputWindows) {
		int k=sampleReadCount;
		int n=sampleReadCount+inputReadCount;
		
		if(n>0){
			double elution=(double)sampleTotalCount/(double)sampleWindows;
			double input=(double) inputTotalCount/(double)inputWindows;
			double p=elution/(elution+input);
			Binomial b=new Binomial(n, p);
			return 1-b.cdf(k);
		}
		return 1.0;
	}



	private Map<Gene, Integer> inputGeneScores(Map<Gene, Collection<SingleInterval>> windowsByGene, Map<SingleInterval, Integer> inputCounts) {
		Map<Gene, Integer> rtrn=new TreeMap<Gene, Integer>();
		
		for(Gene gene: windowsByGene.keySet()) {
			Collection<SingleInterval> windows=windowsByGene.get(gene);
			List<Double> vals=getVals(inputCounts, windows);
			double percentile=Statistics.quantile(vals, percentileToUse);
			int intval=new Double(percentile).intValue()+1;
			rtrn.put(gene, intval);
		}
		
		return rtrn;
	}




	private List<Double> getVals(Map<SingleInterval, Integer> inputCounts, Collection<SingleInterval> windows) {
		List<Double> rtrn=new ArrayList<Double>();
		
		for(SingleInterval window: windows) {
			double val=get(inputCounts, window);
			rtrn.add(val);
		}
		
		return rtrn;
	}




	private int get(Map<SingleInterval, Integer> inputCounts, SingleInterval window) {
		if(inputCounts.containsKey(window)) {return inputCounts.get(window);}
		return 0;
	}




	private Map<SingleInterval, Integer> getCounts(File bam, int binSize, Map<String, IntervalTree<Gene>> genes) {
		Map<SingleInterval, Integer> positionCount=new TreeMap<SingleInterval, Integer>();
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				boolean overlapsGene=overlapsGenes(read, genes);
				if(overlapsGene) {
					Collection<SingleInterval> allBins=SAMFragment.getAllWindows(read, binSize);
					for(SingleInterval bin: allBins) {
						//System.err.println(bin);
						int score=0;
						if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
						score++;
						positionCount.put(bin, score);
					}
				}
			}
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount+" "+read.getReferenceName());}
		}
					
		//this.totalCount=totalCount;
		reads.close();
		inputReader.close();
		//this.numberOfWindows=numWindows(positionCount, genes.keySet());
			
		return positionCount;
	}
		
		
	
	private boolean overlapsGenes(SAMRecord read, Map<String, IntervalTree<Gene>> genes) {
		if(genes.containsKey(read.getReferenceName())) {
			IntervalTree<? extends Annotation> tree=genes.get(read.getReferenceName());
			Iterator<? extends Annotation> iter=tree.overlappingValueIterator(read.getAlignmentStart(), read.getAlignmentEnd());
			while(iter.hasNext()) {
				Annotation g=iter.next();
				boolean overlaps=overlapsExons(g, read);
				boolean overlapsStrand=overlapsStrand(g, read);
				if(overlaps && overlapsStrand) {return true;}
			}
		}
		return false;
	}
		
	private Collection<Annotation> getGenes(SAMRecord read, Map<String, IntervalTree<Annotation>> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		if(genes.containsKey(read.getReferenceName())) {
			IntervalTree<? extends Annotation> tree=genes.get(read.getReferenceName());
			Iterator<? extends Annotation> iter=tree.overlappingValueIterator(read.getAlignmentStart(), read.getAlignmentEnd());
			while(iter.hasNext()) {
				Annotation g=iter.next();
				boolean overlaps=overlapsExons(g, read);
				boolean overlapsStrand=overlapsStrand(g, read);
				if(overlaps && overlapsStrand) {rtrn.add(g);}
			}
		}
		return rtrn;
	}
	
	private boolean overlapsStrand(Annotation g, SAMRecord read) {
		Strand readStrand=getStrand(read);
		if(g.getOrientation().equals(readStrand)) {return true;}
		return false;
	}

	private static Strand getStrand(SAMRecord read) {
		Strand s=Strand.UNKNOWN;
		if(!read.getFirstOfPairFlag()) {
			if(read.getReadNegativeStrandFlag()) {s=Strand.NEGATIVE;}
			else{s=Strand.POSITIVE;}
		}
		else {
			if(read.getReadNegativeStrandFlag()) {s=Strand.POSITIVE;}
			else{s=Strand.NEGATIVE;}
		}
		return s;
	}


	private boolean overlapsExons(Annotation g, SAMRecord read) {
		SAMFragment f=new SAMFragment(read);
		//return f.overlaps(g);
		return g.fullyContained(f);
		
	}




	private void writeSmoothedBedgraph(Map<SingleInterval, Double> windowEnrichment, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval bin: windowEnrichment.keySet()) {
			double score=windowEnrichment.get(bin);
			SingleInterval mid=bin.getMidPoint();
			writer.write(mid.toBedgraph(score)+"\n");
		}
		
		writer.close();
	}




	
	
	
	

	private static Collection<Annotation> makeGeneWindows(Collection<Gene> genes, int binSize) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		for(Gene gene: genes) {
			Collection<SingleInterval> bins=gene.getAllWindows(binSize);
			rtrn.addAll(bins);
		}
		return rtrn;
	}

	
	private static Map<Gene, Collection<SingleInterval>> makeGeneWindows(Map<String, IntervalTree<Gene>> genes, int binSize) {
		Map<Gene, Collection<SingleInterval>> rtrn=new TreeMap<Gene, Collection<SingleInterval>>();
		for(String chr: genes.keySet()) {
			Iterator<Gene> iter=genes.get(chr).valueIterator();
			while(iter.hasNext()) {
				Gene gene=iter.next();
				Collection<SingleInterval> bins=gene.getAllWindows(binSize);
				for(SingleInterval bin: bins) {System.out.println(bin.toBED());}
				rtrn.put(gene, bins);
			}
		}
		return rtrn;
	}



	
	
	
	
	


	public double getWindowNormEnrichment(double sampleCount, double maxInputCount, int sampleTotalCount, int inputTotalCount, int sampleWindows, int inputWindows){
		double numerator=(double)sampleCount/(double)maxInputCount;
		double elution=(double)sampleTotalCount/(double)sampleWindows;
		double input=(double)inputTotalCount/(double)inputWindows;
		double denominator=elution/input;
		return numerator/denominator;
	}


	
	
	private static double getEnrichment(double window1, double window2, int total1, int total2) {
		double num=((double)window1)/(double)total1;
		double denom=((double)window2)/(double)total2;
		return num/denom;
	}
	
	private static double getPValue(int window_sampleCount, int window_inputCount, int total_sampleCount, int total_inputCount) {
		//Return binomial p
		int n=window_inputCount+window_sampleCount;
		if(n>0){
			double p=(double)total_sampleCount/(double)(total_inputCount+total_sampleCount);
			Binomial b=new Binomial(n, p);
			return 1-b.cdf(window_sampleCount);
		}
		return 1.0;
	}
	
	private static int readCount(File bam) {
		SAMFileReader inputReader= new SAMFileReader(bam);
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
		reads.close();
		inputReader.close();
		return totalCount;
	}


	public static void main(String[] args) throws IOException {
		if(args.length>5) {
			File sampleBam=new File(args[0]);
			File inputBam=new File(args[1]);
			int binSize=Integer.parseInt(args[2]);
			
			String save=args[4];
			
			Map<String, IntervalTree<Gene>> genes=BEDFileIO.loadTreePlusIntrons(args[3]);
			double alpha=Double.parseDouble(args[5]);
			//.loadTreeGenomeBin(args[3]);
			
			System.err.println("V6");
			new PlotSmoothedEnrichments(sampleBam, inputBam, binSize, genes, save, alpha);
		}
		else {System.err.println(usage);}
		
		
	}
	
	

	static String usage=" args[0]=sample bam \n args[1]=input bam \n args[2]=bin size \n args[3]=genes (bed) \n args[4]=save \n args[5]=pval cutoff";
}

package guttmanlab.core.spidr;

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
	
	public PlotSmoothedEnrichments(File sampleBAM, File[] controls, double p, int binSize, int numPerm, int smooth, String save, String chr) throws IOException {
		Map<SingleInterval, Double> actualScores=getCounts(sampleBAM, binSize, chr);
		
		Map<SingleInterval, double[]> randomScores=getRandomScores(controls, p, binSize, numPerm, actualScores, chr);
		
		
		writeActual(save, actualScores);
		
		
		/*Map<SingleInterval, Double> smoothed=getWindowEnrichment(actualScores, smooth);
		writeSmoothedBedgraph(smoothed, save+".smoothed.bedgraph");
		
		System.err.println("done writing");*/
	}
	
	
	private Map<SingleInterval, double[]> getRandomScores(File[] bamFiles, double p, int binSize, int numPerm, Map<SingleInterval, Double> actualScores, String chr) {
		Map<SingleInterval, double[]> scores=initialize(actualScores, numPerm);
		
		int counter=0;
		for(int i=0; i<bamFiles.length; i++) {
			SAMFileReader reader=new SAMFileReader(bamFiles[i]);
			SAMRecordIterator reads=reader.iterator();
			while(reads.hasNext()){
				
				SAMRecord read=reads.next();
				
				if(read.getReferenceName().startsWith(chr) && !read.getReadUnmappedFlag() && read.getMappingQuality()>1){
					boolean[] use=sample(p, numPerm);
					Collection<SingleInterval> allBins=SAMFragment.getAllWindows(read, binSize);
					for(SingleInterval bin: allBins) {
						if(scores.containsKey(bin.getMidPoint())) {
							double[] vals=scores.get(bin.getMidPoint());
							vals=increment(vals, use);
							scores.put(bin.getMidPoint(), vals);
						}
					}
				}
				counter++;
				if(counter%1000000 ==0){System.err.println(counter);}
			}
			
			reader.close();
			reads.close();
		}
		return scores;
	}
	
	private boolean[] sample(double p, int numPerm) {
		boolean[] rtrn=new boolean[numPerm];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=sample(p);
		}
		
		return rtrn;
	}
	
	private double[] increment(double[] vals, boolean[] use) {
		double[] rtrn=new double[vals.length];
		
		for(int i=0; i<vals.length; i++) {
			rtrn[i]=vals[i];
			if(use[i]) {rtrn[i]++;}
		}
		
		return rtrn;
	}
	
	private boolean sample(double p) {
		return Math.random()<p;
	}
	
	private Map<SingleInterval, double[]> initialize(Map<SingleInterval, Double> actualScores, int numPerm) {
		Map<SingleInterval, double[]> rtrn=new TreeMap<SingleInterval, double[]>();
		
		for(SingleInterval r: actualScores.keySet()) {
			double[] vals=new double[numPerm];
			rtrn.put(r, vals);
		}
		
		return rtrn;
	}
	
	private void writeActual(String save, Map<SingleInterval, Double> actualScores) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(SingleInterval r: actualScores.keySet()) {
			double score=actualScores.get(r);
			writer.write(r.toBedgraph(score)+"\n");
		}
		writer.close();
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




	/*private Map<SingleInterval, Double> getWindowEnrichment(Map<SingleInterval, Integer> counts1,Map<SingleInterval, Integer> counts2, Map<Gene, Integer> inputGeneScores, Map<Gene, Collection<SingleInterval>> windowsByGene, Pair<Integer> total1, Pair<Integer> total2) {
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
	}*/
	
	
	private Map<SingleInterval, Double> getWindowEnrichment(Map<SingleInterval, Double> counts1, int smoothWindow) {
		Map<SingleInterval, Double> differentialEnrich=new TreeMap<SingleInterval, Double>();
		
		//Map<String, IntervalTree<SingleInterval>> trees=makeTree(counts1);
		//for each window, get all other windows overlapping
		
		//TODO IF count matrix was by midpoint then we can just look for smooth window up and down
		for(SingleInterval window: counts1.keySet()) {
			//int mid=window.getMidPoint().getReferenceStartPosition();
			//SingleInterval newWindow=new SingleInterval(window.getReferenceName(), mid-smoothWindow, mid+smoothWindow);
			List<Double> vals=getWindows(window, counts1, smoothWindow);
			double val=Statistics.mean(vals);
			differentialEnrich.put(window, val);
		}
		

		return differentialEnrich;
	}



	private Map<String, IntervalTree<SingleInterval>> makeTree(Map<SingleInterval, Double> counts1) {
		return BEDFileIO.loadSingleIntervalTree(counts1.keySet());
	}


	
	private List<Double> getWindows(SingleInterval window, Map<SingleInterval, Double> counts1, int size) {
		List<Double> rtrn=new ArrayList<Double>();
		String chr=window.getReferenceName();
		int start=window.getReferenceStartPosition();
		
		for(int i=0; i<=(size/2); i++) {
			SingleInterval w1=new SingleInterval(chr, start+i, start+i+1);
			SingleInterval w2=new SingleInterval(chr, start-i, (start-i)+1);
			double score1=get(counts1, w1);
			double score2=get(counts1, w2);
			rtrn.add(score1);
			rtrn.add(score2);
		}
		
		return rtrn;
	}


	private List<Double> getWindows(SingleInterval window, Map<SingleInterval, Double> counts1, Map<String, IntervalTree<SingleInterval>> trees) {
		List<Double> rtrn=new ArrayList<Double>();
		if(trees.containsKey(window.getReferenceName())) {
			IntervalTree<SingleInterval> tree=trees.get(window.getReferenceName());
			Iterator<SingleInterval> windows=tree.overlappingValueIterator(window.getReferenceStartPosition(), window.getReferenceEndPosition());
			while(windows.hasNext()) {
				rtrn.add(counts1.get(windows.next()));
			}
		}
		return rtrn;
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



	

	private double get(Map<SingleInterval, Double> inputCounts, SingleInterval window) {
		if(inputCounts.containsKey(window)) {return inputCounts.get(window);}
		return 0;
	}




	private Map<SingleInterval, Double> getCounts(File bam, int binSize, String chr) {
		Map<SingleInterval, Double> positionCount=new TreeMap<SingleInterval, Double>();
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(read.getReferenceName().startsWith(chr) && !read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				Collection<SingleInterval> allBins=SAMFragment.getAllWindows(read, binSize);
				for(SingleInterval bin: allBins) {
					double score=0;
					if(positionCount.containsKey(bin.getMidPoint())){score=positionCount.get(bin.getMidPoint());}
					score++;
					positionCount.put(bin.getMidPoint(), score);
				}
			}
			totalCount++;
			if(totalCount%100000 ==0){System.err.println(totalCount+" "+read.getReferenceName());}
		}
					
		reads.close();
		inputReader.close();

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
		if(args.length>4) {
			File sampleBam=new File(args[0]);
			File control=new File(args[1]);
			File[] controls= {control};
			String save=args[2];
			String chr="chrX";
			
			new PlotSmoothedEnrichments(sampleBam, controls, 0.3, 100, 10, 10, save, chr);
		}
		else {System.err.println(usage);}
		
		
	}
	
	

	static String usage=" args[0]=sample bam \n args[1]=bin size \n args[2]=smooth window \n args[2]=save \n args[3]=chromosome to analyze";
}

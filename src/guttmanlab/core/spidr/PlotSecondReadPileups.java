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

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.ScanStat;
import guttmanlab.core.math.Statistics;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PlotSecondReadPileups {
	
	private static double percentile=0.5;
	private double alpha=0.0001;
	private double minEnrichment=3.0;
	private int minCount=10;
	Map<SingleInterval, Double> intraP;
	Map<SingleInterval, Double> interP;
	int numPerm=100;

	
	public PlotSecondReadPileups(String sampleFile, String inputFile, String save, Map<String, IntervalTree<Gene>> genes, int binSize) throws IOException{
		TreeMap<SingleInterval, Integer> sample=score(sampleFile, binSize);
		//write(sample, save+".sample.bedgraph");
		
		Map<SingleInterval, Double> intraScores=computeIntraPeaks(sample, genes, binSize);
		writeScores(intraScores, save+".intraScores");
		
		int totalSample=count(sampleFile);
		int totalInput=count(inputFile);
		
		Map<SingleInterval, Double> interScores=computeInterPeaks(sample, inputFile, genes, binSize, totalSample, totalInput, save);
		writeScores(interScores, save+".interScores");
		
		write(save, intraScores, interScores);
		
		Collection<SingleInterval> bins=getSignificantBins(intraScores, interScores, intraP, interP, sample);
		write(save+".sig.bed", bins);
	}
	
	
	public PlotSecondReadPileups(String sampleFile, String save) throws IOException{
		TreeMap<SingleInterval, Integer> sample=score(sampleFile, 1);
		write(sample, save);
	}
	
	
	public PlotSecondReadPileups(String sampleFile, File[] controlFiles, double p, String save) throws IOException{
		TreeMap<SingleInterval, Integer> sample=score(sampleFile, 1);
		TreeMap<SingleInterval, int[]> controls=score(controlFiles, 1, sample, p);
		write(sample, controls, save);
	}
	
	private TreeMap<SingleInterval, int[]> score(File[] controlFiles, int binSize, TreeMap<SingleInterval, Integer> sample, double p) {
		TreeMap<SingleInterval, int[]> rtrn=new TreeMap<SingleInterval, int[]>();
		
		for(SingleInterval region: sample.keySet()) {
			int[] array=new int[numPerm];
			rtrn.put(region, array);
		}
		
		
		for(int i=0; i<controlFiles.length; i++) {
			if(controlFiles[i].getName().endsWith(".bam")) {
				Map<SingleInterval, int[]> score=score(controlFiles[i], binSize, p);
				for(SingleInterval region: score.keySet()) {
					if(rtrn.containsKey(region)) {
						int[] list=rtrn.get(region);
						int[] list2=score.get(region);
						int[] sum=Statistics.sum(list, list2);
						rtrn.put(region, sum);
					}
				}
			}
		}
		return rtrn;
	}

	
	private boolean sample(double p) {
		return Math.random()<p;
	}

	private boolean[] sample(double p, int numPerm) {
		boolean[] rtrn=new boolean[numPerm];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=sample(p);
		}
		
		return rtrn;
	}
	

	private Map<SingleInterval, int[]> score(File bam, int binSize, double p) {
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<SingleInterval, int[]> positionCount=new TreeMap<SingleInterval, int[]>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			
			SAMFragment frag=new SAMFragment(read);
			
			//Should be the first nucleotide of Read1
			if(!read.getReadUnmappedFlag()){
				boolean passes=read.getFirstOfPairFlag();
				int startPos=read.getAlignmentStart();
				if(frag.getOrientation().equals(Strand.POSITIVE)) {
					passes=read.getSecondOfPairFlag();
					//startPos=read.getAlignmentStart();
					startPos=read.getAlignmentEnd();
				}
				
				if(passes) {
					boolean[] use=sample(p, numPerm);
					SingleInterval start=new SingleInterval(read.getReferenceName(), startPos, startPos+1);
					SingleInterval bin=start.bin(binSize);
					bin.setOrientation(getStrand(read));
					int[] count=new int[numPerm];
					if(positionCount.containsKey(bin)){count=positionCount.get(bin);}
					
					for(int i=0; i<use.length; i++) {
						if(use[i]) {count[i]++;}
					}
					
					positionCount.put(bin, count);
				}
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		return positionCount;
	}


	private Collection<SingleInterval> getSignificantBins(Map<SingleInterval, Double> intraScores,
			Map<SingleInterval, Double> interScores, Map<SingleInterval, Double> intraP2,
			Map<SingleInterval, Double> interP2, TreeMap<SingleInterval, Integer> sample) {
		
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		for(SingleInterval bin: intraScores.keySet()) {
			double intraEnrich=get(intraScores, bin, 0);
			double interEnrich=get(interScores, bin, 0);
			double intraP=get(intraP2, bin, 1.0);
			double interP=get(interP2, bin, 1.0);
			int count=get2(sample, bin, 0);
			if(intraP<alpha && interP<alpha && intraEnrich>this.minEnrichment && interEnrich>this.minEnrichment && count>this.minCount) {
				rtrn.add(bin);
			}
		}
		return rtrn;
	}


	private double get(Map<SingleInterval, Double> intraScores, SingleInterval bin, double i) {
		if(intraScores.containsKey(bin)) {return intraScores.get(bin);}
		return i;
	}
	
	private int get2(Map<SingleInterval, Integer> intraScores, SingleInterval bin, int i) {
		if(intraScores.containsKey(bin)) {return intraScores.get(bin);}
		return i;
	}


	private void write(String string, Collection<SingleInterval> bins) throws IOException {
		BEDFileIO.writeShortBED(bins, string);
	}


	private void writeScores(Map<SingleInterval, Double> interScores, String save) throws IOException {
		FileWriter writerPos=new FileWriter(save+".pos.bedgraph");
		FileWriter writerNeg=new FileWriter(save+".neg.bedgraph");
		
		for(SingleInterval r: interScores.keySet()) {
			double score=interScores.get(r);
			FileWriter writer=writerNeg;
			if(r.getOrientation().equals(Strand.POSITIVE)) {writer=writerPos;}
			writer.write(r.toBedgraph(score)+"\n");
		}
		
		
		writerPos.close();
		writerNeg.close();
	}


	public PlotSecondReadPileups(String sampleFile, String inputFile, String save, int binSize) throws IOException{
		TreeMap<SingleInterval, Integer> sample=score(sampleFile, binSize);
		write(sample, save+".sample.bedgraph");
		
		Map<SingleInterval, Double> intraScores=computeIntraPeaks(sample, binSize);
		
		int totalSample=count(sampleFile);
		int totalInput=count(inputFile);
		
		Map<SingleInterval, Double> interScores=computeInterPeaks(sample, inputFile, binSize, totalSample, totalInput, save);
		
		write(save, intraScores, interScores);
		
	}
	
	
	


	private void write(Map<SingleInterval, Integer> sample, String save) throws IOException {
		BEDFileIO.writeBEDGraphInteger(sample, save);
	}

	
	private void write(Map<SingleInterval, Integer> samples, Map<SingleInterval, int[]> controls, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: samples.keySet()) {
			int score=samples.get(region);
			int[] vals=controls.get(region);
			double quantile=Math.max(1, Statistics.quantile(vals, 0.95));
			double z=score-quantile; //Statistics.mean(vals);
			if(z<0) {z=0;}
			else {
				writer.write(region.toBedgraph(z)+"\n");
			}
		}
		
		writer.close();
	}

	private void write(String save, Map<SingleInterval, Double> intraScores, Map<SingleInterval, Double> interScores) throws IOException {
		FileWriter writerPos=new FileWriter(save+".pos.bedgraph");
		FileWriter writerNeg=new FileWriter(save+".neg.bedgraph");
		
		for(SingleInterval r: intraScores.keySet()) {
			if(interScores.containsKey(r)) {
				double score=Math.min(intraScores.get(r), interScores.get(r));
				//writer.write(r.toBED()+"\n");
				FileWriter writer=writerNeg;
				if(r.getOrientation().equals(Strand.POSITIVE)) {writer=writerPos;}
				writer.write(r.toBedgraph(score)+"\n");
			}
		}
		
		writerPos.close();
		writerNeg.close();
	}


	private Map<SingleInterval, Double> computeInterPeaks(Map<SingleInterval, Integer> sample, String inputFile, Map<String, IntervalTree<Gene>> genes, int binSize, int totalSample, int totalInput, String save) throws IOException {
		Map<SingleInterval, Integer> input=score(inputFile, binSize);
		Map<String, Integer> geneScore=scoreGenesInput(input, genes, percentile);
		input=update(input, sample, geneScore, genes);
		
		//go through each sample bin and compute significance relative to input
		Map<SingleInterval, Double> sigBins=computeSignificantBins(sample, input, totalSample, totalInput);
		return sigBins;
	}

	
	
	private Map<SingleInterval, Double> computeInterPeaks(Map<SingleInterval, Integer> sample, String inputFile, int binSize, int totalSample, int totalInput, String save) throws IOException {
		Map<SingleInterval, Integer> input=score(inputFile, binSize);
		//write(input, save+".input.bedgraph");
		Map<String, Integer> geneScore=scoreGenesInput(input, percentile);
		input=update(input, sample, geneScore);
		
		//go through each sample bin and compute significance relative to input
		Map<SingleInterval, Double> sigBins=computeSignificantBins(sample, input, totalSample, totalInput);
		return sigBins;
	}

	private Map<SingleInterval, Double> computeSignificantBins(Map<SingleInterval, Integer> sample, Map<SingleInterval, Integer> input, int totalSample, int totalInput) {
		this.interP=new TreeMap<SingleInterval, Double>();
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval bin: sample.keySet()) {
			int sampleWindow=sample.get(bin);
			int inputWindow=input.get(bin);
			
			double p=getPValue(sampleWindow, inputWindow, totalSample, totalInput);
			double enrich=getEnrichment(sampleWindow, inputWindow, totalSample, totalInput);
			rtrn.put(bin, enrich);
			this.interP.put(bin, p);
			//if(p<alpha && enrich>minEnrichment) {rtrn.put(bin, enrich);}
			
		}
		
		return rtrn;
	}
	
	private static double getEnrichment(int sampleWindow, int inputWindow, int sampleTotal, int inputTotal) {
		double num=((double)sampleWindow+1)/(double)sampleTotal;
		double denom=((double)inputWindow+1)/(double)inputTotal;
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


	private Map<SingleInterval, Integer> update(Map<SingleInterval, Integer> input, Map<SingleInterval, Integer> sample, Map<String, Integer> geneScore, Map<String, IntervalTree<Gene>> genes) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(SingleInterval bin: input.keySet()) {
			int val1=input.get(bin);
			Collection<String> overlappers=overlappingGenes(genes, bin);
			int max=max(overlappers, geneScore, val1);
			rtrn.put(bin, max);
		}
		
		for(SingleInterval bin: sample.keySet()) {
			if(!rtrn.containsKey(bin)) {
				Collection<String> overlappers=overlappingGenes(genes, bin);
				int max=max(overlappers, geneScore, 1);
				rtrn.put(bin, max);
			}
		}
		return rtrn;
	}
	
	private Map<SingleInterval, Integer> update(Map<SingleInterval, Integer> input, Map<SingleInterval, Integer> sample, Map<String, Integer> geneScore) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(SingleInterval bin: input.keySet()) {
			int val1=input.get(bin);
			int val2=0;
			if(geneScore.containsKey(bin.getReferenceName())) {
				val2=geneScore.get(bin.getReferenceName());
			}
			int max=Math.max(val1, val2);
			rtrn.put(bin, max);
		}
		
		for(SingleInterval bin: sample.keySet()) {
			if(!rtrn.containsKey(bin)) {
				int max=1;
				if(geneScore.containsKey(bin.getReferenceName())) {
					max=geneScore.get(bin.getReferenceName());
				}
				rtrn.put(bin, max);
			}
		}
		return rtrn;
	}

	
	
	private Map<SingleInterval, Double> computeIntraPeaks(Map<SingleInterval, Integer> sample, int binSize) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Map<String, Double> geneScore=scoreGenes(sample);
		
		System.err.println("done parsing length");
		
		int counter=0;
		for(SingleInterval pos: sample.keySet()){
			//Collection<String> overlappingGenes=overlappingGenes(genes, pos);
			double sampleScore=sample.get(pos);
			double expected=geneScore.get(pos.getReferenceName());
			double enrich=sampleScore/expected;
			
			counter++;
			
			if(expected>=0) {
				double pval=ScanStat.poissonPValue(sampleScore, expected);
				if(pval<alpha && enrich>minEnrichment) {
					rtrn.put(pos, enrich);
					//System.out.println(pos.toBedgraph(sampleScore));
				}
			}
			
			if(counter%100000==0) {System.err.println(counter+" "+sample.size());}
		}
		
		return rtrn;
	}

	private Map<SingleInterval, Double> computeIntraPeaks(Map<SingleInterval, Integer> sample,Map<String, IntervalTree<Gene>> genes, int binSize) {
		this.intraP=new TreeMap<SingleInterval, Double>();
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Map<String, Double> geneScore=scoreGenes(sample, genes);
		
		System.err.println("done parsing length");
		
		int counter=0;
		for(SingleInterval pos: sample.keySet()){
			Collection<String> overlappingGenes=overlappingGenes(genes, pos);
			double sampleScore=sample.get(pos);
			double expected=max(overlappingGenes, geneScore);
			double enrich=sampleScore/expected;
			
			/*if(expected>0) {
				System.out.println(pos.toBedgraph(enrich)+"\t"+pos.getOrientation());
			}*/
			
			counter++;
			
			if(expected>0) {
				rtrn.put(pos, enrich);
				
				double pval=ScanStat.poissonPValue(sampleScore, expected);
				this.intraP.put(pos, pval);
				/*if(pval<alpha && enrich>minEnrichment && sampleScore>minCount) {
					rtrn.put(pos, enrich);
					//System.out.println(pos.toBedgraph(sampleScore));
				}*/
			}
			
			if(counter%100000==0) {System.err.println(counter+" "+sample.size());}
		}
		
		return rtrn;
	}


	private double max(Collection<String> overlappingGenes, Map<String, Double> geneScore) {
		double max=-1;
		
		for(String g: overlappingGenes) {
			if(geneScore.containsKey(g)) {
				double score=geneScore.get(g);
				max=Math.max(max, score);
			}
		}
		
		return max;
	}
	
	private int max(Collection<String> overlappingGenes, Map<String, Integer> geneScore, int start) {
		int max=start;
		
		for(String g: overlappingGenes) {
			if(geneScore.containsKey(g)) {
				int score=geneScore.get(g);
				max=Math.max(max, score);
			}
		}
		
		return max;
	}

	
	private Map<String, Double> scoreGenes(Map<SingleInterval, Integer> sample) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		
		Map<String, Integer> length=new TreeMap<String, Integer>();
		Map<String, Integer> scores=new TreeMap<String, Integer>();
		
		int counter=0;
		for(SingleInterval r: sample.keySet()) {
			int binScore=sample.get(r);
			String g=r.getReferenceName();
				int numPositions=0;
				int score=0;
				if(length.containsKey(g)) {
					numPositions=length.get(g);
					score=scores.get(g);
				}
				numPositions++;
				score+=binScore;
				length.put(g, numPositions);
				scores.put(g, score);
			
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+sample.size()+" "+r.toUCSC());}
		}
		
		for(String g: scores.keySet()) {
			double score=scores.get(g);
			double len=length.get(g);
			double norm=score/len;
			rtrn.put(g, norm);
		}
		
		return rtrn;
	}

	private Map<String, Double> scoreGenes(Map<SingleInterval, Integer> sample, Map<String, IntervalTree<Gene>> genes) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		
		Map<String, Integer> length=new TreeMap<String, Integer>();
		Map<String, Integer> scores=new TreeMap<String, Integer>();
		
		int counter=0;
		for(SingleInterval r: sample.keySet()) {
			int binScore=sample.get(r);
			Collection<String> overlappingGenes=overlappingGenes(genes, r);
			for(String g: overlappingGenes) {
				int numPositions=0;
				int score=0;
				if(length.containsKey(g)) {
					numPositions=length.get(g);
					score=scores.get(g);
				}
				numPositions++;
				score+=binScore;
				length.put(g, numPositions);
				scores.put(g, score);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+sample.size()+" "+r.toUCSC());}
		}
		
		for(String g: scores.keySet()) {
			double score=scores.get(g);
			double len=length.get(g);
			double norm=score/len;
			rtrn.put(g, norm);
		}
		
		return rtrn;
	}
	
	private Map<String, Integer> scoreGenesInput(Map<SingleInterval, Integer> sample, Map<String, IntervalTree<Gene>> genes, double percentile) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		
		/*Map<String, Integer> length=new TreeMap<String, Integer>();
		Map<String, Integer> scores=new TreeMap<String, Integer>();*/
		
		
		
		Map<String, List<Integer>> scores=new TreeMap<String, List<Integer>>();
		
		int counter=0;
		for(SingleInterval r: sample.keySet()) {
			int binScore=sample.get(r);
			Collection<String> overlappingGenes=overlappingGenes(genes, r);
			for(String g: overlappingGenes) {
				if(!scores.containsKey(g)) {scores.put(g, new ArrayList<Integer>());}
				List<Integer> list=scores.get(g);
				list.add(binScore);
				scores.put(g, list);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+sample.size()+" "+r.toUCSC());}
		}
		
		for(String g: scores.keySet()) {
			List<Integer> list=scores.get(g);
			int p=Statistics.quantileInt(list, percentile);
			rtrn.put(g, p);
		}
		
		return rtrn;
	}
	
	
	private Map<String, Integer> scoreGenesInput(Map<SingleInterval, Integer> sample, double percentile) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		
		/*Map<String, Integer> length=new TreeMap<String, Integer>();
		Map<String, Integer> scores=new TreeMap<String, Integer>();*/
		
		
		
		Map<String, List<Integer>> scores=new TreeMap<String, List<Integer>>();
		
		int counter=0;
		for(SingleInterval r: sample.keySet()) {
			int binScore=sample.get(r);
			String g=r.getReferenceName();
				if(!scores.containsKey(g)) {scores.put(g, new ArrayList<Integer>());}
				List<Integer> list=scores.get(g);
				list.add(binScore);
				scores.put(g, list);
			
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+sample.size()+" "+r.toUCSC());}
		}
		
		for(String g: scores.keySet()) {
			List<Integer> list=scores.get(g);
			int p=Statistics.quantileInt(list, percentile);
			rtrn.put(g, p);
		}
		
		return rtrn;
	}
	
	
	private static Collection<String> overlappingGenes(Map<String, IntervalTree<Gene>> geneTree, SingleInterval bin) {
		Collection<String> rtrn=new TreeSet<String>();
		
		if(geneTree.containsKey(bin.getReferenceName())) {
			IntervalTree<Gene> tree=geneTree.get(bin.getReferenceName());
			Iterator<Gene> iter=tree.overlappingValueIterator(bin.getReferenceStartPosition(), bin.getReferenceEndPosition());
			while(iter.hasNext()) {
				Gene g=iter.next();
				if(g.getOrientation().equals(bin.getOrientation())) {
					if(bin.overlaps(g)) {
						rtrn.add(g.getName());
					}
				}
			}
		}
		return rtrn;
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
	
	private static TreeMap<SingleInterval, Integer> score(String bam, int binSize){
		SAMFileReader inputReader= new SAMFileReader(new File(bam));
		TreeMap<SingleInterval, Integer> positionCount=new TreeMap<SingleInterval, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			
			SAMFragment frag=new SAMFragment(read);
			
			
			//Should be the first nucleotide of Read1
			if(!read.getReadUnmappedFlag()){
				boolean passes=read.getFirstOfPairFlag();
				int startPos=read.getAlignmentStart();
				if(frag.getOrientation().equals(Strand.POSITIVE)) {
					passes=read.getSecondOfPairFlag();
					//startPos=read.getAlignmentStart();
					startPos=read.getAlignmentEnd();
				}
				
				if(passes) {
					SingleInterval start=new SingleInterval(read.getReferenceName(), startPos, startPos+1);
					SingleInterval bin=start.bin(binSize);
					bin.setOrientation(getStrand(read));
					int count=0;
					if(positionCount.containsKey(bin)){count=positionCount.get(bin);}
					count++;
					positionCount.put(bin, count);
				}
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		return positionCount;
	}
	
	
	
	private static int count(String bam){
		SAMFileReader inputReader= new SAMFileReader(new File(bam));
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getSecondOfPairFlag()){
				totalCount++;
			}
		}
				
				
		reads.close();
		inputReader.close();
		
		return totalCount;
	}
	
	

	

	
	
	public static void main(String[] args) throws IOException{
		if(args.length>2) {
			String sample=args[0];
			File[] controls=new File(args[1]).listFiles();
			String save=args[2];
			double p=Double.parseDouble(args[3]);
			new PlotSecondReadPileups(sample, controls, p, save);
		}
		else {System.err.println(usage);}
	}

	static String usage=" args[0]=sample bam \n args[1]=controls \n args[2]=save \n args[3]=fraction";
	
	
}

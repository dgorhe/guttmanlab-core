package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.ScanStat;
import guttmanlab.core.math.Statistics;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;

public class CLIPCLAPEnrichmentRewrite {
	
	double percentile=0.5;
	int windowSize;
	int stepSize;
	int inputTotalCounts;
	int sampleTotalCounts;
	int inputWindows;
	int sampleWindows;
	int totalGeneSizes;
	//Map<Gene, Integer> geneInputCounts;
	//Map<Gene, Integer> geneSampleCounts;
	
	
	public CLIPCLAPEnrichmentRewrite(File clipBAM, File inputBAM, int windowSize, String save, Map<String, IntervalTree<Gene>> regions, double percentile, int stepSize) throws IOException{
		this.stepSize=stepSize;
		this.windowSize=windowSize;
		this.percentile=percentile;
		//this.geneInputCounts=new TreeMap<Gene, Integer>();
		//this.geneSampleCounts=new TreeMap<Gene, Integer>();
		
		this.totalGeneSizes=getTotalGeneSizes(regions);
		
		Map<SingleInterval, CLAPScore> scores=computeEnrichment(clipBAM, inputBAM, regions, windowSize, stepSize);
		System.err.println("Done with scores");
		
		write(save, scores);
	}
	
	private int getTotalGeneSizes(Map<String, IntervalTree<Gene>> regions) {
		int sum=0;
		
		for(String chr: regions.keySet()){
			Iterator<Gene> iter=regions.get(chr).valueIterator();
			while(iter.hasNext()){
				sum+=iter.next().size();
			}
		}
		
		return sum;
	}

	private void write(String save, Map<SingleInterval, CLAPScore> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write(CLAPScore.getHeader()+"\n");
		
		System.err.println(scores.size());
		int count=0;
		for(SingleInterval region: scores.keySet()){
			CLAPScore score=scores.get(region);
			writer.write(score.toString()+"\n");
			count++;
			if(count%10000==0){System.err.println(region.toUCSC()+" "+count+" "+scores.size());}
		}
		
		
		
		writer.close();
	}

	
	
	private Map<SingleInterval, CLAPScore> computeEnrichment(File sampleBamFile, File inputBamFile, Map<String, IntervalTree<Gene>> geneTree, int windowSize, int stepSize) throws IOException {
		Map<SingleInterval, CLAPScore> rtrn=new TreeMap<SingleInterval, CLAPScore>();
		Map<SingleInterval, CLAPScore> finalRtrn=new TreeMap<SingleInterval, CLAPScore>();
		
		Map<SingleInterval, Integer> inputCounts=score(inputBamFile, geneTree, true);
		Map<SingleInterval, Integer> sampleCounts=score(sampleBamFile, geneTree, false);
		
		this.inputWindows=inputCounts.size();
		this.sampleWindows=sampleCounts.size();
		
		for(SingleInterval region: sampleCounts.keySet()){
			int s=getScore(sampleCounts, region);
			int i=getScore(inputCounts, region);
			CLAPScore score=new CLAPScore(region, s, i, this.sampleTotalCounts, this.inputTotalCounts, this.totalGeneSizes);
			rtrn.put(region, score);
		}
		
		//Add input quantile
		for(String chr: geneTree.keySet()){
			System.err.println(chr);
			Iterator<Gene> iter=geneTree.get(chr).valueIterator();
			while (iter.hasNext()){
				Gene gene=iter.next();
				//int geneInputCount=get(geneInputCounts,gene);
				//int geneSampleCount=get(geneSampleCounts, gene);
				
				Collection<SingleInterval> regions=getRegions(gene, windowSize);
				double quantile=getQuantile(regions, inputCounts);
				int inputSum=sum(regions, inputCounts);
				int sampleSum=sum(regions, sampleCounts);
				int val=new Double(quantile).intValue();
				for(SingleInterval region: regions){
					if(rtrn.containsKey(region)){
						CLAPScore score=rtrn.get(region);
						score.setInputPercentile(val);
						score.setGeneName(gene.getName());
						score.setNumberInputWindows(inputWindows);
						score.setNumberSampleWindows(sampleWindows);
						score.setGeneInputTotal(inputSum);
						score.setGeneSampleTotal(sampleSum);
						score.setGeneLength(gene.size());
						finalRtrn.put(region, score);
					}
				}
			}
		}
		
		return finalRtrn;
		
	}
	

	private int sum(Collection<SingleInterval> regions, Map<SingleInterval, Integer> inputCounts) {
		int sum=0;
		
		for(SingleInterval region: regions){
			if(inputCounts.containsKey(region)){sum+=inputCounts.get(region);}
		}
		
		return sum;
	}

	private int get(Map<Gene, Integer> geneInputCounts2, Gene gene) {
		int rtrn=0;
		
		if(geneInputCounts2.containsKey(gene)){return geneInputCounts2.get(gene);}
		
		return rtrn;
	}

	private Collection<SingleInterval> getRegions(Gene gene, int windowSize) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		int startIndex=gene.getReferenceStartPosition()/windowSize;
		int newStart=startIndex*windowSize;
		
		
		int newEnd=(gene.getReferenceEndPosition()/windowSize)*windowSize;
		
		for(int i=newStart; i<=newEnd; i+=windowSize){
			SingleInterval region=new SingleInterval(gene.getReferenceName(), i, i+windowSize, gene.getOrientation());
			if(gene.overlaps(region)){
				rtrn.add(region);
			}
		}
		
		return rtrn;
	}

	private double getQuantile(Collection<SingleInterval> regions, Map<SingleInterval, Integer> inputCounts) {
		List<Double> list=new ArrayList<Double>();
		for(SingleInterval region: regions){
			if(inputCounts.containsKey(region)){list.add(new Double(inputCounts.get(region)));}
		}
			
		double median=Statistics.quantile(list, 0.5);
			
		return median;
	}

	private Map<SingleInterval, Integer> score(File inputBamFile, Map<String, IntervalTree<Gene>> geneTree, boolean isInput) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		//Iterate over reads
		SAMFileReader inputReader= new SAMFileReader(inputBamFile);
				
		SAMRecordIterator reads=inputReader.iterator();
		int count=0;
		//Map<Gene, Integer> geneCount=this.geneSampleCounts;
		//if(isInput){geneCount=this.geneInputCounts;}
		
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			boolean overlaps=overlapsGene(read, geneTree);
			//Collection<Gene> overlappingGenes=getOverlappingGenes(read, geneTree);
			//if(overlappingGenes.size()>0){
			if(overlaps){	
				SingleInterval bin=bin(read);
				increment(rtrn, bin);
				//incrementGenes(overlappingGenes, geneCount);
			}
			count++;
			if(count%1000000 ==0){System.err.println(count);}
		}
				
		if(isInput){this.inputTotalCounts=count;}
		else{this.sampleTotalCounts=count;}
		reads.close();
		inputReader.close();
				
		//rtrn=normalize(rtrn, count);
	
		return rtrn;
	}

	private Collection<Gene> getOverlappingGenes(SAMRecord read, Map<String, IntervalTree<Gene>> geneTree) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		if(geneTree.containsKey(read.getReferenceName())){
			IntervalTree<Gene> genes=geneTree.get(read.getReferenceName());
			Iterator<Gene> overlappers=genes.overlappingValueIterator(read.getAlignmentStart(), read.getAlignmentEnd());
			if(!overlappers.hasNext()){return rtrn;}
			SAMFragment frag=new SAMFragment(read);
			while(overlappers.hasNext()){
				Gene gene=overlappers.next();
				if(gene.overlaps(frag)){rtrn.add(gene);}
			}
		}
		return rtrn;
	}

	private void incrementGenes(Collection<Gene> overlappingGenes, Map<Gene, Integer> geneCount) {
		for(Gene gene: overlappingGenes){
			int count=0;
			if(geneCount.containsKey(gene)){count=geneCount.get(gene);}
			count++;
			geneCount.put(gene, count);
		}
	}

	/*private Map<SingleInterval, Double> normalize(Map<SingleInterval, Double> rtrn, int count) {
		Map<SingleInterval, Double> normMap=new TreeMap<SingleInterval, Double>();
		for(SingleInterval region: rtrn.keySet()){
			double val=rtrn.get(region);
			double norm=10000*(val/(double)count);
			normMap.put(region, norm);
		}
		return normMap;
	}*/

	private boolean overlapsGene(SAMRecord read, Map<String, IntervalTree<Gene>> geneTree) {
		if(geneTree.containsKey(read.getReferenceName())){
			IntervalTree<Gene> genes=geneTree.get(read.getReferenceName());
			return genes.hasOverlappers(read.getAlignmentStart(), read.getAlignmentEnd());
		}
		return false;
	}

	private SingleInterval bin(SAMRecord read) {
		SAMFragment r=new SAMFragment(read);
		SingleInterval region=new SingleInterval(r.getReferenceName(), r.getReferenceStartPosition(), r.getReferenceEndPosition(), r.getOrientation());
		SingleInterval binned=region.bin(this.windowSize);
		binned.setOrientation(r.getOrientation());
		return binned;
	}

	private void increment(Map<SingleInterval, Integer> rtrn, SingleInterval bin) {
		int count=0;
		if(rtrn.containsKey(bin)){count=rtrn.get(bin);}
		count++;
		rtrn.put(bin, count);
	}

	
	private int getScore(Map<SingleInterval, Integer> geneInputCounts, Annotation region) {
		int score=0;
		if(geneInputCounts.containsKey(region)){score=geneInputCounts.get(region);}
		return score;
	}

	
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File clipBam=new File(args[0]);
			File inputBam=new File(args[1]);
			Map<String, IntervalTree<Gene>> genes=BEDFileIO.loadTree(args[2]);
			String save=args[3];
			int windowSize=100;
			int stepSize=100;
			double percentile=0.5;
			boolean useIntrons=new Boolean(args[4]);
			
			new CLIPCLAPEnrichmentRewrite(clipBam, inputBam, windowSize, save, genes, percentile, stepSize);
		
			if(useIntrons){
				Map<String, IntervalTree<Gene>> introns=BEDFileIO.loadIntronTree(args[2]);
				new CLIPCLAPEnrichmentRewrite(clipBam, inputBam, windowSize, save+".introns", introns, percentile, stepSize); 
			}
		
			
		
		
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clip BAM file \n arg[1]=input BAM file \n  args[2]=gene file \n args[3]=save \n args[4]=use introns";
}

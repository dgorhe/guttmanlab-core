package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.ScanStat;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class RandomDownsampleSignificant {
	double alpha=0.001;
	
	double totalReads;
	double genomeLength;
	double lambda;
	//int binSize=100;
	int extension=10;
	Map<SingleInterval, Double> localEnrichment;
	

	
	
	public RandomDownsampleSignificant(File bam, File[] controlBAMs, int numPerm, String save, double p, int[] binSizes) throws IOException {
		Map<SingleInterval, Integer> actualScores=getSignificantWindows(bam, binSizes);
		Map<SingleInterval, double[]> randomScores=getRandomScores(controlBAMs, p, numPerm, actualScores); //TODO we would need to extend here
		write(save, actualScores, randomScores, this.localEnrichment);
	}

	public RandomDownsampleSignificant(File bam, Map<File, Double> controlFiles, int numPerm, String save, int[] binSizes) throws IOException {
		Map<SingleInterval, Integer> actualScores=getSignificantWindows(bam, binSizes);
		Map<SingleInterval, double[]> randomScores=getRandomScores(controlFiles, numPerm, actualScores);
		write(save, actualScores, randomScores, this.localEnrichment);
	}
	
	
	private Map<SingleInterval, double[]> getMaxLocal(Map<SingleInterval, double[]> randomScores) {
		Map<SingleInterval, double[]> rtrn=new TreeMap<SingleInterval, double[]>();
		
		for(SingleInterval r: randomScores.keySet()) {
			Collection<SingleInterval> extensions=getExtensions(r);
			double[] avg=getAverage(extensions, randomScores, randomScores.get(r).length);
			double[] max=getMax(avg, randomScores.get(r));
			rtrn.put(r, max);
		}
		
		return rtrn;
	}

	private Collection<SingleInterval> getExtensions(SingleInterval r) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		int start=r.getReferenceStartPosition();
		String chr=r.getReferenceName();
		
		int minStart=start;
		int maxEnd=r.getReferenceEndPosition();
		int binSize=r.size();
		for(int i=-this.extension; i<this.extension; i++) {
			int newStart=start+(i*binSize);
			int newEnd=newStart+binSize;
			rtrn.add(new SingleInterval(chr, newStart, newEnd));
			minStart=Math.min(minStart, newStart);
			maxEnd=Math.max(maxEnd, newEnd);
		}
		
		//System.out.println(new SingleInterval(chr, minStart, maxEnd).toShortBED(r.toUCSC()));
		return rtrn;
	}

	private double[] getAverage(Collection<SingleInterval> extensions, Map<SingleInterval, double[]> randomScores, int length) {
		double[] sum=new double[length];
		
		int count=0;
		for(SingleInterval r: extensions) {
			sum=sum(sum, randomScores.get(r));
			count++;
		}
		
		return divide(sum, count);
		
	}

	private double[] getMax(double[] avg, double[] ds) {
		double avgMean=Statistics.mean(avg);
		double originalMean=Statistics.mean(ds);
		
		if(avgMean>originalMean) {return avg;}
		return ds;
	}

	private Map<SingleInterval, Integer> getSignificantWindows(File file, int[] binSizes) {
		Map<SingleInterval, Integer> counts=getScores(file, binSizes);
		this.localEnrichment=new TreeMap<SingleInterval, Double>();
		
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		this.lambda=this.totalReads/genomeLength;
		//int w=binSize;
		
		Map<Integer, Double> cutoffs=getCutoffs(binSizes);
		//double cutoff=getCutoff(w, lambda, genomeLength);
		
		
		for(SingleInterval region: counts.keySet()) {
			int score=counts.get(region);
			double cutoff=cutoffs.get(region.size());
			if(cutoff!=-1) {
				if(score>=cutoff) {
					rtrn.put(region, counts.get(region));
					double enrichment=localEnrichment(region, counts);
					this.localEnrichment.put(region, enrichment);
				}
			}
			else {
				System.err.println(lambda+" "+genomeLength+" "+region.size()+" "+cutoff);
				throw new IllegalStateException("no cutoff");
			}
		}
		
		
		return rtrn;
	}
	
	private Map<Integer, Double> getCutoffs(int[] binSizes) {
		Map<Integer, Double> rtrn=new TreeMap<Integer, Double>();
		
		for(int i=0; i<binSizes.length; i++) {
			double cutoff=getCutoff(binSizes[i], lambda, genomeLength);
			rtrn.put(binSizes[i], cutoff);
		}
		
		return rtrn;
	}

	private double localEnrichment(SingleInterval region, Map<SingleInterval, Integer> counts) {
		Collection<SingleInterval> extensions=this.getExtensions(region);
		
		extensions.remove(region);
		
		double[] vals=getVals(extensions, counts);
		
		
		double observed=counts.get(region);
		double expected=Statistics.mean(vals);
		double enrichment=observed/expected;
		return enrichment;
	}

	private SingleInterval collapse(Collection<SingleInterval> extensions) {
		return BEDFileIO.collapse(extensions).iterator().next();
	}

	private double[] getVals(Collection<SingleInterval> extensions, Map<SingleInterval, Integer> counts) {
		double[] vals=new double[extensions.size()];
		
		int i=0;
		for(SingleInterval r: extensions) {
			vals[i]=get(counts, r);
			i++;
		}
		
		
		return vals;
	}

	private double getCutoff(int w, double lambda2, double genomeLength2) {
		
		int start=(int) Math.rint(lambda2*w);
		start=Math.max(1, start);
		for(int i=start; i<10000*start; i++) {
			double pval=ScanStat.getPValue(i, lambda2, w, genomeLength2);
			System.err.println(w+" "+i+" "+pval);
			if(pval<alpha) {return i;}
		}
		
		return -1;
	}

	private Map<SingleInterval, Integer> getScores(File file, int[] binSizes) {
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		this.genomeLength=CoordinateSpace.getGenomeLength(reader.getFileHeader());
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			Collection<SingleInterval> allBins=SAMFragment.getSingleInterval(record).allBins(binSizes);	
			for(SingleInterval binned: allBins) {
				int score=0;
				if(counts.containsKey(binned)) {score=counts.get(binned);}
				score++;
				counts.put(binned, score);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		this.totalReads=counter;
		reader.close();
		reads.close();
		return counts;	
	}

	private double countReads(File file) {
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		int counter=0;
		while(reads.hasNext()){
			reads.next();
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
			
		reader.close();
		reads.close();
			
		return counter;
		
	}

	public RandomDownsampleSignificant(File bam, Collection<SingleInterval> intervals, String save, double alpha) throws IOException {
		Map<SingleInterval, Integer> actualScores=getScores(bam, intervals);
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: actualScores.keySet()) {
			double pval=ScanStat.getPValue(actualScores.get(r), lambda, r.size(), genomeLength);
			if(pval<alpha) {writer.write(r.toBedgraph(1.0)+"\n");}
		}
		
		writer.close();
	}
	
	

	private void write(String save, Map<SingleInterval, Integer> actualScores, Map<SingleInterval, double[]> randomScores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: actualScores.keySet()) {
			double[] rand=randomScores.get(r);
			double observed=actualScores.get(r);
			double p=Statistics.percentGreaterThan(observed, rand);
			double expected=Statistics.mean(rand);
			double enrichment=observed/expected;
			if(p<alpha && enrichment>2.0) {
				writer.write(r.toBedgraph(enrichment)+"\n");
			}
		}
		
		writer.close();
	}

	
	private void write(String save, Map<SingleInterval, Integer> actualScores, Map<SingleInterval, double[]> randomScores, Map<SingleInterval, Double> localScores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: actualScores.keySet()) {
			double[] rand=randomScores.get(r);
			double observed=actualScores.get(r);
			double p=Statistics.percentGreaterThan(observed, rand);
			double expected=Statistics.mean(rand);
			double enrichment=observed/expected;
			double localEnrich=localScores.get(r);
			double minEnrichment=Math.min(enrichment, localEnrich);
			if(p<alpha) {
				writer.write(r.toBedgraph(minEnrichment)+"\t"+enrichment+"\t"+localEnrich+"\n");
			}
		}
		
		writer.close();
	}







	private Map<SingleInterval, Integer> getScores(File bam, Collection<SingleInterval> regions) {
		Map<SingleInterval, Integer> rtrn=score(bam, regions);
		
		this.lambda=this.totalReads/genomeLength;
		
		
		return rtrn;
	}
	
	private Map<SingleInterval, Integer> score(File file, Collection<SingleInterval> regions) {
		Map<String, IntervalTree<SingleInterval>> tree=ChIPUtils.makeTree(regions);
		
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		this.genomeLength=CoordinateSpace.getGenomeLength(reader.getFileHeader());
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			Collection<SingleInterval> overlappers=ChIPUtils.getRegions(tree, record);
			for(SingleInterval o: overlappers) {
				increment(counts, o);
			}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		this.totalReads=counter;
		
		reader.close();
		reads.close();
		
		return counts;
	}
	
	private static void increment(Map<SingleInterval, Integer> counts, SingleInterval o) {
		int count=0;
		if(counts.containsKey(o)){
			count=counts.get(o);
		}
		count++;
		counts.put(o,  count);
	}


	private Map<SingleInterval, double[]> getRandomScores(File[] bamFiles, double p, int numPerm, Map<SingleInterval, Integer> actualScores) {
		Map<SingleInterval, double[]> scores=initialize(actualScores, numPerm);
		
		
		Map<String, IntervalTree<SingleInterval>> tree=ChIPUtils.makeTree(scores.keySet());
		
		int counter=0;
		for(int i=0; i<bamFiles.length; i++) {
			System.err.println(bamFiles[i].getName());
			SAMFileReader reader=new SAMFileReader(bamFiles[i]);
			
			SAMRecordIterator reads=reader.iterator();
			
			while(reads.hasNext()){
				boolean[] use=sample(p, numPerm);
				SAMRecord record=reads.next();
				
				Collection<SingleInterval> overlappers=ChIPUtils.getRegions(tree, record);
				
				for(SingleInterval binned: overlappers) {
					if(scores.containsKey(binned)) {
						double[] vals=scores.get(binned);
						vals=increment(vals, use);
						scores.put(binned, vals);
					}
				}
				
					
				counter++;
				if(counter%1000000 ==0){System.err.println(counter);}
			}
			
			reader.close();
			reads.close();
		}
		
		
		scores=getMaxLocal(scores);
		
		return scores;
	}
	
	private Map<SingleInterval, double[]> getRandomScores(File bamFile, double p, int numPerm, Map<SingleInterval, Integer> actualScores) {
		File[] files=new File[1];
		files[0]=bamFile;
		return getRandomScores(files, p, numPerm, actualScores);
	}
	
	
	private Map<SingleInterval, double[]> getRandomScores(Map<File, Double> controlFiles, int numPerm, Map<SingleInterval, Integer> actualScores) {
		File[] files=new File[controlFiles.size()];
		
		double sum=0;
		int index=0;
		for(File f: controlFiles.keySet()) {
			files[index]=f;
			sum+=controlFiles.get(f);
			index++;
		}
		
		double p=this.totalReads/sum;
		System.err.println(this.totalReads+" "+sum+" "+p);
		
		return getRandomScores(files, p, numPerm, actualScores);
	}
	
	private Map<SingleInterval, double[]> initialize(Map<SingleInterval, Integer> actualScores, int numPerm) {
		Map<SingleInterval, double[]> rtrn=new TreeMap<SingleInterval, double[]>();
		
		for(SingleInterval r: actualScores.keySet()) {
			Collection<SingleInterval> regions=getExtensions(r);
			for(SingleInterval bin: regions) {
				if(!rtrn.containsKey(bin)) {
					double[] vals=new double[numPerm];
					rtrn.put(bin, vals);
				}
			}
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


	/*private Map<SingleInterval, Integer>[] getRandomScores(File bamFile, double p, int binSize, int numPerm) {
		Map<SingleInterval, Integer>[] scores=new Map[numPerm];
		for(int i=0; i<scores.length; i++) {scores[i]=new TreeMap<SingleInterval, Integer>();}
		
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		
		SAMRecordIterator reads=reader.iterator();
		
		
		int counter=0;
		while(reads.hasNext()){
			boolean[] use=sample(p, numPerm);
			SAMRecord record=reads.next();
			SAMFragment f=new SAMFragment(record);
			SingleInterval binned=f.getSingleInterval().bin(binSize);
			add(scores, use, binned);
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		
		
		reader.close();
		reads.close();
		return scores;
	}*/
	
	
	private void add(Map<SingleInterval, double[]> scores, boolean[] use, SingleInterval binned) {
		double[] vals=new double[use.length];
		for(int i=0; i<use.length; i++) {
			if(use[i]) {vals[i]=1;}
		}
		
		if(scores.containsKey(binned)) {
			double[] other=scores.get(binned);
			vals=sum(vals, other);
		}
		
		scores.put(binned, vals);
		
		
	}
	
	private double[] sum(double[] vals, double[] other) {
		if(other==null) {return vals;}
		double[] sum=new double[vals.length];
		
		for(int i=0; i<vals.length; i++) {
			sum[i]=vals[i]+other[i];
		}
		
		return sum;
	}
	
	private double[] divide(double[] vals, double count) {
		double[] rtrn=new double[vals.length];
		
		for(int i=0; i<vals.length; i++) {
			rtrn[i]=vals[i]/count;
		}
		
		return rtrn;
	}


	/*private void add(Map<SingleInterval, Integer>[] scores, boolean[] use, SingleInterval binned) {
		for(int i=0; i<use.length; i++) {
			if(use[i]) {
				int score=0;
				if(scores[i].containsKey(binned)) {
					score=scores[i].get(binned);
				}
				score++;
				scores[i].put(binned, score);
			}
		}
	}*/


	private Map<SingleInterval, Integer> getRandomScores(File bamFile, double p, int binSize) {
		SAMFileReader reader=new SAMFileReader(bamFile);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		int collected=0;
		while(reads.hasNext()){
			boolean use=sample(p);
			SAMRecord record=reads.next();
				if(use) {
				
				SAMFragment f=new SAMFragment(record);
				SingleInterval binned=f.getSingleInterval().bin(binSize);
							
				int score=0;
				if(scores.containsKey(binned)) {
					score=scores.get(binned);
				}
				score++;
				scores.put(binned, score);
				collected++;
					
				
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter+" "+collected);}
		}
		
		
		
		reader.close();
		reads.close();
		return scores;
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
	

	private void write(String save, Map<SingleInterval, Integer>[] randomScores, Collection<SingleInterval> allIntervals) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: allIntervals) {
			writer.write(r.toUCSC());
			for(int i=0; i<randomScores.length; i++) {
				double score=get(randomScores[i], r);
				writer.write("\t"+score);
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	
	private double get(Map<SingleInterval, Integer> map, SingleInterval key) {
		if(map.containsKey(key)) {return map.get(key);}
		return 0;
	}
	
	private static Map<File, Double> parseFiles(String input) throws IOException {
		Map<File, Double> rtrn=new TreeMap<File, Double>();
		
		List<String> lines=BEDFileIO.loadLines(input);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			rtrn.put(new File(tokens[0]), Double.parseDouble(tokens[1]));
		}
		
		return rtrn;
	}
	
	private static File[] parseFiles2(String input) throws IOException {
		String[] tokens=input.split(",");
		File[] rtrn=new File[tokens.length];
		
		for(int i=0; i<tokens.length; i++) {
			rtrn[i]=new File(tokens[i]);
		}
		
		return rtrn;
	}
	
	private static int[] parseSizes(String input) throws IOException {
		String[] tokens=input.split(",");
		int[] rtrn=new int[tokens.length];
		
		for(int i=0; i<tokens.length; i++) {
			rtrn[i]=Integer.parseInt(tokens[i]);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		if(args.length>5) {
			File sampleBam=new File(args[0]);
			File[] controlFiles=parseFiles2(args[1]);
			double p=Double.parseDouble(args[2]);
			int numPerm=Integer.parseInt(args[3]);
			String save=args[4];
			int[] binSizes=parseSizes(args[5]);
			
			new RandomDownsampleSignificant(sampleBam, controlFiles, numPerm, save, p, binSizes);
		
		}
		else {System.err.println(usage);}
		
		
		
	}
	
	

	static String usage=" args[0]=sample bam \n args[1]=control bams (, sep) \n args[2]=fraction to sample \n args[3]=num perm \n args[4]=save \n args[5]=binSizes (, sep)";
}

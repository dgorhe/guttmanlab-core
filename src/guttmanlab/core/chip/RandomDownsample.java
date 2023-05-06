package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class RandomDownsample {

	public RandomDownsample(File bam, File controlBAM, double p, int numPerm, int binSize, String save, File sizesFile) throws IOException, InterruptedException {
		
		Map<SingleInterval, Integer> actualScores=getScores(bam, binSize);
		
		Map<SingleInterval, double[]> randomScores=getRandomScores(controlBAM, p, binSize, numPerm, actualScores);
		
		writeDiff(save+".subtract", actualScores, randomScores, sizesFile);
	}
	
	public RandomDownsample(File bam, File[] controlBAMs, double p, int numPerm, int binSize, String save, File sizesFile) throws IOException, InterruptedException {
		Map<SingleInterval, Integer> actualScores=getScores(bam, binSize);
		Map<SingleInterval, double[]> randomScores=getRandomScores(controlBAMs, p, binSize, numPerm, actualScores);
		writeDiff(save, actualScores, randomScores, sizesFile);
	}
	
	
	
	private void writeDiff(String save, Map<SingleInterval, Integer> actualScores, Map<SingleInterval, double[]> randomScores, File sizesFile) throws IOException, InterruptedException {
		FileWriter writer=new FileWriter(save);
		for(SingleInterval r: actualScores.keySet()) {
			double score=actualScores.get(r);
			double[] vals=randomScores.get(r);
			double z=score-Statistics.mean(vals);
			if(z<0) {z=0;}
			else {
				writer.write(r.toBedgraph(z)+"\n");
			}
		}
		writer.close();
		
		String cmd="/central/groups/guttman/software/kentUtils/bin/bedGraphToBigWig "+save+" "+sizesFile+" "+save+".bw";
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
	}


	


	private Map<SingleInterval, Integer> getScores(File bam, int binSize) {
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment f=new SAMFragment(record);
			//SingleInterval binned=f.getSingleInterval().bin(binSize);
			Collection<SingleInterval> allBins=f.getSingleInterval().allBins(binSize);
			
			for(SingleInterval binned: allBins) {
				int score=0;
				if(scores.containsKey(binned)) {
					score=scores.get(binned);
				}
				score++;
				scores.put(binned, score);
			}
			
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		
		
		reader.close();
		reads.close();
		return scores;
	}

	
	private Map<SingleInterval, double[]> getRandomScores(File bamFiles, double p, int binSize, int numPerm, Map<SingleInterval, Integer> actualScores) {
		File[] files=new File[1];
		files[0]=bamFiles;
		return getRandomScores(files, p, binSize,numPerm, actualScores);
	}

	private Map<SingleInterval, double[]> getRandomScores(File[] bamFiles, double p, int binSize, int numPerm, Map<SingleInterval, Integer> actualScores) {
		Map<SingleInterval, double[]> scores=initialize(actualScores, numPerm);
		
		//for(int i=0; i<scores.length; i++) {scores[i]=new TreeMap<SingleInterval, Integer>();}
		
		int counter=0;
		for(int i=0; i<bamFiles.length; i++) {
			SAMFileReader reader=new SAMFileReader(bamFiles[i]);
			SAMRecordIterator reads=reader.iterator();
			while(reads.hasNext()){
				boolean[] use=sample(p, numPerm);
				SAMRecord record=reads.next();
				//SingleInterval binned=SAMFragment.getSingleInterval(record).bin(binSize);
				Collection<SingleInterval> allbins=SAMFragment.getSingleInterval(record).allBins(binSize);
				for(SingleInterval binned: allbins) {
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
		return scores;
	}
	
	
	private Map<SingleInterval, double[]> initialize(Map<SingleInterval, Integer> actualScores, int numPerm) {
		Map<SingleInterval, double[]> rtrn=new TreeMap<SingleInterval, double[]>();
		
		for(SingleInterval r: actualScores.keySet()) {
			double[] vals=new double[numPerm];
			rtrn.put(r, vals);
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

	private boolean[] sample(double p, int numPerm) {
		boolean[] rtrn=new boolean[numPerm];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=sample(p);
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
	
	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length>5) {
		File sampleBam=new File(args[0]);
		File[] controlBam=parseFiles2(args[1]);
		double p=Double.parseDouble(args[2]);
		int numPerm=Integer.parseInt(args[3]);
		int binSize=10;
		String save=args[4];
		String sizesFile=args[5];
		
		new RandomDownsample(sampleBam, controlBam, p, numPerm, binSize, save, new File(sizesFile));
		
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]= sample bam \n args[1]=control bam \n args[2]=p \n args[3]=num perm \n args[4]=save \n args[5]=sizes file";
}

package guttmanlab.core.spidr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class RandomDownsample_Features {

	
	public RandomDownsample_Features(File bam, File[] controlBAMs, double p, int numPerm, String save) throws IOException, InterruptedException {
		Map<String, Integer> actualScores=getScores(bam);
		Map<String, double[]> randomScores=getRandomScores(controlBAMs, p, numPerm, actualScores);
		writeDiff(save+".diff.bedgraph", actualScores, randomScores);
	}
	
	private void writeActual(String save, Map<SingleInterval, Integer> actualScores) throws IOException, InterruptedException {
		FileWriter writer=new FileWriter(save);
		for(SingleInterval r: actualScores.keySet()) {
			double score=actualScores.get(r);
			writer.write(r.toBedgraph(score)+"\n");
		}
		writer.close();
	}
	
	
	private void writePVal(String save, Map<SingleInterval, Integer> actualScores, Map<SingleInterval, double[]> randomScores) throws IOException, InterruptedException {
		FileWriter writer=new FileWriter(save);
		for(SingleInterval r: actualScores.keySet()) {
			double score=actualScores.get(r);
			double[] vals=randomScores.get(r);
			double p=1-Statistics.percentLessThan(score, vals);
			writer.write(r.toBedgraph(p)+"\n");
		}
		writer.close();
	}
	
	private void writeDiff(String save, Map<String, Integer> actualScores, Map<String, double[]> randomScores) throws IOException, InterruptedException {
		FileWriter writer=new FileWriter(save);
		for(String r: actualScores.keySet()) {
			double score=actualScores.get(r);
			double[] vals=randomScores.get(r);
			double z=score-Statistics.quantile(vals, 0.95); //Statistics.mean(vals);
			if(z<0) {z=0;}
			else {
				writer.write(r+"\t"+z+"\n");
			}
		}
		writer.close();
	}
	

	private Map<String, Integer> getScores(File bam) {
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<String, Integer> scores=new TreeMap<String, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			int score=0;
			String chr=record.getReferenceName();
			if(scores.containsKey(chr)) {
				score=scores.get(chr);
			}
			score++;
			scores.put(chr, score);
			
			
			
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		
		
		reader.close();
		reads.close();
		return scores;
	}

	
	

	private Map<String, double[]> getRandomScores(File[] bamFiles, double p, int numPerm, Map<String, Integer> actualScores) {
		Map<String, double[]> scores=initialize(actualScores, numPerm);
		
		int counter=0;
		for(int i=0; i<bamFiles.length; i++) {
			SAMFileReader reader=new SAMFileReader(bamFiles[i]);
			SAMRecordIterator reads=reader.iterator();
			while(reads.hasNext()){
				boolean[] use=sample(p, numPerm);
				SAMRecord record=reads.next();
				String chr=record.getReferenceName();
				
				
				if(scores.containsKey(chr)) {
					double[] vals=scores.get(chr);
					vals=increment(vals, use);
					scores.put(chr, vals);
				}
				
					
				counter++;
				if(counter%1000000 ==0){System.err.println(counter);}
			}
			
			reader.close();
			reads.close();
		}
		return scores;
	}
	
	
	
	
	
	private Map<String, double[]> initialize(Map<String, Integer> actualScores, int numPerm) {
		Map<String, double[]> rtrn=new TreeMap<String, double[]>();
		
		for(String r: actualScores.keySet()) {
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
	
	
	
	
	private static double getP(File sampleBam, Map<File, Double> counts) {
		double sum=sum(counts);
		double num=counts.get(sampleBam);
		double denom=sum-num;
		return num/denom;
	}


	private static double sum(Map<File, Double> counts) {
		double rtrn=0;
		
		for(File f: counts.keySet()) {rtrn+=counts.get(f);}
		
		return rtrn;
	}


	private static double getCount(File[] controlBAMs) {
		double sum=0;
		
		for(int i=0; i<controlBAMs.length; i++) {
			sum+=getCount(controlBAMs[i]);
		}
		
		return sum;
	}
	

	/*private static Map<File, Double> getCount(File[] controlBAMs) {
		Map<File, Double> rtrn=new TreeMap<File, Double>();
		
		for(int i=0; i<controlBAMs.length; i++) {
			rtrn.put(controlBAMs[i], getCount(controlBAMs[i]));
		}
		
		return rtrn;
	}*/


	private static double getCount(File bam) {
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			counter++;
		}
		
		reader.close();
		reads.close();
		
		return counter;
	}


	private static File[] subset(File[] allBams, int index) {
		File[] rtrn=new File[allBams.length-1];
		
		int counter=0;
		for(int i=0; i<allBams.length; i++) {
			if(i!=index) {
				rtrn[counter]=allBams[i];
				counter++;
			}
		}
		
		return rtrn;
	}
	

	private static File[] getFiles(String[] args) {
		File[] rtrn=new File[args.length-2];
		for(int i=1; i<args.length-1; i++) {
			rtrn[i-1]=new File(args[i]);
		}
		return rtrn;
	}

	
	private static double getP(File sampleBam, File[] controlBAMs) {
		double num=getCount(sampleBam);
		double denom=getCount(controlBAMs);
		return num/denom;
	}


	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length>1) {
			int numPerm=100;
			
			File sampleBam=new File(args[0]);
			File[] controlBAMs=getFiles(args);
			double p=getP(sampleBam, controlBAMs);
			System.err.println(sampleBam.getName()+"\t"+p);
			String save=args[args.length-1];
			new RandomDownsample_Features(sampleBam, controlBAMs, p, numPerm, save);
			System.err.println("done "+sampleBam.getAbsolutePath());
			
			
			
		}
		else {System.err.println(usage);}
	}



	static String usage=" args[0]= sample bam \n args[1,...n-1]=control bams \n args[n]=save";
}

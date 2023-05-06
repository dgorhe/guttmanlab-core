package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class RandomDownsample_CLIP {

	
	public RandomDownsample_CLIP(File bam, File[] controlBAMs, double p, int numPerm, int binSize, String save, String chr) throws IOException, InterruptedException {
		Map<SingleInterval, Integer> actualScores=getScores(bam, binSize, chr);
		Map<SingleInterval, double[]> randomScores=getRandomScores(controlBAMs, p, binSize, numPerm, actualScores, chr);
		writeDiff(save+".diff.bedgraph", actualScores, randomScores);
		writeRandom(save+".random.bedgraph", actualScores, randomScores);
		writeActual(save+".actual.bedgraph", actualScores);
		writePVal(save+".pval.bed", actualScores, randomScores);
	}
	
	private void writeRandom(String save, Map<SingleInterval, Integer> actualScores, Map<SingleInterval, double[]> randomScores) throws IOException, InterruptedException {
		FileWriter writer=new FileWriter(save);
		for(SingleInterval r: randomScores.keySet()) {
			//double score=actualScores.get(r);
			double[] vals=randomScores.get(r);
			double quantile=Statistics.quantile(vals, 0.95); //Statistics.mean(vals);
			
			writer.write(r.toBedgraph(quantile)+"\t"+Statistics.quantile(vals, 0.05)+"\t"+Statistics.quantile(vals, 0.5)+"\t"+Statistics.mean(vals));
			
			for(int i=0; i<vals.length; i++) {
				writer.write("\t"+vals[i]);
			}
			writer.write("\n");
			
		}
		writer.close();
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
	
	private void writeDiff(String save, Map<SingleInterval, Integer> actualScores, Map<SingleInterval, double[]> randomScores) throws IOException, InterruptedException {
		FileWriter writer=new FileWriter(save);
		for(SingleInterval r: actualScores.keySet()) {
			double score=actualScores.get(r);
			double[] vals=randomScores.get(r);
			double z=score-Statistics.quantile(vals, 0.95); //Statistics.mean(vals);
			if(z<0) {z=0;}
			if(z>0) {
			//else {
				if(r.getOrientation().equals(Strand.NEGATIVE)) {z=-z;}
				writer.write(r.toBedgraph(z)+"\n");
			}
		}
		writer.close();
	}
	

	private Map<SingleInterval, Integer> getScores(File bam, int binSize, String chr) {
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			if(sameChr(record, chr)) {
			
				SAMFragment f=new SAMFragment(record);
				//SingleInterval binned=f.getSingleInterval().bin(binSize);
				Collection<SingleInterval> allBins=f.getSingleInterval().allBins(binSize);
				
				for(SingleInterval binned: allBins) {
					binned.setOrientation(f.getOrientation());
					int score=0;
					if(scores.containsKey(binned)) {
						score=scores.get(binned);
					}
					score++;
					scores.put(binned, score);
				}
			}
				
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		
		
		
		reader.close();
		reads.close();
		return scores;
	}

	
	private boolean sameChr(SAMRecord record, String chr) {
		if(chr.equalsIgnoreCase("all")) {return true;}
		return record.getReferenceName().equals(chr);
	}

	private Map<SingleInterval, double[]> getRandomScores(File bamFiles, double p, int binSize, int numPerm, Map<SingleInterval, Integer> actualScores, String chr) {
		File[] files=new File[1];
		files[0]=bamFiles;
		return getRandomScores(files, p, binSize,numPerm, actualScores, chr);
	}

	private Map<SingleInterval, double[]> getRandomScores(File[] bamFiles, double p, int binSize, int numPerm, Map<SingleInterval, Integer> actualScores, String chr) {
		//Map<SingleInterval, double[]> scores=initialize(actualScores, numPerm);
		
		Map<SingleInterval, double[]> scores=new TreeMap<SingleInterval, double[]>();
		
		int counter=0;
		for(int i=0; i<bamFiles.length; i++) {
			System.err.println(bamFiles[i].getAbsolutePath());
			SAMFileReader reader=new SAMFileReader(bamFiles[i]);
			SAMRecordIterator reads=reader.iterator();
			while(reads.hasNext()){
				SAMRecord record=reads.next();
				if(sameChr(record, chr)) {
					boolean[] use=sample(p, numPerm);
					
					SAMFragment f=new SAMFragment(record);
					Collection<SingleInterval> allbins=f.getSingleInterval().allBins(binSize);
					
					for(SingleInterval binned: allbins) {
						binned.setOrientation(f.getOrientation()); //TODO
						//binned.setOrientation(SAMFragment.getOrientation(record));
						if(!scores.containsKey(binned)) {scores.put(binned, new double[use.length]);}
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
	
	
	
	private Map<SingleInterval, double[]> getRandomScores(File[] bamFiles, double p, int numPerm, Map<SingleInterval, Integer> actualScores) {
		Map<SingleInterval, double[]> scores=initialize(actualScores, numPerm);
		
		int counter=0;
		for(int i=0; i<bamFiles.length; i++) {
			SAMFileReader reader=new SAMFileReader(bamFiles[i]);
			SAMRecordIterator reads=reader.iterator();
			while(reads.hasNext()){
				boolean[] use=sample(p, numPerm);
				SAMRecord read=reads.next();
				
				SingleInterval start=new SingleInterval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentStart()+1);
					if(scores.containsKey(start)) {
						double[] vals=scores.get(start);
						vals=increment(vals, use);
						scores.put(start, vals);
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
		
		//System.err.println(vals.length+" "+use.length);
		
		for(int i=0; i<vals.length; i++) {
			if(use[i]) {
				//System.err.println("increment");
				rtrn[i]=vals[i]+1;
			}
		}
		
		return rtrn;
	}


	
	

	private boolean sample(double p) {
		double rand=Math.random();
		boolean lessThan=rand<p;
		//System.err.println(rand+" "+p+" "+lessThan);
		return lessThan;
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
	
	private static File[] getFiles(File[] files, File sampleBam) {
		List<File> rtrn=new ArrayList<File>();
		for(int i=0; i<files.length; i++) {
			if(!files[i].equals(sampleBam)) {
				if(files[i].getName().endsWith("bam")) {rtrn.add(files[i]);}
			}
		}
		
		File[] list=new File[rtrn.size()];
		for(int i=0; i<rtrn.size(); i++) {
			list[i]=rtrn.get(i);
		}
		
		return list;
	}

	
	private static double getP(File sampleBam, File[] controlBAMs) {
		double num=getCount(sampleBam);
		double denom=getCount(controlBAMs);
		return num/denom;
	}


	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length>5) {
			
			File sampleBam=new File(args[0]);
			File[] controlBAMs=getFiles(new File(args[1]).listFiles(), sampleBam);
			int binSize=Integer.parseInt(args[2]);
			String save=args[3];
			String chr=args[4];
			
			int numPerm=100;
			//double p=getP(sampleBam, controlBAMs);
			double p=Double.parseDouble(args[5]);
			System.err.println(sampleBam.getName()+"\t"+p);
			
			new RandomDownsample_CLIP(sampleBam, controlBAMs, p, numPerm, binSize, save, chr);
			System.err.println("done "+sampleBam.getAbsolutePath());
			
			
			
		}
		else {System.err.println(usage);}
	}



	static String usage=" args[0]= sample bam \n args[1]=directory of all bams \n args[2]=bin size \n args[3]=save \n args[4]=chr \n args[5]=fraction";
}

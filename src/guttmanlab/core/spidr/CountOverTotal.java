package guttmanlab.core.spidr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class CountOverTotal {

	
	public CountOverTotal(File[] bams, String save) throws IOException, InterruptedException {
		
		Map<String, Double>[] scores=new Map[bams.length];
		
		for(int i=0; i<bams.length; i++) {
			System.err.println(bams[i]+" "+i+" "+bams.length);
			scores[i]=getScores(bams[i]);
		}
		
		MatrixWithHeaders matrix= makeMatrix(scores, bams);
		matrix.write(save);
		
	}
	
	private MatrixWithHeaders makeMatrix(Map<String, Double>[] scores, File[] bams) {
		List<String> columns=new ArrayList<String>();
		Collection<String> rows=new TreeSet<String>();
		for(int i=0; i<bams.length; i++) {
			columns.add(bams[i].getName());
			rows.addAll(scores[i].keySet());
		}
		
		List<String> rowList=new ArrayList<String>();
		rowList.addAll(rows);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rowList, columns);
		
		for(int i=0; i<scores.length; i++) {
			String column=bams[i].getName();
			for(String row: scores[i].keySet()) {
				double score=scores[i].get(row);
				rtrn.set(row, column, score);
			}
		}
		
		return rtrn;
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
	

	private Map<String, Double> getScores(File bam) {
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
		
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		for(String chr: scores.keySet()) {
			double score=scores.get(chr);
			double norm=1000000*(score/(double)counter);
			rtrn.put(chr, norm);
		}
		
		return rtrn;
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

	
	private static File[] getBams(File[] listFiles) {
		List<File> files=new ArrayList<File>();
		
		
		for(int i=0; i<listFiles.length; i++) {
			if(listFiles[i].getName().endsWith(".bam")) {files.add(listFiles[i]);}
		}
		
		
		File[] rtrn=new File[files.size()];
		
		for(int i=0; i<files.size(); i++) {rtrn[i]=files.get(i);}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length>1) {
			File[] sampleBams=getBams(new File(args[0]).listFiles());
			String save=args[1];
			new CountOverTotal(sampleBams, save);;
			
			
			
		}
		else {System.err.println(usage);}
	}



	



	static String usage=" args[0]= sample bams \n args[1]=save";
}

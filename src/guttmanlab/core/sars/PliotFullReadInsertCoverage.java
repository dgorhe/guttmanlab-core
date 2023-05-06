package guttmanlab.core.sars;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PliotFullReadInsertCoverage {

	int totalCounts;
	int inputCounts;
	
	public PliotFullReadInsertCoverage(File bam, int binSize, String save) throws IOException{
		Map<SingleInterval, Integer> scores=score(bam, binSize);
		
		write(save, scores);
		
	}
	
	public PliotFullReadInsertCoverage(File bam, File input, int binSize, String save) throws IOException{
		Map<SingleInterval, Integer> scores=score(bam, binSize, false);
		Map<SingleInterval, Integer> inputScores=score(input, binSize, true);
		
		write(save, scores, inputScores);
		
	}
	
	private void write(String save, Map<SingleInterval, Integer> scores, Map<SingleInterval, Integer> inputScores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: scores.keySet()){
			double score=((double)scores.get(region)/(double)totalCounts);
			double inputscore=((double)inputScores.get(region)/(double)inputCounts);
			double ratio=score/inputscore;
			writer.write(region.tobedgraph(ratio)+"\n");
		}
		
		writer.close();
	}
	
	private void write(String save, Map<SingleInterval, Integer> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: scores.keySet()){
			double score=1000000*((double)scores.get(region)/(double)totalCounts);
			writer.write(region.tobedgraph(score)+"\n");
		}
		
		writer.close();
	}
	
	private Collection<SingleInterval> bin(SAMRecord read, int binSize) {
		//get read and mate and fill in all bins in between
		
		String chr=read.getReferenceName();
		int mateStart=read.getMateAlignmentStart();
		int readStart=read.getAlignmentStart();
		
		
		Collection<SingleInterval> rtrn=binnedInterval(chr, mateStart, readStart, binSize);
		return rtrn;
	}
	
	
	

	private Collection<SingleInterval> binnedInterval(String chr, int start1, int start2, int binSize) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		int startIndex=Math.min(start1, start2)/binSize;
		int endIndex=Math.max(start1, start2)/binSize;
				
		for(int i=startIndex; i<=endIndex; i++){
			int newStart=i*binSize;
			int newEnd=newStart+binSize;
			SingleInterval newInterval=new SingleInterval(chr, newStart, newEnd);
			rtrn.add(newInterval);
			//System.err.println(start1+" "+start2+" "+newInterval.toUCSC());
		}
		
		return rtrn;
	}
	
	private TreeMap<SingleInterval, Integer> score(File bam, int binSize){
		return score(bam, binSize, false);
	}
	
	private TreeMap<SingleInterval, Integer> score(File bam, int binSize, boolean input){
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<SingleInterval, Integer> positionCount=new TreeMap<SingleInterval, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1 && read.getFirstOfPairFlag() && !read.getMateUnmappedFlag() && read.getReferenceName().equals(read.getMateReferenceName())){
				Collection<SingleInterval> bins=bin(read, binSize);
				
				for(SingleInterval bin: bins){
					int score=0;
					if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
					score++;
					positionCount.put(bin, score);
				}
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		if(input){this.inputCounts=totalCount;}
		else{this.totalCounts=totalCount;}
		
		return positionCount;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File file=new File(args[0]);
			//File input=new File(args[1]);
			String save=args[1];
			int bin=new Integer(args[2]);
			
			
			new PliotFullReadInsertCoverage(file, bin, save);
		}
		else{System.err.println(usage);}
	}
	
	
	static String usage=" args[0]=bam file \n args[1]=save \n args[2]=binSize";
	
}

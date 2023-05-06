package guttmanlab.core.sars;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;

public class BinCounts {

	int totalCounts;
	public BinCounts(File file, String save, int bin) throws IOException{
		Map<String, Integer> scores=score(file, bin);
		long totalSize=getSize(file);
		write(save, scores, bin, totalSize);
	}
	
	private long getSize(File file) {
		SAMFileReader inputReader= new SAMFileReader(file);
		
		SAMSequenceDictionary dict=inputReader.getFileHeader().getSequenceDictionary();
					
		inputReader.close();
		return dict.getReferenceLength();
		
	}
	
	private void write(String save, Map<String, Integer> scores, int bin, long totalSize) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String region: scores.keySet()){
			writer.write(region+"\t"+scores.get(region)+"\t"+totalCounts+"\t"+bin+"\t"+totalSize+"\n");
		}
		
		writer.close();
	}
	
	private String bin(SAMRecord read, int binSize) {
		int startIndex=read.getAlignmentStart()/binSize;
		int newStart=startIndex*binSize;
		int newEnd=newStart+binSize;
		SingleInterval newInterval=new SingleInterval(read.getReferenceName(), newStart, newEnd);
		return newInterval.toUCSC();
	}

	private TreeMap<String, Integer> score(File bam, int binSize){
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<String, Integer> positionCount=new TreeMap<String, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				String bin=bin(read, binSize);
				
				int score=0;
				if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
				score++;
				positionCount.put(bin, score);
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		this.totalCounts=totalCount;
		
		return positionCount;
	}




	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File file=new File(args[0]);
			String save=args[1];
			int bin=new Integer(args[2]);
			
			
			new BinCounts(file, save, bin);
		}
		else{System.err.println(usage);}
	}
	
	
	static String usage=" args[0]=bam file \n args[1]=save \n args[2]=binSize";
}

package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.rnasprite.Cluster;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class XToA {

	public XToA(File bamFile,File inputFile, String save, int binResolution) throws IOException{
		//TODO Compute per Kb of sequence
		
		Map<SingleInterval, Double> CHIPRatio=getRatio(bamFile, binResolution);
		Map<SingleInterval, Double> InputRatio=getRatio(inputFile, binResolution);
		
		
		write(save+"."+binResolution+".input.bedgraph", InputRatio);
		write(save+"."+binResolution+".ChIP.bedgraph", CHIPRatio);
		write(save+"."+binResolution+".ratio.bedgraph", CHIPRatio, InputRatio);
		
	}
	
	
	private Map<SingleInterval, Double> getRatio(File bamFile, int binResolution) {
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Integer> regionCount=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(!record.getNotPrimaryAlignmentFlag()){
				SingleInterval region=bin(record, binResolution);
				
				int count=0;
				if(regionCount.containsKey(region)){count=regionCount.get(region);}
				count++;
				regionCount.put(region, count);
				
				
				
				counter++;
				if(counter%1000000 ==0){System.err.println(counter);}
			}
		}
		
		
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		for(SingleInterval region: regionCount.keySet()){
			double score=regionCount.get(region);
			double ratio=1000000 *(score/(double)counter);
			rtrn.put(region, ratio);
		}
		
		
		reads.close();
		reader.close();
		
		return rtrn;
	}


	private SingleInterval bin(SAMRecord record, int resolution) {
		int startIndex=record.getAlignmentStart()/resolution;
		int newStart=startIndex*resolution;
		int newEnd=newStart+Math.max(record.getReadLength(), resolution);
		SingleInterval newInterval=new SingleInterval(record.getReferenceName(), newStart, newEnd);
		return newInterval;
	}


	private void write(String string, Map<SingleInterval, Double> regionCount) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(SingleInterval region: regionCount.keySet()){
			writer.write(region.toBedgraph(regionCount.get(region))+"\n");
		}
		
		writer.close();
	}
	
	private void write(String string, Map<SingleInterval, Double> regionCount, Map<SingleInterval, Double> denomCounts) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(SingleInterval region: regionCount.keySet()){
			double num=regionCount.get(region);
			if(denomCounts.containsKey(region)){
				double denom=denomCounts.get(region);
				double ratio=num/denom;
				writer.write(region.toBedgraph(ratio)+"\n");
			}
		}
		
		writer.close();
	}


	private void write(String save, Map<String, Integer> chrCounts, SAMSequenceDictionary samSequenceDictionary) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<String, Integer> chrLength=getLength(samSequenceDictionary);
		
		for(String chr: chrCounts.keySet()){
			double count=chrCounts.get(chr);
			double length=chrLength.get(chr);
			double ratio=count/length;
			writer.write(chr+"\t"+count+"\t"+length+"\t"+ratio+"\n");
		}
		
		
		
		writer.close();
	}


	private Map<String, Integer> getLength(SAMSequenceDictionary samSequenceDictionary) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(SAMSequenceRecord r: samSequenceDictionary.getSequences()){
			String chr=r.getSequenceName();
			int length=r.getSequenceLength();
			rtrn.put(chr,  length);
			System.err.println(chr+" "+length);
		}
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException{
		File bamFile=new File(args[0]);
		File inputFile=new File(args[1]);
		String save=args[2];
		int binResolution=new Integer(args[3]);
		new XToA(bamFile, inputFile, save, binResolution);
	}
	
}

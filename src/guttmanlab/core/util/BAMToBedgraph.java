package guttmanlab.core.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BAMToBedgraph {

	
	public static void convertBAMToBedgraphPairs(File bam, String save, int binSize) throws IOException {
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator reads=inputReader.iterator();
		
		Map<String, Pair<SAMRecord>> pairedReads=new TreeMap<String, Pair<SAMRecord>>();
		
		
		int totalCount=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			
			String name=read.getReadName();
			Pair<SAMRecord> pair=new Pair<SAMRecord>();
			if(pairedReads.containsKey(name)) {pair=pairedReads.get(name);}
			pair=update(pair, read);
			if(complete(pair)) {
				pairedReads.remove(name);
				add(pair, counts, binSize);
			}
			else {
				pairedReads.put(name, pair);
			}
			
			
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		BEDFileIO.writeBEDGraphInteger(counts, save);
		
	}
	
	
	public static void convertBAMToBedgraphSingle(File bam, String save, int binSize) throws IOException {
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator reads=inputReader.iterator();
		
		
		
		int totalCount=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			String name=read.getReadName();
			add(read, counts, binSize);
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		BEDFileIO.writeBEDGraphInteger(counts, save);
		
	}
	
	
	private static void add(SAMRecord read, Map<SingleInterval, Integer> counts, int binSize) {
		
		Collection<SingleInterval> bins=SAMFragment.allBins(read, binSize);
		
		for(SingleInterval bin: bins) {
			int count=0;
			if(counts.containsKey(bin)) {count=counts.get(bin);}
			count++;
			counts.put(bin, count);
		}
	}
	
	private static void add(Pair<SAMRecord> pair, Map<SingleInterval, Integer> counts, int binSize) {
		SAMFragment read1=new SAMFragment(pair.getValue1());
		SAMFragment read2=new SAMFragment(pair.getValue2());
		
		
		SingleInterval fragment=new SingleInterval(read1.getReferenceName(), Math.min(read1.getReferenceStartPosition(), read2.getReferenceStartPosition()), Math.max(read1.getReferenceEndPosition(), read2.getReferenceEndPosition()));
		fragment.setOrientation(read1.getOrientation());
		
		Collection<SingleInterval> bins= fragment.allBins(binSize);
		
		for(SingleInterval bin: bins) {
			int count=0;
			if(counts.containsKey(bin)) {count=counts.get(bin);}
			count++;
			counts.put(bin, count);
		}
	}

	private static boolean complete(Pair<SAMRecord> pair) {
		return pair.isComplete();
	}
	
	private static Pair<SAMRecord> update(Pair<SAMRecord> pair, SAMRecord read) {
		if(read.getSecondOfPairFlag()) {pair.setValue2(read);}
		if(read.getFirstOfPairFlag()) {pair.setValue1(read);}
		return pair;
	}
	
	public static void main(String[] args) throws IOException {
		File bam=new File(args[0]);
		String save=args[1];
		int binSize=Integer.parseInt(args[2]);
		convertBAMToBedgraphSingle(bam, save, binSize);
	}
	
}

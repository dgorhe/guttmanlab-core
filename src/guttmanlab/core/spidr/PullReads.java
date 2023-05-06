package guttmanlab.core.spidr;

import java.io.File;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.Pair;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PullReads {

	public static void pullReads(File bam, String save) {
		
		SAMFileReader reader=new SAMFileReader(bam);
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save));
		
		
		SAMRecordIterator reads=reader.iterator();
		
		int total=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(record.getReferenceName().equals("18S")) {
				if(record.getProperPairFlag() && Math.abs(record.getInferredInsertSize())>70 && record.getMappingQuality()>10) {
					writer.addAlignment(record);
				}
			}
				
			total++;
			if(total%1000000 ==0){System.err.println(total);}
		}
		
		
		
		reader.close();
		reads.close();
		writer.close();
	}
	
	public static void plotInsert(File bam, int binSize) {
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<String, Pair<SAMRecord>> pairs=new TreeMap<String, Pair<SAMRecord>>();
		
		int total=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(record.getReferenceName().equals("18S")) {
				if(record.getProperPairFlag() && Math.abs(record.getInferredInsertSize())>70 && record.getMappingQuality()>10) {
					String name=record.getReadName();
					if(!pairs.containsKey(name)) {pairs.put(name, new Pair<SAMRecord>());}
					Pair<SAMRecord> pair=pairs.get(name);
					if(record.getFirstOfPairFlag()) {pair.setValue1(record);}
					else {pair.setValue2(record);}
					pairs.put(name, pair);
				}
			}
				
			total++;
			if(total%1000000 ==0){System.err.println(total);}
		}
		
		reader.close();
		reads.close();
		
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		for(String name: pairs.keySet()) {
			Pair<SAMRecord> pair=pairs.get(name);
			if(pair.isComplete()) {
				SingleInterval region=getFragment(pair);
				//System.out.println(region);
				Collection<SingleInterval> bins= region.allBins(binSize);
				for(SingleInterval bin: bins) {
					int score=0;
					if(scores.containsKey(bin)) {score=scores.get(bin);}
					score++;
					scores.put(bin, score);
				}
			}
		}
		
		write(scores);
		
	}
	
	private static void write(Map<SingleInterval, Integer> scores) {
		for(SingleInterval bin: scores.keySet()) {
			System.out.println(bin.toBedgraph(scores.get(bin)));
		}
		
	}

	private static SingleInterval getFragment(Pair<SAMRecord> pair) {
		String chr=pair.getValue1().getReferenceName();
		int start=Math.min(pair.getValue1().getAlignmentStart(), pair.getValue2().getAlignmentStart());
		int end=Math.max(pair.getValue1().getAlignmentEnd(), pair.getValue2().getAlignmentEnd());
		return new SingleInterval(chr, start, end);
	}

	public static void main(String[] args) {
		File bam=new File(args[0]);
		//String save=args[1];
		plotInsert(bam, 10);
	}
	
}

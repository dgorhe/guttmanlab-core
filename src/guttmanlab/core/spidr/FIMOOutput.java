package guttmanlab.core.spidr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class FIMOOutput {

	public FIMOOutput(File bam, Collection<SingleInterval> regions, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		Map<SingleInterval, Integer> counts=score(bam);
		
		for(SingleInterval region: regions) {
			int[] scores=getScore(region, counts);
			writer.write(region.toUCSCStrand());
			for(int i=0; i<scores.length; i++) {writer.write("\t"+scores[i]);}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	private int[] getScore(SingleInterval region, Map<SingleInterval, Integer> counts) {
		int[] scores=new int[region.size()];
		for(int i=0; i<scores.length; i++) {
			int start=i+region.getReferenceStartPosition();
			SingleInterval newInterval=new SingleInterval(region.getReferenceName(), start, start+1);
			newInterval.setOrientation(Strand.antisense(region.getOrientation()));
			if(counts.containsKey(newInterval)) {
				scores[i]=counts.get(newInterval);
			}
		}
		return scores;
	}

	private Map<SingleInterval, Integer> score(File bam) {
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		
		int total=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment f=new SAMFragment(record);
			if(!record.getReadUnmappedFlag()){
				boolean passes=record.getFirstOfPairFlag();
				int startPos=record.getAlignmentStart();
				if(f.getOrientation().equals(Strand.POSITIVE)) {
					passes=record.getSecondOfPairFlag();
					startPos=record.getAlignmentEnd();
				}
					
				if(passes) {
					SingleInterval start=new SingleInterval(record.getReferenceName(), startPos, startPos+1);
					start.setOrientation(getStrand(record));
					if(!counts.containsKey(start)) {
						counts.put(start, 0);
					}
					int score=counts.get(start);
					score=score+1;
					counts.put(start, score);
					total++;
				}
			}
	
			if(total%1000000 ==0){System.err.println(total);}
		}
		
		
		
		reader.close();
		reads.close();
		return counts;
	}
	
	private static Strand getStrand(SAMRecord read) {
		Strand s=Strand.UNKNOWN;
		if(!read.getFirstOfPairFlag()) {
			if(read.getReadNegativeStrandFlag()) {s=Strand.NEGATIVE;}
			else{s=Strand.POSITIVE;}
		}
		else {
			if(read.getReadNegativeStrandFlag()) {s=Strand.POSITIVE;}
			else{s=Strand.NEGATIVE;}
		}
		return s;
	}
	
	public static Collection<SingleInterval> makeBedgraph(File input, int size) throws NumberFormatException, IOException {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
			String[] tokens=nextLine.split(" ");
			SingleInterval interval=new SingleInterval(tokens[0]);
			//int center=interval.getReferenceStartPosition()+Integer.parseInt(tokens[1]);
			//int start=center-(size/2);
			//int end=center+(size/2);
			//int end=Integer.parseInt(tokens[2]);
			
			int start=Integer.parseInt(tokens[1])-1;
			int end=Integer.parseInt(tokens[2]);
			
			SingleInterval newInterval=interval.trim(start, end);
			newInterval=new SingleInterval(newInterval.getReferenceName(), newInterval.getReferenceStartPosition()-100, newInterval.getReferenceEndPosition()+100);
			newInterval.setOrientation(interval.getOrientation());
			System.out.println(newInterval.toBED());
			rtrn.add(newInterval);
			
			/*if(interval.getOrientation().equals(Strand.POSITIVE)) {
				SingleInterval newInterval=new SingleInterval(interval.getReferenceName(), interval.getReferenceStartPosition()+start-1, interval.getReferenceStartPosition()+end);
				
			}*/
			
			}
			
		}
		reader.close();
		return rtrn;
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		File input=new File(args[0]);
		int size=100;
		Collection<SingleInterval> regions=makeBedgraph(input, size);
		
		File bam=new File(args[1]);
		String save=args[2];
		new FIMOOutput(bam, regions, save);
	}
	
}

package guttmanlab.core.spidr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class AllSecondReads {

	public AllSecondReads(File[] bams, String save) throws IOException {
		Map<SingleInterval, int[]> counts=new TreeMap<SingleInterval, int[]>();
		int[] totals=new int[bams.length];
		
		for(int i=0; i<bams.length; i++) {
			System.err.println(bams[i].getName()+" "+i+" "+bams.length);
			score(bams[i], i, counts, totals); 
		}
		
		write(save, counts, bams, totals);
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
	
	private void score(File bam, int i, Map<SingleInterval, int[]> counts, int[] totals) {
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
						counts.put(start, new int[totals.length]);
					}
					int[] scores=counts.get(start);
					scores[i]=scores[i]+1;
					counts.put(start, scores);
					total++;
				}
			}
	
			if(total%1000000 ==0){System.err.println(total);}
		}
		
		totals[i]=total;
		
		
		reader.close();
		reads.close();
		
	}

	private void write(String save, Map<SingleInterval, int[]> counts, File[] bams, int[] totals) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Region");
		for(int i=0; i<bams.length; i++) {
			writer.write("\t"+bams[i].getName());
		}
		writer.write("\n");
		
		writer.write("total");
		for(int i=0; i<totals.length; i++) {
			writer.write("\t"+totals[i]);
		}
		writer.write("\n");
		
		for(SingleInterval region: counts.keySet()) {
			int[] vals=counts.get(region);
			writer.write(region.toUCSCStrand());
			for(int i=0; i<vals.length; i++) {writer.write("\t"+vals[i]);}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	
	private static File[] getFiles(File[] files) {
		List<File> rtrn=new ArrayList<File>();
		for(int i=0; i<files.length; i++) {
			if(files[i].getName().endsWith("bam")) {rtrn.add(files[i]);}
		}
		
		File[] list=new File[rtrn.size()];
		for(int i=0; i<rtrn.size(); i++) {
			list[i]=rtrn.get(i);
		}
		
		return list;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
		File[] bams=getFiles(new File(args[0]).listFiles());
		String save=args[1];
		new AllSecondReads(bams, save);
		}
		else {System.err.println(usage);}
		
	}
	
	static String usage="args[0]=bams \n args[1]=save";
	
}

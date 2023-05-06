package guttmanlab.core.util;

import java.io.File;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PullReadsOverCoordinates {

	public static void main(String[] args) {
		File bam=new File(args[0]);
		SingleInterval coordinate=new SingleInterval(args[1]);
		String save=args[2];
		pullOverlapping(bam, coordinate, save);
		
	}

	private static void pullOverlapping(File bam, SingleInterval coordinate, String save) {
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(inputReader.getFileHeader(), false, new File(save+".bam"));
		
		int minLength=100;
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getNotPrimaryAlignmentFlag()) {
				if(overlaps(read, coordinate)) {
					if(readLength(read)>minLength) {
						writer.addAlignment(read);
					}
				}
			}
			
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount+" "+read.getReferenceName());}
		}
		reads.close();
		inputReader.close();
		writer.close();
	}

	private static int readLength(SAMRecord read) {
		return new SAMFragment(read).size();
	}

	private static boolean overlaps(SAMRecord read, SingleInterval coordinate) {
		if(read.getReferenceName().equals(coordinate.getReferenceName()) && read.getMateReferenceName().equals(coordinate.getReferenceName())) {
			if(read.getAlignmentStart()>=coordinate.getReferenceStartPosition()) {
				if(read.getAlignmentEnd()<=coordinate.getReferenceEndPosition()) {
					return true;
				}
			}
		}
		return false;
	}
	
}

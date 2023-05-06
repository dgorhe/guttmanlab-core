package guttmanlab.core.test;

import java.io.File;

import guttmanlab.core.annotation.SingleInterval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;

public class PullReadsOverRegion {

	public PullReadsOverRegion(File bamFile, SingleInterval region) {
		
		SAMFileReader inputReader= new SAMFileReader(bamFile);
		SAMRecordIterator iter=inputReader.queryOverlapping(region.getReferenceName(), region.getReferenceStartPosition(), region.getReferenceEndPosition());
		
		int counter=0;
		while(iter.hasNext()) {
			iter.next();
			counter++;
		}
		
		System.out.println(region.toUCSC()+"\t"+counter);
		inputReader.close();
	}
	
	
	public static void main(String[] args) {
		File bam=new File(args[0]);
		SingleInterval region=new SingleInterval(args[1]);
		new PullReadsOverRegion(bam, region);
	}
}

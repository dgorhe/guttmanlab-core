package guttmanlab.core.spidr;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class MakeCountTable {

	
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
	
	public static void main(String[] args) {
		File[] files=new File(args[0]).listFiles();
		for(int i=0; i<files.length; i++) {
			if(files[i].getName().endsWith(".bam")) {
				double count=getCount(files[i]);
				System.out.println(files[i].getName()+"\t"+count);
			}
		}
	}
	
}

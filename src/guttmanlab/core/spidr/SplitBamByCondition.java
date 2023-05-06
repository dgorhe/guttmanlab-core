package guttmanlab.core.spidr;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SplitBamByCondition {

	private static void filterBam(File input, String output) {
		SAMFileReader reader=new SAMFileReader(input);
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(output+".ALS.bam"));
		SAMFileWriter writer2=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(output+".ISO.bam"));
		
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			String name=record.getReadName();//[ALS]	[ISO]
			if(name.contains("[ALS")) {writer1.addAlignment(record);}
			if(name.contains("[ISO]")) {writer2.addAlignment(record);}
			
			
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		writer1.close();
		writer2.close();
	}
	
	
	public static void main(String[] args) {
		File[] files=new File(args[0]).listFiles();
		String saveDir=args[1];
		
		for(int i=0; i<files.length; i++) {
			String save=saveDir+"/"+files[i].getName();
			System.err.println(files[i].getName()+"\t"+i+"\t"+files.length+"\t"+save);
			filterBam(files[i], save);
		}
		
	}
	
}

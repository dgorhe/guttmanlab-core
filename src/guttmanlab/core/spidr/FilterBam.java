package guttmanlab.core.spidr;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class FilterBam {

	private static void filterBam(File input, String output) {
		SAMFileReader reader=new SAMFileReader(input);
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(output));
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			if(record.getReferenceName().startsWith("chr")) {writer.addAlignment(record);}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		writer.close();
	}
	
	
	public static void main(String[] args) {
		File[] files=new File(args[0]).listFiles();
		String saveDir=args[1];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName()+"\t"+i+"\t"+files.length);
			String save=saveDir+"/"+files[i].getName();
			filterBam(files[i], save);
		}
		
	}
	
}

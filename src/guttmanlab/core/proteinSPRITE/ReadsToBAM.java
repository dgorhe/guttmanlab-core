package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import guttmanlab.core.annotation.io.BEDFileIO;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class ReadsToBAM {

	public ReadsToBAM(File bamFile, String readFile, String save) throws IOException {
		Set<String> list=BEDFileIO.loadLineSet(readFile);
		System.err.println(list.size());
		
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
	
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(save));
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			if(list.contains(record.getReadName())){
				writer.addAlignment(record);
			}
			
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
		
		
		reader.close();
		reads.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		File bamFile=new File(args[0]);
		String readFile=args[1];
		String save=args[2];
		new ReadsToBAM(bamFile, readFile, save);
	}
	
	
}
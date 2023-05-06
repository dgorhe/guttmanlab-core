package guttmanlab.core.chip;

import java.io.File;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SubsampleBAM {

	
	public SubsampleBAM(File bam, String save, String chr) {
		
		SAMFileReader reader=new SAMFileReader(bam);
		
		
		SAMFileHeader fileHeader=reader.getFileHeader();
		
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(save));
		
		
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			if(record.getReferenceName().equalsIgnoreCase(chr)) {writer1.addAlignment(record);}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		
		reader.close();
		reads.close();
		writer1.close();
	}
	
	public static void main(String[] args) {
		File file=new File(args[0]);
		String save=args[1];
		String chr=args[2];
		new SubsampleBAM(file, save, chr);
	}
	
	
}

package guttmanlab.core.sars;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class FilterByReadQuality {

	public FilterByReadQuality(File bam, String save){
		SAMFileReader bamData= new SAMFileReader(bam);
		
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(bamData.getFileHeader(), false, new File(save));
		
		SAMRecordIterator reads=bamData.iterator();
		int counter=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1 && !read.getMateUnmappedFlag() && read.getReferenceName().equals(read.getMateReferenceName())){
				writer.addAlignment(read);
			}
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		writer.close();
		bamData.close();
	}
	
	public static void main(String[] args){
		new FilterByReadQuality(new File(args[0]), args[1]);
	}
	
}

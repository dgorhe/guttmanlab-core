package guttmanlab.core.trip;

import java.io.File;
import java.util.Collection;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SAMFragment;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class GrabSplicedReads {

	public static void main(String[] args) {
		
		SAMFileReader inputReader= new SAMFileReader(new File(args[0]));
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(inputReader.getFileHeader(), false, new File(args[1]));
		
		int totalCount=0;
		int spliced=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(SAMFragment.isSpliced(read)) {
				writer.addAlignment(read);
				spliced++;
			}
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount+" "+spliced);}
		}
				
		
		System.out.println(spliced+" "+totalCount);
		reads.close();
		inputReader.close();
		writer.close();
	}
	
}

package guttmanlab.core.smit;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.annotation.SAMFragment;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class NanoporeAnalysis {

	
	public NanoporeAnalysis(File bamFile, String save) throws IOException {
		//SAM to BED
		parse(bamFile, save);
		
	}

	private void parse(File sam, String save) throws IOException {
		SAMFileReader reader=new SAMFileReader(sam);
		SAMRecordIterator reads=reader.iterator();
		SAMFileHeader header=reader.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(header, false, new File(save));
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			record=SAMFragment.compressCIGAR(record);
			//SAMFragment frag=new SAMFragment(record);
			//writer.write(frag.toBED()+"\n");
			
			writer1.addAlignment(record);
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reads.close();
		reader.close();
		writer1.close();
		
	}
	
	public static void main(String[] args) throws IOException {
		File bam=new File(args[0]);
		String save=args[1];
		new NanoporeAnalysis(bam, save);
	}
	
}

package guttmanlab.core.spidr;

import java.io.File;
import java.util.Iterator;

import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;

public class ABCParsing {

	public static void main(String[] args) {
		//String barcode="ATCAC";
		String barcode=args[2];
		String rcBarcode=Sequence.reverseComplement(barcode);
		
		FastqReader reader=new FastqReader(new File(args[0]));
		Iterator<FastqRecord> iter=reader.iterator();
		
		FastqWriterFactory factory = new FastqWriterFactory();
		FastqWriter writer=factory.newWriter(new File(args[1]));
		
		
		int counter=0;
		int found=0;
		int rcFound=0;
		while(iter.hasNext()) {
			FastqRecord record=iter.next();
			String seq=record.getReadString();
			/*if(seq.contains(barcode)) {found++;}
			if(seq.contains(rcBarcode)) {
				writer.write(record);
				rcFound++;
			}*/
			
			if(seq.startsWith(rcBarcode)) {
				writer.write(record);
				rcFound++;
			}
			
			else if(seq.endsWith(barcode)) {
				//writer.write(record);
				found++;
			}
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter+" "+found+" "+rcFound);}
		}
		reader.close();
		writer.close();
	}
	
	
}

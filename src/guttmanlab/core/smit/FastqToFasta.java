package guttmanlab.core.smit;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class FastqToFasta {

	public static void main(String[] args) throws IOException {
		FastqReader p1=new FastqReader(new File(args[0]));
		FileWriter writer=new FileWriter(args[1]);
		
		while(p1.hasNext()){
			FastqRecord record=p1.next();
			writer.write(">"+record.getReadHeader()+"\n"+record.getReadString()+"\n");
		}
		
		p1.close();
		writer.close();
		
	}
	
}

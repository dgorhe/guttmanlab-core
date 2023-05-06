package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BAMToBarcodeString {

	
	public BAMToBarcodeString(File bam, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
	
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			writer.write(record.getReadName()+"\n");
			
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
		
		
		reader.close();
		reads.close();
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		File b=new File(args[0]);
		String save=args[1];
		new BAMToBarcodeString(b, save);
	}
	
}

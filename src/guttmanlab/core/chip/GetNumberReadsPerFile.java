package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;

public class GetNumberReadsPerFile {

	
	public static void main(String[] args) throws IOException {
		File[] files=new File(args[0]).listFiles();
		FileWriter writer=new FileWriter(args[1]);
		
		for(int i=0; i<files.length; i++) {
			if(files[i].getName().endsWith(".bam")) {
				double numberOfReads=countReads(files[i]);
				writer.write(files[i].getAbsolutePath()+"\t"+numberOfReads+"\t"+files[i].getName()+"\n");
			}
		}
		writer.close();
	}

	private static double countReads(File file) {
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		int counter=0;
		while(reads.hasNext()){
			reads.next();
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
			
		reader.close();
		reads.close();
			
		return counter;
	}
	
	
}

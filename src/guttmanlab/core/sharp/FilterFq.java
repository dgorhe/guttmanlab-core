package guttmanlab.core.sharp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class FilterFq {

	public FilterFq(File fastq1, File fastq2, File samFile, String save) throws IOException {
		
		Collection<String> readNames=new TreeSet<String>();
		SAMFileReader inputReader= new SAMFileReader(samFile);
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				readNames.add(read.getReadName());
			}
		}
				
		reads.close();
		inputReader.close();
		
		remove(fastq1, readNames, save+".1.fq");
		remove(fastq2, readNames, save+".2.fq");
		
	}

	private void remove(File fastq2, Collection<String> readNames, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		FastqReader p=new FastqReader(fastq2);
		
		Iterator<FastqRecord> iter=p.iterator();
		
		while(iter.hasNext()) {
			FastqRecord record=iter.next();
			String name=record.getReadHeader().split(" ")[0];
			
			if(!readNames.contains(name)) {
				//w.write(record);
				writer.write(record.getReadHeader()+"\n"+record.getReadString()+"\n+"+"\n"+record.getBaseQualityString()+"\n");
			}
		}
		p.close();
		writer.close();
		
	}
	
	
	public static void main(String[] args) throws IOException {
		File fastq1=new File(args[0]);
		File fastq2=new File(args[1]);
		File sam=new File(args[2]);
		String save=args[3];
		new FilterFq(fastq1, fastq2, sam, save);
	}
	
	
}

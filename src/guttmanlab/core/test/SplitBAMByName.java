package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SplitBAMByName {
	
	public SplitBAMByName(File bamFile) throws IOException{
		SAMFileReader reader=new SAMFileReader(bamFile);
		
		SAMFileHeader header=reader.getFileHeader();
		SAMRecordIterator reads=reader.iterator();
		
		Map<String, SAMFileWriter> fileWriters=new TreeMap<String, SAMFileWriter>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFileWriter writer;
			String name=getName(record);
			if(fileWriters.containsKey(name)){
				writer=fileWriters.get(name);
			}
			else{
				writer=makeWriter(header, bamFile, name);
				fileWriters.put(name, writer);
			}
			
			writer.addAlignment(record);
			
			counter++;
			if(counter%100000 ==0){System.err.println(counter +" "+fileWriters.size());}
		}
		
		close(fileWriters);
		
		reader.close();
		reads.close();
	}
	
	
	private String getName(SAMRecord record) {
		//HISEQ:1144:H3G5TBCX3:1:1101:1094:2234::[Rbm25_B]
		
		String name=record.getReadName();
		String truncated=name.split("\\[")[1].replaceAll("\\]", "");
		//System.err.println(name+" "+truncated);
		return truncated;
	}


	private SAMFileWriter makeWriter(SAMFileHeader header, File bamFile, String name) {
		File outFile=new File(bamFile.getAbsolutePath()+"."+name+".bam");
		
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(header, false, outFile);
		return writer1;
	}


	private void close(Map<String, SAMFileWriter> fileWriters) {
		for(String name: fileWriters.keySet()){
			SAMFileWriter writer=fileWriters.get(name);
			writer.close();
		}
		
	}




	
	public static void main(String[] args) throws IOException{
		if(args.length>0){
			File file=new File(args[0]);
			new SplitBAMByName(file);
	
		}
		else{System.err.println(usage);}
		
		
		
	}
	
	static String usage=" args[0]=bam file";


	
	
}

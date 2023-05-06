package guttmanlab.core.smit;

import java.io.File;
import java.io.FileWriter;
import java.util.Iterator;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;


public class ProcessNanoporeFastq {

	ProcessNanoporeFastq(File fastq, String primer1, String primer2, String save, int minLength){
		
		//go through fastq file and find primer 1 and primer 2... trim sequence in between
		FastqReader parser=new FastqReader(fastq);
		FastqWriter writer=new BasicFastqWriter(new File(save));
		
		int counter=0;
		Iterator<FastqRecord> iter=parser.iterator();
		while(iter.hasNext()){
			FastqRecord record=iter.next();
			FastqRecord newRecord=parse(record, primer1, primer2);
			write(writer, newRecord, minLength);
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		parser.close();
		writer.close();
	}

	private FastqRecord parse(FastqRecord record, String primer1, String primer2) {
		String read=record.getReadString();
		int startPos=read.indexOf(primer1);
		int endPos=read.indexOf(primer2);
		if(startPos>=0 && endPos>=0 && endPos>startPos){
			startPos+=primer1.length();
			//System.err.println(read+" "+startPos+" "+endPos);
			String subRead=read.substring(startPos, endPos);
			String subQuality=record.getBaseQualityString().substring(startPos, endPos);
			FastqRecord rtrn=new FastqRecord(record.getReadHeader(), subRead, record.getBaseQualityHeader(), subQuality);
			//FastqRecord(java.lang.String seqHeaderPrefix, java.lang.String seqLine, java.lang.String qualHeaderPrefix, java.lang.String qualLine) 
			//System.err.println(sub);
			return rtrn;
		}
		
		return null;
	}

	private void write(FastqWriter writer, FastqRecord newRecord, int minLength) {
		if(newRecord!=null && newRecord.getReadString().length()>minLength){
			writer.write(newRecord);
		}
	}
	
	public static void main(String[] args){
		if(args.length>3){
		File fastq=new File(args[0]);
		String primer1=args[1];
		String primer2=args[2];
		String save=args[3];
		int minLength=new Integer(args[4]);
		new ProcessNanoporeFastq(fastq, primer1, primer2, save, minLength);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fastq file \n args[1]=primer 1 \n args[2]=primer 2 \n args[3]=output \n args[4]=min length";
}

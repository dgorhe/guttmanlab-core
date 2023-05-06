package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class PullReads {

	public PullReads(File fq, File file2, String save) throws IOException{
		FastqReader reader=new FastqReader(fq);
		
		SamReader lines=SamReaderFactory.makeDefault().open((file2));
		
		//List<String> lines=BEDFileIO.loadLines(file2.getAbsolutePath());
		
		Collection<String> names=new TreeSet<String>();
		
		for(SAMRecord line: lines){
			String name=line.getReadName();
			names.add(name);
		}
		
		Collection<FastqRecord> reads=new ArrayList<FastqRecord>();
		Iterator<FastqRecord> iter=reader.iterator();
		while(iter.hasNext()){
			FastqRecord record=iter.next();
			String name=record.getReadHeader();
			//System.err.println(name);
			if(names.contains(name)){
				reads.add(record);
			}
		}
		
		write(save, reads);
		reader.close();
		lines.close();
	}

	private void write(String save, Collection<FastqRecord> reads) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(FastqRecord record: reads){
			writer.write(record.toString()+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		File fq=new File(args[0]);
		File file2=new File(args[1]);
		String save=args[2];
		new PullReads(fq, file2, save);
	}
	
}

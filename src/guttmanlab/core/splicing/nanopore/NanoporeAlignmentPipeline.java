package guttmanlab.core.splicing.nanopore;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;

public class NanoporeAlignmentPipeline {

	String rRNAs="/central/groups/guttman/mguttman/bowtieIndex/mouse/fasta/musmusculusRN45S.fa";
	String minimap="/groups/guttman/mguttman/minimap2/minimap2";
	String genome="/groups/guttman/genomes/mm10/mm10withchr.fa";
	String geneBed="/groups/guttman/mguttman/Refseq.mm10.bed";
	
	public NanoporeAlignmentPipeline(String fastqFile, String save) throws IOException, InterruptedException {
		//Step 1: Align to repetative RNAs
		align(rRNAs, fastqFile, save+".rRNA.sam");
		
		//Step 2: Filter fasta to remove reads that map to repeat RNAs
		String filteredFasta=filter(fastqFile, save+".rRNA.sam");
		
		//Step 3: Align to genome with genes
		align(genome, filteredFasta, save+".genome.sam", geneBed);
		
		//Step 4: Cleanup sam file
		cleanup(save+".genome.sam", save+".genome.compressedCigar.bam");
	}
	
	private void cleanup(String sam, String save) throws IOException {
		SAMFileReader reader=new SAMFileReader(new File(sam));
		SAMRecordIterator reads=reader.iterator();
		SAMFileHeader header=reader.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(header, false, new File(save));
			
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
				
			record=SAMFragment.compressCIGAR(record);
			writer1.addAlignment(record);
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
			
		reads.close();
		reader.close();
		writer1.close();	
	}

	private String filter(String fastqFile, String alignmentFile) throws IOException {
		Collection<String> mapped=new TreeSet<String>();
		
		SAMFileReader reader=new SAMFileReader(new File(alignmentFile));
		SAMRecordIterator reads=reader.iterator();
		
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			mapped.add(record.getReadName());
		}
		
		reader.close();
		reads.close();
		
		String output=fastqFile+".filteredByrRNA.fa";
		FastqReader p1=new FastqReader(new File(fastqFile));
		FileWriter writer=new FileWriter(output);
		while(p1.hasNext()){
			FastqRecord record=p1.next();
			String name=record.getReadHeader().split(" ")[0];
			if(!mapped.contains(name)){
				writer.write(">"+name+"\n"+record.getReadString()+"\n");
			}
		}
		
		writer.close();
		p1.close();
		return output;
	}

	private void align(String referenceFasta, String queryFasta, String output, String geneBed) throws IOException, InterruptedException {
		String cmd=minimap+ " -o "+output+" --sam-hit-only -ax splice --junc-bed "+geneBed+ " "+referenceFasta+" "+queryFasta;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
	}
	
	
	private void align(String referenceFasta, String queryFasta, String output) throws IOException, InterruptedException {
		String cmd=minimap+ " -o "+output+" --sam-hit-only -a "+referenceFasta+" "+queryFasta;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		String fastq=args[0];
		String save=args[1];
		new NanoporeAlignmentPipeline(fastq, save);
	}
	
}

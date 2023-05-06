package guttmanlab.core.test;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

public class SplitFastqByBarcode {

	public SplitFastqByBarcode(File fastq1, File fastq2, Map<String, String> barcodeToName, String save) throws IOException, InterruptedException{
		FastqReader r1=new FastqReader(fastq1);
		FastqReader r2=new FastqReader(fastq2);
		
		Map<String, Pair<FastqWriter>> writers=initialize(save, barcodeToName);
		
		while(r1.hasNext()){
			FastqRecord read1=r1.next();
			FastqRecord read2=r2.next();
			
			String seq=read1.getReadString().substring(0, 7);
			if(barcodeToName.containsKey(seq)){
				String name=barcodeToName.get(seq);
				Pair<FastqWriter> temp=writers.get(name);
				temp.getValue1().write(read1);
				temp.getValue2().write(read2);
			}
		}
		
		r1.close();
		r2.close();
		close(writers);
		
		//Align
		Collection<String> samFiles=align(writers, save);
		
		//Sort
		sort(samFiles);
		
		//Index
		
	}

	private void sort(Collection<String> samFiles) throws InterruptedException, IOException {
		for(String sam: samFiles){
			String cmd="java -jar /groups/guttman/software/picard.2.18.7/picard.jar SortSam I="+sam+" O="+sam+".sorted.bam SO=coordinate";
			System.err.println(cmd);
			Process p=Runtime.getRuntime().exec(cmd);
			p.waitFor();
			
			cmd="java -jar /groups/guttman/software/picard.2.18.7/picard.jar BuildBamIndex I="+sam+".sorted.bam";
			System.err.println(cmd);
			p=Runtime.getRuntime().exec(cmd);
			p.waitFor();
		}
		
	}

	private Collection<String> align(Map<String, Pair<FastqWriter>> writers, String save) throws IOException, InterruptedException {
		Collection<String> rtrn=new TreeSet<String>();
		for(String name: writers.keySet()){
			String cmd="bowtie2 --local -x /groups/guttman/genomes/mus_musculus/GRCm38/bowtie2/GRCm38.primary_assembly.genome -1 "+save+"."+name+".R1.fastq -2 "+save+"."+name+".R2.fastq -S "+name+".sam";
			System.err.println(name+" "+cmd);
			Process p=Runtime.getRuntime().exec(cmd);
			p.waitFor();
			rtrn.add(name+".sam");
		}
		return rtrn;
	}

	private Map<String, Pair<FastqWriter>> initialize(String save, Map<String, String> barcodeToName) {
		Map<String, Pair<FastqWriter>> rtrn=new TreeMap<String, Pair<FastqWriter>>();
		
		for(String barcode: barcodeToName.keySet()){
			String name=barcodeToName.get(barcode);
			File r1=new File(save+"."+name+".R1.fastq");
			File r2=new File(save+"."+name+".R2.fastq");
			FastqWriter f1=new BasicFastqWriter(r1);
			FastqWriter f2=new BasicFastqWriter(r2);
			Pair<FastqWriter> pair=new Pair<FastqWriter>(f1, f2);
			rtrn.put(name, pair);
		}
		
		return rtrn;
	}

	private void close(Map<String, Pair<FastqWriter>> writers) {
		for(String name: writers.keySet()){
			writers.get(name).getValue1().close();
			writers.get(name).getValue2().close();
		}
		
	}
	
	private static Map<String, String> parse(String save) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		List<String> lines=BEDFileIO.loadLines(save);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			if(tokens.length>2){
				String name=tokens[1];
				String seq=tokens[2];
				rtrn.put(seq, name);
			}
		}
		
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		File f1=new File(args[0]);
		File f2=new File(args[1]);
		Map<String, String> barcodeToName=parse(args[2]);
		String save=args[3];
		new SplitFastqByBarcode(f1, f2, barcodeToName, save);
	}

	
}

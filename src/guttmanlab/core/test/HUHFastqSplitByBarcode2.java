package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
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

public class HUHFastqSplitByBarcode2 {

	public HUHFastqSplitByBarcode2(File file1, File file2, Map<String, String> barcodeName, String save) throws IOException{
		
		FastqReader p1=new FastqReader(file1);
		FastqReader p2=new FastqReader(file2);
		
		Map<String, Pair<FileWriter>> writers=new TreeMap<String, Pair<FileWriter>>();
		
		//Collection<String> writers=new TreeSet<String>();
		
		while(p1.hasNext()){
			FastqRecord read1=p1.next();
			FastqRecord read2=p2.next();
			
			String barcode=getBarcode(read1);
			String name="other";
			
			if(barcodeName.containsKey(barcode)){name=barcodeName.get(barcode);}
			
			String fileName=save+"."+name;
			
			if(!writers.containsKey(fileName)){
				Pair<FileWriter> pair=new Pair<FileWriter>(new FileWriter(fileName+".R1.fq"), new FileWriter(fileName+".R2.fq"));
				writers.put(fileName, pair);
			}
			
			FileWriter writer1=writers.get(fileName).getValue1();
			FileWriter writer2=writers.get(fileName).getValue2();
				
			//Pair<FastqWriter> writer=map.get(name);
			writer1.write(read1.toString()+"\n");
			writer2.write(read2.toString()+"\n");
		
		}
		
		p1.close();
		p2.close();
		close(writers);
		
	}

	private void close(Map<String, Pair<FileWriter>> writers) throws IOException {
		for(String name: writers.keySet()){
			writers.get(name).getValue1().close();
			writers.get(name).getValue2().close();
		}
		
	}

	private String getBarcode(FastqRecord read1) {
		String seq=read1.getReadString();
		String barcode=seq.substring(8, 25);
		
		//System.err.println(seq+" "+barcode);
		
		return barcode;
	}

	private Pair<FastqWriter> makeWriter(String barcode, String save) {
		FastqWriter writer1=new BasicFastqWriter(new File(save+"."+barcode+"_R1.fq"));
		FastqWriter writer2=new BasicFastqWriter(new File(save+"."+barcode+"_R2.fq"));
		
		return new Pair<FastqWriter>(writer1, writer2);
	}

	/*private void close(Map<String, Pair<FastqWriter>> map) {
		for(String barcode: map.keySet()){
			Pair<FastqWriter> pair=map.get(barcode);
			pair.getValue1().close();
			pair.getValue2().close();
		}
		
	}*/
	
	private static Map<String, String> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(String line: lines){
			rtrn.put(line.split("\t")[0], line.split("\t")[1]);
		}
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File file1=new File(args[0]);
			File file2=new File(args[1]);
			Map<String, String> barcodeName=parse(args[2]);
			String save=args[3];
			new HUHFastqSplitByBarcode2(file1, file2, barcodeName, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	

	static String usage=" args[0]=fastq read1 \n args[1]=fastq read2 \n args[2]=barcode to name \n args[3]=save";
}

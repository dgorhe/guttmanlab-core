package guttmanlab.core.test;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;

public class FastqSplitByBarcode {

	public FastqSplitByBarcode(File file1, File file2, Map<String, String> barcodeToName, String save){
		
		FastqReader p1=new FastqReader(file1);
		FastqReader p2=new FastqReader(file2);
		
		Map<String, Pair<FastqWriter>> map=new TreeMap<String, Pair<FastqWriter>>();
		
		while(p1.hasNext()){
			FastqRecord read1=p1.next();
			FastqRecord read2=p2.next();
			
			String name=getBarcode(read2, barcodeToName);
			/*String name="unmatched";
			if(barcodeToName.containsKey(barcode)){
				name=barcodeToName.get(barcode);
			}*/
			
			if(!map.containsKey(name)){
				Pair<FastqWriter> writer=makeWriter(name, save);
				map.put(name, writer);
			}
				
			Pair<FastqWriter> writer=map.get(name);
			writer.getValue1().write(read1);
			writer.getValue2().write(read2);
		
		}
		
		p1.close();
		p2.close();
		close(map);
		
	}

	private String getBarcode(FastqRecord read2, Map<String, String> barcodeToName) {
		String seq=read2.getReadString();
		for(int i=5; i<=8; i++){
			String barcode=seq.substring(i, i+6);
			if(barcodeToName.containsKey(barcode)){return barcodeToName.get(barcode);}
		}
		return "unmatched";
		
		/*String barcode=seq.substring(7, 13);
		
		
		//System.err.println(seq+" "+barcode);
		return barcode;*/
	}

	private Pair<FastqWriter> makeWriter(String barcode, String save) {
		FastqWriter writer1=new BasicFastqWriter(new File(save+"."+barcode+"_R1.fq"));
		FastqWriter writer2=new BasicFastqWriter(new File(save+"."+barcode+"_R2.fq"));
		
		return new Pair<FastqWriter>(writer1, writer2);
	}

	private void close(Map<String, Pair<FastqWriter>> map) {
		for(String barcode: map.keySet()){
			Pair<FastqWriter> pair=map.get(barcode);
			pair.getValue1().close();
			pair.getValue2().close();
		}
		
	}
	
	private static Map<String, String> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(String line: lines){
			rtrn.put(line.split("\t")[0], line.split("\t")[1]);
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File file1=new File(args[0]);
			File file2=new File(args[1]);
			Map<String, String> barcodeToName=parse(args[2]);
			String save=args[3];
			new FastqSplitByBarcode(file1, file2, barcodeToName, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	

	static String usage=" args[0]=fastq read1 \n args[1]=fastq read2 \n args[2]=barcode to name file \n args[3]=save";
}

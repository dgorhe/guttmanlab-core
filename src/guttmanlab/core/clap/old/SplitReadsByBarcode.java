package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class SplitReadsByBarcode {

	
	public SplitReadsByBarcode(File fq1, String save) throws IOException{
		
		FastqReader reader1=new FastqReader(fq1);
		
		Map<String, Collection<String>> barcodeMap=new TreeMap<String, Collection<String>>();
		
		for(FastqRecord record1: reader1){
			String barcode=getBarcode(record1);
			String UMI=getUMI(record1);
			
			Collection<String> list=new TreeSet<String>();
			if(barcodeMap.containsKey(barcode)){list=barcodeMap.get(barcode);}
			list.add(UMI);
			barcodeMap.put(barcode, list);
				
		}
		
		
		write(save, barcodeMap);
		reader1.close();
	}
	
	
	private void write(String save, Map<String, Collection<String>> barcodeMap) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: barcodeMap.keySet()){
			writer.write(barcode+"\t"+barcodeMap.get(barcode).size()+"\n");
		}
		
		writer.close();
	}


	//With ORFs
	public SplitReadsByBarcode(File fq1, File fq2, String save) throws IOException{
		
		FastqReader reader1=new FastqReader(fq1);
		FastqReader reader2=new FastqReader(fq2);
		
		FileWriter writer=new FileWriter(save);
		
		Iterator<FastqRecord> iter=reader2.iterator();
		for(FastqRecord record1: reader1){
			FastqRecord record2=iter.next();
			String barcode=getBarcode(record1);
			String ORF=getORF(record2);
			String quality=getQuality(record2);
			writer.write("@"+record1.getReadHeader().split(" ")[0]+"_"+barcode+"\n"+ORF+"\n+\n"+quality+"\n");	
		}
		
		reader1.close();
		reader2.close();
		writer.close();
	}
	
	
	private String getQuality(FastqRecord record2) {
		//5 Ns + Constant region + ORF - Read 2
		//Full constant region (from sequencer):  GTTAGGGATAGGCTTACCGATATCAACCACTTTGTACAAGAAAGTTGGG
				
		String rtrn=record2.getBaseQualityString().substring(54, record2.getReadString().length());
		//System.err.println(record2.getReadString()+" "+rtrn);
		return rtrn;
	}


	private String getORF(FastqRecord record2) {
		//5 Ns + Constant region + ORF - Read 2
		//Full constant region (from sequencer):  GTTAGGGATAGGCTTACCGATATCAACCACTTTGTACAAGAAAGTTGGG
		
		String rtrn=record2.getReadString().substring(54, record2.getReadString().length());
		//System.err.println(record2.getReadString()+" "+rtrn);
		return rtrn;
	}

	
	private String getUMI(FastqRecord record1) {
		//CGC CTG CGA GAG GGA AAT CCA
		//5 Ns, Constant region, 7 Ns, Constant (BoxB) — Read 1 (edited)
		
		//5 Ns, CGCCTGCGAGAGGGAAATCCA, 7 Ns, BoxB (new).
		
		String rtrn=record1.getReadString().substring(0, 5);
		return rtrn;
	}

	private String getBarcode(FastqRecord record1) {
		//CGC CTG CGA GAG GGA AAT CCA
		//5 Ns, Constant region, 7 Ns, Constant (BoxB) — Read 1 (edited)
		
		//5 Ns, CGCCTGCGAGAGGGAAATCCA, 7 Ns, BoxB (new).
		
		String rtrn=record1.getReadString().substring(26, 26+7);
		return rtrn;
	}


	public static void main(String[] args)throws IOException{
		File fq1=new File(args[0]);
		File fq2=new File(args[1]);
		String save=args[2];
		new SplitReadsByBarcode(fq1, save);
	}
	
}

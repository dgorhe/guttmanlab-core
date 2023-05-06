package guttmanlab.core.test;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SplitProteinBarcodes {

	static String usage=" args[0]= RNA reads (BAM file) \n args[1]=protein reads (Fastq) \n args[2]=save";
	
	public static void main(String[] args){
		if(args.length<2){System.err.println(usage);}
		SAMFileReader rnaReads=new SAMFileReader(new File(args[0]));
		FastqReader reader=new FastqReader(new File(args[1]));
		String save=args[2];
		
		
		
		/*SAMFileReader reader=new SAMFileReader(new File(args[0]));
		SAMRecordIterator iter=reader.iterator();
		while(iter.hasNext()){
			SAMRecord record=iter.next();
			record.getReadName();
			
		}
		iter.close();
		reader.close();*/
		
		Map<String, Collection<String>> proteinToBarcodes=new TreeMap<String, Collection<String>>();
		
		while(reader.hasNext()){
			FastqRecord read=reader.next();
			String name=read.getReadHeader();
			//System.err.println(name);
			
			List<String> barcodes=splitBarcodes(name);
			merge(barcodes, proteinToBarcodes);
		}
		
		Map<String, String> barcodeToProtein=new TreeMap<String, String>();
		for(String protein: proteinToBarcodes.keySet()){
			System.err.println(protein+" "+proteinToBarcodes.get(protein).size());
			Collection<String> list=proteinToBarcodes.get(protein);
			for(String barcode: list){
				barcodeToProtein.put(barcode, protein);
				if(barcodeToProtein.containsKey(barcode)){
					String protein2=barcodeToProtein.get(barcode);
					if(!protein.equals(protein2)){
						System.err.println("ALREADY HAS "+barcode+" "+protein2+" "+protein);
					}
					}
			}
		}
		
		reader.close();
		
		SAMRecordIterator iter=rnaReads.iterator();
		SAMFileHeader header=rnaReads.getFileHeader();
		Map<String, SAMFileWriter> writers=new TreeMap<String, SAMFileWriter>();
		for(String protein: proteinToBarcodes.keySet()){
			String fileName=save+"."+protein+".bam";
			SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(header, false, new File(fileName));
			writers.put(protein, writer);
		}
		
		int count=0;
		while(iter.hasNext()){
			SAMRecord record=iter.next();
			String name=record.getReadName();
			String barcode=getBarcode(name);
			String protein=barcodeToProtein.get(barcode);
			//System.err.println(barcode +" ");
			if(protein!=null){
				SAMFileWriter writer=writers.get(protein);
				writer.addAlignment(record);
			}
			else{
				count++;
				//System.err.println(barcode);
			}
		}
		
		iter.close();
		
		for(String protein: writers.keySet()){
			writers.get(protein).close();
		}
		
		rnaReads.close();
		System.err.println("Skipped "+count);
	}

	private static String getBarcode(String name) {
		List<String> list=splitBarcodes(name);
		String barcodeString="";
		for(String barcode: list){
			if(barcode.startsWith("R1")){
				name=barcode;
			}
			else if(barcode.startsWith("R") && !barcode.startsWith("RPM")){
				barcodeString+="["+barcode+"]";
			}
		}
		
		return barcodeString;
		
	}

	private static void merge(List<String> barcodes, Map<String, Collection<String>> proteinToBarcodes) {
		String name="";
		String barcodeString="";
		for(String barcode: barcodes){
			if(barcode.startsWith("R1")){
				name=barcode;
			}
			else if(barcode.startsWith("R")){
				barcodeString+="["+barcode+"]";
			}
		}
		
		Collection<String> list=new TreeSet<String>();
		if(proteinToBarcodes.containsKey(name)){
			list=proteinToBarcodes.get(name);
		}
		list.add(barcodeString);
		proteinToBarcodes.put(name, list);
		
	}

	private static List<String> splitBarcodes(String name) {
		List<String> list=new ArrayList<String>();
		String barcode=name.split("::")[1];
		String[] tokens=barcode.split("\\[");
		for(int i=0; i<tokens.length; i++){
			String tag=tokens[i].replaceAll("\\]", "");
			list.add(tag);
			//System.err.println(tag);
		}
		return list;
	}
	
	
	
}

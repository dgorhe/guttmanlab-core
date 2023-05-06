package guttmanlab.core.trip;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


public class BarcodesToGenomicPosition {

	
	public BarcodesToGenomicPosition(File fastq1, File fastq2, File bam) throws IOException {
		Map<String, String> nameToBarcode=parseBarcodes(fastq2);
		//Map<String, String> nameToSequence=parseSequence(fastq1);
		
		Map<String, SingleInterval> nameToAlignment=parseAlignment(bam);
		
		Map<String, SingleInterval> barcodeToAlignment=new TreeMap<String, SingleInterval>();
		
		Collection<String> notunique=new TreeSet<String>();
		
		for(String name: nameToAlignment.keySet()) {
			SingleInterval alignment=nameToAlignment.get(name);
			String barcode=nameToBarcode.get(name);
			if(barcodeToAlignment.containsKey(barcode)) {
				SingleInterval temp=barcodeToAlignment.get(barcode);
				if(!temp.equals(alignment)) {notunique.add(barcode);}
			}
			barcodeToAlignment.put(barcode, alignment);
			//System.out.println(alignment.getReferenceName()+"\t"+alignment.getReferenceStartPosition()+"\t"+alignment.getReferenceEndPosition()+"\t"+barcode);
		}
		
		System.err.println(barcodeToAlignment.keySet().size()+" "+notunique.size());
		for(String barcode: barcodeToAlignment.keySet()) {
			if(!notunique.contains(barcode)) {
				SingleInterval alignment=barcodeToAlignment.get(barcode);
				System.out.println(alignment.getReferenceName()+"\t"+alignment.getReferenceStartPosition()+"\t"+alignment.getReferenceEndPosition()+"\t"+barcode);
			}
		}
		
	}

	private Map<String, SingleInterval> parseAlignment(File bam) {
		Map<String, SingleInterval> rtrn=new TreeMap<String, SingleInterval>();
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag()) {
				SingleInterval align=SAMFragment.getSingleInterval(read);
				rtrn.put(read.getReadName(), align);
			//System.err.println(read.getReadName());
			}
		}
		
		reads.close();
		inputReader.close();
		
		
		return rtrn;
	}

	private Map<String, String> parseSequence(File fastq1) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		FastqReader f=new FastqReader(fastq1);
		
		Iterator<FastqRecord> iter=f.iterator();
		
		
		while(iter.hasNext()) {
			FastqRecord seq=iter.next();
			String name=seq.getReadHeader();
			String seqBases=seq.getReadString();
			String subsequences=seqBases.substring(24, seqBases.length());
			//System.out.println(name+" "+subsequences);
			rtrn.put(name,  subsequences);
		}
		
		
		f.close();
		
		return rtrn;
	}

	private Map<String, String> parseBarcodes(File fastq2) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		FastqReader f=new FastqReader(fastq2);
		
		Iterator<FastqRecord> iter=f.iterator();
		
		
		while(iter.hasNext()) {
			FastqRecord seq=iter.next();
			String name=seq.getReadHeader().split(" ")[0];
			String seqBases=seq.getReadString();
			String barcode=seqBases.substring(138, 154);
			//System.err.println(name+" "+barcode);
			rtrn.put(name,  barcode);
		}
		
		
		f.close();
		
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException {
		File file1=new File(args[0]);
		File file2=new File(args[1]);
		File bam=new File(args[2]);
		new BarcodesToGenomicPosition(file1, file2, bam);
	}
	
}

package guttmanlab.core.trip;

import java.io.File;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.datastructures.Pair;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class ProcessCDNAReads {

	public ProcessCDNAReads(File cDNAFile, File splicedBAM, File unsplicedBAM) {
		//align to spliced reporter
		//align unaligned to unspliced reporter
		
		//de-duplicate and count by barcode
		
		Map<String, Collection<String>> cDNABarcode=pullBarcodes(cDNAFile); //barcode to names
		Map<String, Pair<Integer>> spliceState=getSpliceState(splicedBAM, unsplicedBAM); //name to splice state
		Map<String, String> readToUMI=getReadToUMI(splicedBAM, unsplicedBAM);
		Map<String, Pair<Integer>> barcodeToSpliceCounts=merge(cDNABarcode, spliceState, readToUMI);
		
		print(barcodeToSpliceCounts);
	}
	
	private static void print(Map<String, Pair<Integer>> barcodeToSpliceCounts) {
		for(String barcode: barcodeToSpliceCounts.keySet()) {
			Pair<Integer> pair=barcodeToSpliceCounts.get(barcode);
			System.out.println(barcode+"\t"+pair.getValue1()+"\t"+pair.getValue2());
		}
		
	}
	
	private Map<String, String> getReadToUMI(File splicedBAM, File unsplicedBAM) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		rtrn.putAll(getReadToUMI(splicedBAM));
		rtrn.putAll(getReadToUMI(unsplicedBAM));
		
		return rtrn;
	}

	private Map<String, String> getReadToUMI(File splicedBam) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		
		SAMFileReader inputReader= new SAMFileReader(splicedBam);
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			String name=read.getReadName();
			String seq=read.getReadString();
			String umi=getUMI(seq);
			rtrn.put(name, umi);
		}
				
		reads.close();
		inputReader.close();
		
		
		return rtrn;
	}

	private String getUMI(String seq) {
		String commonSeq="GGTGGTGAGGCCCTGGGCAG";
		int commonPos=seq.indexOf(commonSeq);
		if(commonPos!=-1) {
			String barcode=seq.substring(0, commonPos);
			//System.err.println(barcode+" "+seq+" "+barcode.length());
			return barcode;
		}
		return "empty";
	}

	private static Map<String, Pair<Integer>> merge(Map<String, Collection<String>> cDNABarcode, Map<String, Pair<Integer>> spliceState, Map<String, String> readToUMI) {
		Map<String, Pair<Integer>> rtrn=new TreeMap<String, Pair<Integer>>();
		for(String barcode: cDNABarcode.keySet()) {
			Collection<String> readNames=cDNABarcode.get(barcode);
			readNames=getUniqueNames(readNames, readToUMI);
			Pair<Integer> sum=sum(spliceState, readNames);
			//System.out.println(barcode+"\t"+sum.getValue1()+"\t"+sum.getValue2());
			rtrn.put(barcode, sum);
		}
		return rtrn;
	}
	
	private static Collection<String> getUniqueNames(Collection<String> readNames, Map<String, String> readToUMI) {
		Collection<String> rtrn=new TreeSet<String>();
		
		Collection<String> umis=new TreeSet<String>();
		
		for(String name: readNames) {
			if(readToUMI.containsKey(name)) {
				String UMI=readToUMI.get(name);
				if(!umis.contains(UMI)) {
					umis.add(UMI);
					rtrn.add(name);
				}
			}
		}
		
		return rtrn;
	}

	private static Pair<Integer> sum(Map<String, Pair<Integer>> spliceState, Collection<String> readNames) {
		Pair<Integer> rtrn=new Pair<Integer>(0,0);
		for(String name: readNames) {
			//System.err.println(name);
			if(spliceState.containsKey(name)) {
				Pair<Integer> val=spliceState.get(name);
				rtrn.setValue1(rtrn.getValue1()+val.getValue1());
				rtrn.setValue2(rtrn.getValue2()+val.getValue2());
			}
		}
		return rtrn;
	}
	
	private static Map<String, Pair<Integer>> getSpliceState(File splicedBam, File unsplicedBam) {
		Map<String, Pair<Integer>> rtrn=new TreeMap<String, Pair<Integer>>();
		int counter=0;
		
		SAMFileReader inputReader= new SAMFileReader(splicedBam);
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(read.getAlignmentStart()<93 && read.getAlignmentEnd()>93) {
				String name=read.getReadName();
				Pair<Integer> pair=new Pair<Integer>(1,0);
				rtrn.put(name, pair);
			}
			else {System.err.println(SAMFragment.getSingleInterval(read).toUCSC());}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
				
		reads.close();
		inputReader.close();
		
		
		counter=0;
		inputReader= new SAMFileReader(unsplicedBam);
		reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			String name=read.getReadName();
			Pair<Integer> pair=new Pair<Integer>(0,1);
			rtrn.put(name, pair);
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
				
		reads.close();
		inputReader.close();
		
		
		
		
		return rtrn;
	}
	
	
	private static Map<String, Collection<String>> pullBarcodes(File fastq2) {
		Map<String, Collection<String>> barcodeCount=new TreeMap<String, Collection<String>>();
		FastqReader f=new FastqReader(fastq2);
			
		Iterator<FastqRecord> iter=f.iterator();
			
		int counter=0;
		while(iter.hasNext()) {
			FastqRecord seq=iter.next();
			String seqBases=seq.getReadString();
			String barcode=seqBases;
			if(!barcodeCount.containsKey(barcode)) {barcodeCount.put(barcode, new TreeSet<String>());}
			Collection<String> list=barcodeCount.get(barcode);
			String name=seq.getReadHeader().split(" ")[0];
			list.add(name);
			barcodeCount.put(barcode, list);
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		f.close();
		return barcodeCount;
	}
	
	public static void main(String[] args) {
		File r2=new File(args[0]);
		File spliced=new File(args[1]);
		File unspliced=new File(args[2]);
		new ProcessCDNAReads(r2, spliced, unspliced);
	}
	
}

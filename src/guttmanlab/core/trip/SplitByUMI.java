package guttmanlab.core.trip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.pipeline.util.FastqSequence;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class SplitByUMI {

	private static Map<String, Collection<FastqRecord>> parseByUMI(File fastq){
		FastqReader f=new FastqReader(fastq);
		Map<String, Collection<FastqRecord>> rtrn=new TreeMap<String, Collection<FastqRecord>>();
		Iterator<FastqRecord> iter=f.iterator();
		
		int counter=0;
		while(iter.hasNext()) {
			FastqRecord seq=iter.next();
			String seqBases=seq.getReadString();
			String UMI=getUMI(seqBases);
			
			if(!rtrn.containsKey(UMI)) {
				rtrn.put(UMI, new ArrayList<FastqRecord>());
			}
			
			Collection<FastqRecord> list=rtrn.get(UMI);
			list.add(seq);
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		f.close();
		return rtrn;
		
	}
	
	
	private static void write(Map<String, Collection<String>> map, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String umi: map.keySet()) {
			Collection<String> list=map.get(umi);
			System.out.println(umi+" "+list.size());
			//If read1 and read2 are equal, only write once
			/*for(String name: list) {
				String mergedRead=r.getValue1().getReadString()+"_"+r.getValue2().getReadString();
				if(!used.contains(mergedRead)) {
					used.add(mergedRead);
					writer.write(toFastq(r)+"\n");
				}
				
			}*/
		}
		
		
		writer.close();
	}
	
	
	private static String toFastq(Pair<FastqRecord> r) {
		String name=r.getValue1().getReadHeader().split(" ")[0]+"_"+r.getValue2().getReadString();
		String sequence=r.getValue1().getReadString();
		String description="+";
		String quality=r.getValue1().getBaseQualityString();
		
		FastqSequence seq=new FastqSequence(name, sequence, description, quality);
		
		return seq.toFastq();
	}


	private static Map<String, Collection<String>> makeUnique(File fastq, File fastq2) {
		FastqReader f=new FastqReader(fastq);
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		Iterator<FastqRecord> iter=f.iterator();
		
		FastqReader f2=new FastqReader(fastq2);
		Iterator<FastqRecord> iter2=f2.iterator();
		
		int counter=0;
		while(iter.hasNext()) {
			FastqRecord seq=iter.next();
			FastqRecord seq2=iter2.next();
			
			String seqBases=seq.getReadString()+"_"+seq2.getReadString();
			
			
			if(!rtrn.containsKey(seqBases)) {
				rtrn.put(seqBases, new ArrayList<String>());
			}
			
			
			Collection<String> list=rtrn.get(seqBases);
			list.add(seq.getReadHeader());
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		f.close();
		f2.close();
		return rtrn;
		
	}
	
	/*private static Map<String, Collection<Pair<FastqRecord>>> makeUnique(File fastq, File fastq2) {
		FastqReader f=new FastqReader(fastq);
		Map<String, Collection<Pair<FastqRecord>>> rtrn=new TreeMap<String, Collection<Pair<FastqRecord>>>();
		Iterator<FastqRecord> iter=f.iterator();
		
		FastqReader f2=new FastqReader(fastq2);
		Iterator<FastqRecord> iter2=f2.iterator();
		
		int counter=0;
		while(iter.hasNext()) {
			FastqRecord seq=iter.next();
			String seqBases=seq.getReadString();
			FastqRecord seq2=iter2.next();
			
			
			if(!rtrn.containsKey(seqBases)) {
				rtrn.put(seqBases, new ArrayList<Pair<FastqRecord>>());
			}
			
			Pair<FastqRecord> pair=new Pair<FastqRecord>(seq, seq2);
			
			Collection<Pair<FastqRecord>> list=rtrn.get(seqBases);
			list.add(pair);
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		f.close();
		f2.close();
		return rtrn;
		
	}*/
	
	private static void writeSize(Map<String, Collection<FastqRecord>> map, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String umi: map.keySet()) {
			Collection<FastqRecord> reads=map.get(umi);
			writer.write(umi+"\t"+reads.size()+"\n");
		}
		
		writer.close();
	}
	
	
	private static String getUMI(String seq) {
		String commonSeq="GGTGGTGAGGCCCTGGGCAG";
							//GTGGTGAGGCCCTGGGCAGGCTGCTGGTTGT
		
		int commonPos=seq.indexOf(commonSeq);
		if(commonPos!=-1) {
			String barcode=seq.substring(0, commonPos);
			//System.err.println(barcode+" "+seq+" "+barcode.length());
			return barcode;
		}
		return "empty";
	}


	public static void main(String[] args) throws IOException {
		if(args.length>1) {
			//File fastq1=new File(args[0]);
			//File fastq2=new File(args[1]);
			//String save=args[2];
			//Map<String, Collection<FastqRecord>> map=parseByUMI(fastq);
			//Map<String, Collection<String>> map=makeUnique(fastq1, fastq2);
			//write(map, save);
			
			Map<String, Collection<FastqRecord>> map=parseByUMI(new File(args[0]));
			writeSize(map, args[1]);
		}
		else {System.err.println(usage);}
	}

	

	static String usage=" args[0]=fastq1 \n args[1]=save";
}

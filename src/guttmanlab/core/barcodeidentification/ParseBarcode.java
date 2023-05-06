package guttmanlab.core.barcodeidentification;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import htsjdk.samtools.fastq.FastqRecord;

public class ParseBarcode {

	static int gap=7;
	
	public static String parseBarcodes(FastqRecord record2, Map<String, String> b8, Map<String, String> b7, Map<String, String> b6, Map<String, String> b5, Map<String, String> b4, Map<String, String> b3, Map<String, String> b2) {
		String seq=record2.getReadString();
		String rtrn="NF";
		
		
		Kmer R8=get(seq, 0, b8);
		
		if(R8.isFound()){
			Kmer R7=get(seq, R8.end+gap, b7);
				
			Kmer R6=get(seq, R7.end+gap, b6);
				
			Kmer R5=get(seq,R6.end+gap, b5);
				
			Kmer R4=get(seq, R5.end+gap, b4);
				
			Kmer R3=get(seq, R4.end+gap, b3);
				
			Kmer R2=get(seq, R3.end+gap, b2);
				
			
			rtrn=R8+"_"+R7+"_"+R6+"_"+R5+"_"+R4+"_"+R3+"_"+R2;
		}
		
		
		return rtrn;
	}

	

	private static Kmer get(String seq, int start, Map<String, String> b2) {
		Collection<Integer> sizes=getSizes(b2);
		
		Kmer rtrn=new Kmer();
		
		for(int size: sizes){
			String substring=seq.substring(start, start+size);
			if(b2.containsKey(substring)){
				rtrn.setIsFound(true);
				rtrn.name=b2.get(substring);
				rtrn.seq=substring;
				rtrn.start=start;
				rtrn.end=start+size;
			}
		}
		
		
		return rtrn;
	}
	
	private static Collection<Integer> getSizes(Map<String, String> b2) {
		Collection<Integer> sizes=new TreeSet<Integer>();
		for(String b: b2.keySet()){sizes.add(b.length());}
		return sizes;
	}

	private static class Kmer{
		int start;
		int end;
		String seq;
		String name;
		boolean isFound;
		
		public Kmer(){
			this.isFound=false;
		}
		

		public void setIsFound(boolean b) {
			this.isFound=b;
		}

		
		
		public boolean isFound(){return isFound;}
		
		public String toString(){
			return "["+name+"]";
		}
	}
	
}

package guttmanlab.core.barcodeidentification;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class FindRepeats {

	public FindRepeats(File file1) {
		
		FastqReader f1=new FastqReader(file1);
		
		Iterator<FastqRecord> iter1=f1.iterator();
		
		while(iter1.hasNext()){
			FastqRecord record1=iter1.next();
			ArrayList<String> barcodes=parse(record1.getReadHeader());
			if(repeats(barcodes)) {
				System.out.println(record1.getReadHeader());
				System.out.println(record1.getReadString());
			}
		}
		
		f1.close();
	}

	private boolean repeats(ArrayList<String> barcodes) {
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		
		for(String b: barcodes) {
			int count=0; 
			if(counts.containsKey(b)) {count=counts.get(b);}
			count++;
			counts.put(b, count);
			if(count>1) {return true;}
		}
		return false;
	}

	private ArrayList<String> parse(String readHeader) {
		ArrayList<String> rtrn=new ArrayList<String>();
		
		//::[ROUND8_C3][ROUND7_B4][ROUND6_B10][ROUND5_B4][ROUND4_B4][ROUND3_C9][ROUND2_B5][NOT_FOUND]
		
		String sub=readHeader.split("::")[1];
		String[] tokens=sub.split("\\[");
		
		//System.err.println(sub);
		for(int i=1; i<tokens.length; i++) {
			String val=tokens[i].substring(0, tokens[i].length()-1);
			if(!val.equals("NOT_FOUND")) {
				rtrn.add(val);
			}
			//System.err.println(i+" "+val);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) {
		File file=new File(args[0]);
		new FindRepeats(file);
	}
	
}

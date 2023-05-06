package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.datastructures.Pair;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class CountNanobodySequences {
	
	//String adaptor="AATGATACGGCGACCAC";
	String adaptor="TGGTTCCGTCA";

	
	public CountNanobodySequences(FastqReader reader1, FastqReader reader2, String save) throws IOException{
		Map<String, Integer> counts1=getCounts(reader1);
		Map<String, Integer> counts2=getCounts(reader2);
		
		
		FileWriter writer=new FileWriter(save);
		
		Collection<String> seqs=new TreeSet<String>();
		seqs.addAll(counts1.keySet());
		seqs.addAll(counts2.keySet());
		
		for(String seq: seqs){
			int count1=0;
			int count2=0;
			if(counts2.containsKey(seq)){count2=counts2.get(seq);}
			if(counts1.containsKey(seq)){count1=counts1.get(seq);}
			writer.write(seq+"\t"+count1+"\t"+count2+"\n");
		}
		
		writer.close();
		
	}
	
	private Map<String, Integer> getCounts(FastqReader reader1) {
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		
		while(reader1.hasNext()){
			FastqRecord read=reader1.next();
			String seq=read.getReadString();
			if(seq.contains(adaptor)){
				seq=seq.substring(21, 51);
				int count=0;
				if(counts.containsKey(seq)){
					count=counts.get(seq);
				}
				count++;
				
				counts.put(seq, count);
			}
		}
		reader1.close();
		
		return counts;
	}

	public static void main(String[] args) throws IOException{
		FastqReader reader=new FastqReader(new File(args[0]));
		FastqReader reader2=new FastqReader(new File(args[1]));
		String save=args[2];
		new CountNanobodySequences(reader, reader2, save);
	}
	
}

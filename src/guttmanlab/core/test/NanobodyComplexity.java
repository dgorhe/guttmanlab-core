package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class NanobodyComplexity {
	
	public NanobodyComplexity(){}
	
	static String start="GCCTGCGAGAGGGAAATCCA";
	static String end="CTTCAGGGCCC";
	
	
	//BoxB
	public static void main(String[] args)throws IOException{
		
		//Read 2: NNNNN GGC CCT TCC GGG CAG CTG NNN NNN N GGGCCC , read 2
		
		FastqReader reader1=new FastqReader(new File(args[0]));
		Iterator<FastqRecord> iter1=reader1.iterator();
		
		FastqReader reader2=new FastqReader(new File(args[1]));
		Iterator<FastqRecord> iter2=reader2.iterator();
		
		FileWriter writer=new FileWriter(args[2]);
		
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		
		Map<String, Integer> barcodeCounts=new TreeMap<String, Integer>();
		
		while(iter2.hasNext()){
			FastqRecord read2=iter2.next();
			
			FastqRecord read1=iter1.next();
			
			String read=read2.getReadString();
			
			String UMI2=read1.getReadString().substring(0,5);
			String UMI=read.substring(0, 5);
			String barcode=read.substring(23, 30);
			
			String key=UMI2+"\t"+UMI+"\t"+barcode;
			int counter=0;
			if(counts.containsKey(key)){
				counter=counts.get(key);
			}
			counter++;
			counts.put(key, counter);
			
		}
			
		double total=0;
		for(String key: counts.keySet()){
			String barcode=key.split("\t")[2];
			int counter=0;
			if(barcodeCounts.containsKey(barcode)){
				counter=barcodeCounts.get(barcode);
			}
			counter++;
			total++;
			barcodeCounts.put(barcode, counter);
			writer.write(key+"\t"+counts.get(key)+"\n");
		}
			
		for(String barcode: barcodeCounts.keySet()){
			if(barcodeCounts.get(barcode)>10){
				double ratio=(double)barcodeCounts.get(barcode)/total;
				System.out.println(barcode+"\t"+barcodeCounts.get(barcode)+"\t"+ratio);
			}
		}
		
		reader1.close();
		reader2.close();
		writer.close();
		
	}
	
	
	//Nanobodies
	/*public static void main(String[] args)throws IOException{
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		
		FastqReader reader1=new FastqReader(new File(args[0]));
		FastqReader reader2=new FastqReader(new File(args[1]));
		String save=args[2];
		
		Iterator<FastqRecord> iter1=reader1.iterator();
		Iterator<FastqRecord> iter2=reader2.iterator();
		
		while(iter1.hasNext()){
			FastqRecord read1=iter1.next();
			FastqRecord read2=iter2.next();
			
			int count=0;
			String read=read1.getReadString()+read2.getReadString();
			if(counts.containsKey(read)){
				count=counts.get(read);
			}
			count++;
			counts.put(read,  count);
			
		}
		
		reader1.close();
		reader2.close();
		
		Map<Integer, Integer> dist=new TreeMap<Integer, Integer>();
		FileWriter writer=new FileWriter(save);
		for(String read: counts.keySet()){
			int count=counts.get(read);
			int score=0;
			if(dist.containsKey(count)){score=dist.get(count);}
			score++;
			dist.put(count, score);
			if(count>1){
				writer.write(read+"\t"+counts.get(read)+"\n");
			}
		}
		writer.close();
		
		for(Integer val: dist.keySet()){
			System.err.println(val+"\t"+dist.get(val));
		}
		
	}*/
}

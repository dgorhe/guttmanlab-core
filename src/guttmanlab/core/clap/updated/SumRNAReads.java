package guttmanlab.core.clap.updated;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SumRNAReads {

	private static Map<SingleInterval, Double> binBam(File bam, int binSize) {
		Map<SingleInterval, Double> map=new TreeMap<SingleInterval, Double>();
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator reads=inputReader.iterator();
		int totalCount=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			
			SAMFragment frag=new SAMFragment(read);
			SingleInterval bin=frag.bin(binSize);
			
			double count=0;
			if(map.containsKey(bin)){count=map.get(bin);}
			count++;
			map.put(bin, count);
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
		reads.close();
		inputReader.close();
		
		map=normalize(map, totalCount);
		//write(map, totalCount);
		return map;
	}
	
	
	private static void binBedgraph(File bedgraph, int binSize) throws IOException {
		Map<SingleInterval, Double> map=new TreeMap<SingleInterval, Double>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(bedgraph)));
		String nextLine;
		double totalCount=0;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens = nextLine.split("\t");
			String chr=(tokens[0]);
			int start=Integer.parseInt(tokens[1]);
			int end=Integer.parseInt(tokens[2]);
			SingleInterval annotation=new SingleInterval(chr, start, end);
			double score=Double.parseDouble(tokens[3]);
			
			SingleInterval bin=annotation.bin(binSize);
			
			double count=0;
			if(map.containsKey(bin)){count=map.get(bin);}
			count+=score;
			map.put(bin, count);
			totalCount+=score;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
		reader.close();
		
		normalize(map, totalCount);
	}
	
	private static Map<SingleInterval, Double> normalize(Map<SingleInterval, Double> map, double totalCount) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		for(SingleInterval r: map.keySet()) {
			double score=map.get(r);
			double norm=1000000*(double)score/(double)totalCount;
			rtrn.put(r, norm);
		}
		return rtrn;
	}

	private static void write(Map<SingleInterval, Double> map1, Map<SingleInterval, Double> map2, double cutoff) {
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		regions.addAll(map1.keySet());
		regions.addAll(map2.keySet());
		
		for(SingleInterval region: regions) {
			double val1=get(map1, region);
			double val2=get(map2, region);
			if(val1>cutoff || val2>cutoff) {
				double diff=val1-val2;
				System.out.println(region.toBedgraph(diff));
				//System.out.println(region.toUCSC()+"\t"+val1+"\t"+val2);
			}
		}
		
	}
	
	private static void write(Map<SingleInterval, Double> map1, Map<SingleInterval, Double> map2) {
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		regions.addAll(map1.keySet());
		regions.addAll(map2.keySet());
		
		for(SingleInterval region: regions) {
			double val1=get(map1, region);
			double val2=get(map2, region);
			System.out.println(region.toBedgraph(val1)+"\t"+val2);
		}
		
	}
	
	private static void write(Map<SingleInterval, Double> map1) {
		
		for(SingleInterval region: map1.keySet()) {
			double val1=map1.get(region);
			System.out.println(region.toBedgraph(val1));
		}
		
	}
	
	private static void writeBed(Map<SingleInterval, Double> map2, double cutoff) {
		
		for(SingleInterval region: map2.keySet()) {
			double val2=get(map2, region);
			if(val2>cutoff) {
				System.out.println(region.toBedgraph(val2));
			}
		}
		
		
	}

	
	private static double get(Map<SingleInterval, Double> map1, SingleInterval region) {
		if(map1.containsKey(region)) {return map1.get(region);}
		return 0;
	}


	public static void main(String[] args) throws IOException {
		File bam=new File(args[0]);
		int binSize=Integer.parseInt(args[1]);
		Map<SingleInterval, Double> map1=binBam(bam, binSize);
		//Map<SingleInterval, Double> map2=binBam(bam2, binSize);
		
		//double cutoff=Double.parseDouble(args[3]);
		//writeBed(map2, cutoff);
		
		write(map1);
		
		//binBedgraph(bam, binSize);
	}


	

	
	
}

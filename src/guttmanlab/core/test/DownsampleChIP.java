package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;
import guttmanlab.core.annotation.SingleInterval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class DownsampleChIP {

	public DownsampleChIP(File bamFile, double fractionToRetain, int resolution, String save) throws IOException{
		SAMFileReader inputReader= new SAMFileReader(bamFile);
		
		Map<SingleInterval, Integer> subCounts=new TreeMap<SingleInterval, Integer>();
		Map<SingleInterval, Integer> allCounts=new TreeMap<SingleInterval, Integer>();
		
		int totalCount=0;
		int subCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			double random=Math.random();
			SAMRecord read=reads.next();
			SingleInterval region=new SingleInterval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
			SingleInterval binned=region.bin(resolution);
			if(random<fractionToRetain){
				add(binned, subCounts);
				subCount++;
			}
			add(binned, allCounts);
			totalCount++;
		}
		
		double ratio=(double)subCount/(double)totalCount;
		System.err.println(subCount+" "+totalCount+" "+ratio);
		
		reads.close();
		inputReader.close();
		
		write(save+".all.bedgraph", allCounts);
		write(save+".downsample.bedgraph", subCounts);
		
		
	}
	
	private void add(SingleInterval binned, Map<SingleInterval, Integer> subCounts) {
		int count=0;
		if(subCounts.containsKey(binned)){count=subCounts.get(binned);}
		count++;
		subCounts.put(binned, count);
	}

	private static void write(String save, Map<SingleInterval, Integer> data) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: data.keySet()){
			writer.write(region.toBedgraph(data.get(region))+"\n");
		}
		
		writer.close();
	}
	
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File bam=new File(args[0]);
			double fractionToRetain=new Double(args[1]);
			int resolution=new Integer(args[2]);
			String save=args[3];
			new DownsampleChIP(bam, fractionToRetain, resolution, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=fraction to retain \n args[2]=bin resolution \n args[3]=save";
}

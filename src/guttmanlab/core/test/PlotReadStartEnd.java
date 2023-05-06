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

public class PlotReadStartEnd {

	
	public PlotReadStartEnd(File bamFile, String save) throws IOException{
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Integer> end=new TreeMap<SingleInterval, Integer>();
		Map<SingleInterval, Integer> start=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			SingleInterval region=new SingleInterval(read.getReferenceName(), read.getAlignmentEnd(), read.getAlignmentEnd());
			SingleInterval region2=new SingleInterval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentStart());
			int count=0;
			int count2=0;
			if(end.containsKey(region)){
				count=end.get(region);
			}
			count++;
			if(start.containsKey(region2)){
				count2=start.get(region2);
			}
			count2++;
			start.put(region2, count2);
			end.put(region, count);
			counter++;
			if(counter%10000==0){System.err.println(counter);}
		}
		
		reads.close();
		reader.close();
		
		write(save+".end.bedgraph", end);
		write(save+".start.bedgraph", start);
		
	}

	private void write(String save, Map<SingleInterval, Integer> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			int count=map.get(region);
			writer.write(region.toBedgraph(count)+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new PlotReadStartEnd(file, save);
	}
	
}

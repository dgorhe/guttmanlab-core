package guttmanlab.core.sars;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PlotSecondReadPileups {

	public PlotSecondReadPileups(String file1, String save) throws IOException{
		TreeMap<SingleInterval, Double> sample=score(file1);
		
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval pos: sample.keySet()){
			double sampleScore=sample.get(pos);
			writer.write(pos.getReferenceName()+"\t"+pos.getReferenceStartPosition()+"\t"+(pos.getReferenceStartPosition()+1)+"\t"+sampleScore+"\n");
			//writer.write(pos.getReferenceName()+"\t"+pos.getReferenceStartPosition()+"\t"+sampleScore+"\t"+inputScore+"\t"+conservativeInputScore+"\t"+ratio+"\t"+conservativeRatio+"\n");
		}
		
		writer.close();
	}
	
	
	private static TreeMap<SingleInterval, Double> score(String bam){
		SAMFileReader inputReader= new SAMFileReader(new File(bam));
		TreeMap<SingleInterval, Integer> positionCount=new TreeMap<SingleInterval, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getSecondOfPairFlag()){
				SingleInterval start=new SingleInterval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentStart()+1);
				int count=0;
				if(positionCount.containsKey(start)){count=positionCount.get(start);}
				count++;
				positionCount.put(start, count);
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		TreeMap<SingleInterval, Double> norm=new TreeMap<SingleInterval, Double>();
		for(SingleInterval pos: positionCount.keySet()){
			double ratio=1000000*((double)positionCount.get(pos)/(double)totalCount);
			norm.put(pos, ratio);
		}
		
		
		return norm;
	}
	
	private static double get(TreeMap<SingleInterval, Double> input, SingleInterval pos, boolean conservative) {
		double rtrn=0;
		
		if(input.containsKey(pos)){rtrn=input.get(pos);}
		else if(conservative){
			rtrn=getScores(input, pos.getReferenceName());
		}
		
		return rtrn;
	}

	private static double getScores(Map<SingleInterval, Double> input, String referenceName) {
		List<Double> list=new ArrayList<Double>();
		for(SingleInterval region: input.keySet()){
			if(region.getReferenceName().equals(referenceName)){list.add(input.get(region));}
		}
		
		double median=Statistics.quantile(list, 0.5);
		
		return median;
	}

	private static double get(TreeMap<SingleInterval, Double> input, SingleInterval pos) {
		return get(input, pos, false);
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>1){
			String sample=args[0];
			String save=args[1];
			new PlotSecondReadPileups(sample, save);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=sample bam \n args[1]=bedgraph";
	
	
}

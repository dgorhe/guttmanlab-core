package guttmanlab.core.clap.old;

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

public class PloitByReference {

	
	private static TreeMap<SingleInterval, Double> score(String bam){
		SAMFileReader inputReader= new SAMFileReader(new File(bam));
		TreeMap<SingleInterval, Integer> positionCount=new TreeMap<SingleInterval, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag()){
			//int start=read.getAlignmentStart();
				SingleInterval start=new SingleInterval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentStart()+1);
				int count=0;
				if(positionCount.containsKey(start)){count=positionCount.get(start);}
				count++;
				positionCount.put(start, count);
			}
				totalCount++;
			
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
			/*SingleInterval previous=input.floorKey(pos);
			SingleInterval next=input.higherKey(pos);
			//System.err.println(pos.toUCSC()+" "+previous.toUCSC()+" "+next.toUCSC());
			rtrn=getAverage(input, pos, previous, next);*/
			rtrn=getScores(input, pos.getReferenceName());
		}
		
		return rtrn;
	}
	
	private static double getAverage(TreeMap<SingleInterval, Double> input, SingleInterval pos, SingleInterval previous, SingleInterval next) {
		if(previous==null || !previous.getReferenceName().equals(pos.getReferenceName())){return input.get(next);}
		if(next==null || !next.getReferenceName().equals(pos.getReferenceName())){return input.get(previous);}
		
		double val1=input.get(previous);
		double val2=input.get(next);
		return (val1+val2)/2.0;
		
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
		TreeMap<SingleInterval, Double> sample=score(args[0]);
		TreeMap<SingleInterval, Double> input=score(args[1]);
		
		FileWriter writer=new FileWriter(args[2]);
		
		for(SingleInterval pos: sample.keySet()){
			double sampleScore=sample.get(pos);
			double inputScore=get(input, pos);
			double conservativeInputScore=get(input, pos, true);
			double ratio=(sampleScore+1)/(inputScore+1);
			double conservativeRatio=sampleScore/conservativeInputScore;
			writer.write(pos.getReferenceName()+"\t"+pos.getReferenceStartPosition()+"\t"+sampleScore+"\t"+inputScore+"\t"+conservativeInputScore+"\t"+ratio+"\t"+conservativeRatio+"\n");
		}
		
		writer.close();
	}

	
	
	
}

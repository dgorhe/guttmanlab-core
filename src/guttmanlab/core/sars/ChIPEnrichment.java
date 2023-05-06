package guttmanlab.core.sars;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class ChIPEnrichment {

	public ChIPEnrichment(String file1, String file2, String save) throws IOException{
		TreeMap<SingleInterval, Double> sample=score(file1);
		TreeMap<SingleInterval, Double> input=score(file2);
		
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval pos: sample.keySet()){
			double sampleScore=get(sample, pos);
			double inputScore=get(input, pos);
			double conservativeInputScore=get(input, pos, true);
			
			double scoreToUse=Math.max(inputScore, conservativeInputScore);
			if(inputScore==0){System.err.println(inputScore+" "+conservativeInputScore);}
			
	
			double ratio=Math.log((sampleScore)/(scoreToUse))/Math.log(2);
			writer.write(pos.getReferenceName()+"\t"+pos.getReferenceStartPosition()+"\t"+(pos.getReferenceEndPosition())+"\t"+ratio+"\n");
			
			//writer.write(pos.getReferenceName()+"\t"+pos.getReferenceStartPosition()+"\t"+sampleScore+"\t"+inputScore+"\t"+conservativeInputScore+"\t"+ratio+"\t"+conservativeRatio+"\n");
		}
		
		writer.close();
	}
	
	
	private TreeMap<SingleInterval, Double> score(String file1) throws IOException {
		TreeMap<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();

		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file1)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			SingleInterval region=new SingleInterval(tokens[0], new Integer(tokens[1]), new Integer(tokens[2]));
			Double val=new Double(tokens[3]);
			rtrn.put(region, val);
		}
		reader.close();
		return rtrn;
		
	}


	/*private static TreeMap<SingleInterval, Double> score(String bam){
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
	}*/
	
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
		
		double median=Statistics.min(list);
		
		return median;
	}

	private static double get(TreeMap<SingleInterval, Double> input, SingleInterval pos) {
		return get(input, pos, false);
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			String sample=args[0];
			String input=args[1];
			String save=args[2];
			System.err.println("new");
			new ChIPEnrichment(sample, input, save);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=sample bam \n args[1]=input bam \n args[2]=bedgraph";
	
	
}

package guttmanlab.core.sars;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PlotFullReadCounts {

	private int totalSampleCounts;

	public PlotFullReadCounts(File sampleBam, String save) throws IOException{
		Map<SingleInterval, Double> sampleScores=score(sampleBam);
		
		write(save, sampleScores);
	}
	
	private void write(String save, Map<SingleInterval, Double> sampleScores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: sampleScores.keySet()){
			double sample=get(sampleScores, region);
			double norm=1000000* (sample/this.totalSampleCounts);
			writer.write(region.toBedgraph(norm)+"\n");
		}
		
		writer.close();
	}
	
	private double get(Map<SingleInterval, Double> sampleScores, SingleInterval region) {
		if(sampleScores.containsKey(region)){return sampleScores.get(region);}
		return 0;
	}
	
	private TreeMap<SingleInterval, Double> score(File bam){
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<SingleInterval, Double> positionCount=new TreeMap<SingleInterval, Double>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				int start=read.getAlignmentStart();
				int end=read.getAlignmentEnd();
				//int mateStart=read.getMateAlignmentStart();
				//System.err.println(start+" "+mateStart);
				Collection<SingleInterval> bin=bin(read.getReferenceName(), start, end);
				for(SingleInterval r: bin){
					double score=0;
					if(positionCount.containsKey(r)){score=positionCount.get(r);}
					score++;
					positionCount.put(r, score);
				}
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		this.totalSampleCounts=totalCount;
		
		return positionCount;
	}
	
	

	private Collection<SingleInterval> bin(String referenceName, int start, int end) {
		int min=Math.min(start, end);
		int max=Math.max(start, end);
		
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(int i=min; i<max; i++){
			SingleInterval region=new SingleInterval(referenceName, i, i+1);
			rtrn.add(region);
		}
		return rtrn;
	}

	private Collection<SingleInterval> bin(SAMRecord read, int binSize) {
		
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		int startIndex=read.getAlignmentStart();
		
		for(int i=0; i<binSize; i++){
			int newStart=startIndex+i;
			SingleInterval region=new SingleInterval(read.getReferenceName(), newStart, newStart+1);
			rtrn.add(region);
			newStart=startIndex-i;
			region=new SingleInterval(read.getReferenceName(), newStart, newStart+1);
			rtrn.add(region);
		}
		return rtrn;
		
		
		
		
		
	}

	public static void main(String[] args) throws IOException{
		if(args.length>1){
			File sample=new File(args[0]);
			String save=args[1];
			new PlotFullReadCounts(sample, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sample bam \n args[1]=save";
	
	
}

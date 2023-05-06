package guttmanlab.core.sars;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PlotSlidingWindow {

	private int totalSampleCounts;
	private int totalInputCounts;

	public PlotSlidingWindow(File sampleBam, File inputBam, String save) throws IOException{
		Map<SingleInterval, Double> sampleScores=score(sampleBam, true);
		Map<SingleInterval, Double> inputScores=score(inputBam, false);
		Map<String, Double> geneScore=average(inputScores);
		
		write(save, sampleScores, inputScores, geneScore);
	}
	
	private Map<String, Double> average(Map<SingleInterval, Double> inputScores) {
		Map<String, List<Double>> lists=new TreeMap<String, List<Double>>();
		
		for(SingleInterval region: inputScores.keySet()){
			String gene=region.getReferenceName();
			List<Double> list=new ArrayList<Double>();
			if(lists.containsKey(gene)){list=lists.get(gene);}
			list.add(inputScores.get(region));
			lists.put(gene, list);
		}
		
		
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(String gene: lists.keySet()){
			double median=Statistics.quantile(lists.get(gene), 0.5);
			rtrn.put(gene, median);
		}
		
		return rtrn;
	}

	public PlotSlidingWindow(File sampleBam, String save) throws IOException{
		Map<SingleInterval, Double> sampleScores=score(sampleBam, true);
		
		write(save, sampleScores);
	}
	
	
	private void write(String save, Map<SingleInterval, Double> sampleScores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: sampleScores.keySet()){
			double sample=get(sampleScores, region);
			//double input=get(inputScores, region);
			//double geneInput=get(geneInputCount, region.getReferenceName());
			//double maxInput=Math.max(input, 1);
			double ratio=1000000*(sample/this.totalSampleCounts);
			writer.write(region.toBedgraph(ratio)+"\n");
		}
		
		writer.close();
	}
	
	private void write(String save, Map<SingleInterval, Double> sampleScores, Map<SingleInterval, Double> inputScores, Map<String, Double> geneInputScore) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: sampleScores.keySet()){
			double sample=get(sampleScores, region);
			double input=Math.max(1, get(inputScores, region));
			double geneInput=get(geneInputScore, region.getReferenceName());
			double maxInput=Math.max(input, geneInput);
			double ratio=(sample/this.totalSampleCounts)/(maxInput/this.totalInputCounts);
			writer.write(region.toBedgraph(ratio)+"\n");
		}
		
		writer.close();
	}
	
	private double get(Map<String, Double> geneInputScore, String referenceName) {
		if(geneInputScore.containsKey(referenceName)){return geneInputScore.get(referenceName);}
		return 0;
	}

	private double get(Map<SingleInterval, Double> sampleScores, SingleInterval region) {
		if(sampleScores.containsKey(region)){return sampleScores.get(region);}
		return 0;
	}
	
	private TreeMap<SingleInterval, Double> score(File bam, boolean sample){
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<SingleInterval, Double> positionCount=new TreeMap<SingleInterval, Double>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag()){
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
		
		if(sample){this.totalSampleCounts=totalCount;}
		else{this.totalInputCounts=totalCount;}
		
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
		if(args.length>2){
			File sample=new File(args[0]);
			File control=new File(args[1]);
			String save=args[2];
			new PlotSlidingWindow(sample, control, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sample bam \n args[1]=control bam \n args[2]=save";
	
	
}

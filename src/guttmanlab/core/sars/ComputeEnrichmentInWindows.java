package guttmanlab.core.sars;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.math.Statistics;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class ComputeEnrichmentInWindows {

	private int totalSampleCounts;
	private int totalInputCounts;
	private double minCount=10;
	private double minEnrichment=2;
	
	public ComputeEnrichmentInWindows(File sample, File input, String save, int binSize) throws IOException{
		Map<SingleInterval, Double> sampleScores=score(sample, binSize, true);
		Map<SingleInterval, Double> inputScores=score(input, binSize, false);
		Map<String, Double> geneInputCount=medianCount(inputScores);
		
		//write(save+"."+sample.getName()+".bedgraph", sampleScores);
		//write(save+"."+input.getName()+".bedgraph", inputScores, geneInputCount);
		writeRatio(save+".ratio.bedgraph", sampleScores, inputScores, geneInputCount);
		
	}
	
	/**
	 * P values based on window normalization
	 * @return window normalized p-value
	 */
	/*public double getWindowPVal() {
		int k=this.getSampleReadCount();
		int n=getSampleReadCount()+getMaxInputCount();
		
		if(n>0){
			double elution=(double)this.getSampleTotalCount()/(double)this.getNumberSampleWindows();
			double input=(double) this.getInputTotalCount()/(double)this.getNumberInputWindows();
			double p=elution/(elution+input);
			Binomial b=new Binomial(n, p);
			return 1-b.cdf(k);
		}
		return 1.0;
		
	}*/
	
	public double getPValue(double sampleCounts, double inputCounts) {
		//Return binomial p
		
		int n=new Double(sampleCounts+inputCounts).intValue();
		if(n>0){
			double p=(double)this.totalSampleCounts/((double)this.totalSampleCounts+(double)this.totalInputCounts);
			Binomial b=new Binomial(n, p);
			return 1-b.cdf(sampleCounts);
		}
		return 1.0;
	}
	
	
	private void write(String save, Map<SingleInterval, Double> inputScores, Map<String, Double> geneInputCount) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: inputScores.keySet()){
			String name=region.getReferenceName();
			double regionScore=inputScores.get(region);
			double geneScore=geneInputCount.get(name);
			writer.write(region.toUCSC()+"\t"+regionScore+"\t"+geneScore+"\t"+this.totalInputCounts+"\n");
		}
		
		writer.close();
	}



	private Map<String, Double> medianCount(Map<SingleInterval, Double> inputScores) {
		Map<String, List<Double>> vals=new TreeMap<String, List<Double>>();
		
		for(SingleInterval region: inputScores.keySet()){
			String name=region.getReferenceName();
			List<Double> list=new ArrayList<Double>();
			if(vals.containsKey(name)){list=vals.get(name);}
			list.add(inputScores.get(region));
			vals.put(name, list);
		}
		
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(String name: vals.keySet()){
			List<Double> list=vals.get(name);
			double median=Statistics.quantile(list, 0.5);
			rtrn.put(name, median);
		}
		
		
		return rtrn;
	}



	private void write(String save, Map<SingleInterval, Double> inputScores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: inputScores.keySet()){
			double score=get(inputScores, region);
			writer.write(region.tobedgraph(score)+"\n");
		}
		
		writer.close();
	}



	private void writeRatio(String save, Map<SingleInterval, Double> sampleScores, Map<SingleInterval, Double> inputScores, Map<String, Double> geneInputCount) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: sampleScores.keySet()){
			double sample=get(sampleScores, region);
			double input=get(inputScores, region);
			double geneInput=get(geneInputCount, region.getReferenceName());
			double maxInput=Math.max(input, geneInput);
			maxInput=Math.max(1, maxInput);
			double ratio=(sample/this.totalSampleCounts)/(maxInput/this.totalInputCounts);
			double p=this.getPValue(sample, maxInput);
			//writer.write(region.toBedgraph(ratio)+"\n");
			if(p<0.001 && sample>minCount && ratio>minEnrichment){
				System.err.println(region.toUCSC()+"\t"+sample+"\t"+input+"\t"+ratio+"\t"+p);
				writer.write(region.toBedgraph(ratio)+"\n");
				//writer.write(region.toUCSC()+"\t"+sample+"\t"+input+"\t"+geneInput+"\t"+ratio+"\t"+p+"\n");
			}
		}
		
		writer.close();
	}



	private double get(Map<String, Double> geneInputCount, String referenceName) {
		if(geneInputCount.containsKey(referenceName)){return geneInputCount.get(referenceName);}
		return 0;
	}

	private double get(Map<SingleInterval, Double> sampleScores, SingleInterval region) {
		if(sampleScores.containsKey(region)){return sampleScores.get(region);}
		return 0;
	}



	private TreeMap<SingleInterval, Double> score(File bam, int binSize, boolean sample){
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<SingleInterval, Double> positionCount=new TreeMap<SingleInterval, Double>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				SingleInterval bin=bin(read, binSize);
				double score=0;
				if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
				score++;
				positionCount.put(bin, score);
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
	
	private static SingleInterval bin(SAMRecord read, int resolution) {
		int startIndex=read.getAlignmentStart()/resolution;
		int newStart=startIndex*resolution;
		int genomicLength=read.getAlignmentEnd()-read.getAlignmentStart();
		int newEnd=newStart+resolution;
		SingleInterval newInterval=new SingleInterval(read.getReferenceName(), newStart, newEnd);
		return newInterval;
	}



	public static void main(String[] args) throws IOException{
		if(args.length>3){
		File sample=new File(args[0]);
		File input=new File(args[1]);
		String save=args[2];
		int binSize=new Integer(args[3]);
		new ComputeEnrichmentInWindows(sample, input, save, binSize);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sample bam \n args[1]=input bam \n args[2]=save \n args[3]=bin size";
		
}

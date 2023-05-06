package guttmanlab.core.sharp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.barcodeidentification.PeakCalling;
import guttmanlab.core.datastructures.IntervalTree;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class CLAPAnalysis {
	
	Map<Annotation, Double> differentialPVal;
	Map<Annotation, Double> differentialEnrich;
	Map<Annotation, String> windowsToGene;
	double alpha=0.001;
	double minEnrichment=3.0;
	int minCount=20;
	Map<Annotation, Integer> maxInputScore;
	
	public CLAPAnalysis(File sampleBAM, File inputBAM, int binSize, Map<String, IntervalTree<Annotation>> genes, String save) throws IOException {
		PeakCalling samplePeaks= new PeakCalling(sampleBAM, binSize, genes); 
		PeakCalling inputPeaks= new PeakCalling(inputBAM, binSize, genes); 
		
		Map<Annotation, Double> geneScores=inputPeaks.sampleGeneScores();
		
		System.err.println("done sampling");
		
		//compute differential
		differential(samplePeaks, inputPeaks, geneScores);
		
		System.err.println("done differential");
		
		writeScores(samplePeaks, inputPeaks, genes, save);
		
		System.err.println("done writing");
	}
	
	

	private void writeScores(PeakCalling samplePeaks, PeakCalling inputPeaks, Map<String, IntervalTree<Annotation>> genes, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("window\tnumber of exons\tstrand\tfeature name\tsample count\tinput count\tenrichment (intra)\tenrichment (inter)\tp-val (intra)\tp-val (inter)\t window normalized p-value\twindow normalized enrichment\n");
		
		for(Annotation a: this.differentialEnrich.keySet()) {
			String name=this.windowsToGene.get(a);
			double enrichment2=this.differentialEnrich.get(a);
			double p2=this.differentialPVal.get(a);
			if(samplePeaks.getEnrichments().containsKey(a) && samplePeaks.getPvalues().containsKey(a)) {
				double enrichment1=samplePeaks.getEnrichments().get(a);
				double p1=samplePeaks.getPvalues().get(a);
				int count1=samplePeaks.getCounts().get(a);
				int count2=Math.max(1,maxInputScore.get(a));
				double windowEnrichment=getWindowNormEnrichment(count1, count2, samplePeaks.getTotalCount(), inputPeaks.getTotalCount(), samplePeaks.getNumberOfWindows(), inputPeaks.getNumberOfWindows());
				double windowP=getWindowNormPValue(count1, count2, samplePeaks.getTotalCount(), inputPeaks.getTotalCount(), samplePeaks.getNumberOfWindows(), inputPeaks.getNumberOfWindows());
				writer.write(a.toBlocks()+"\t"+a.getNumberOfBlocks()+"\t"+a.getOrientation()+"\t"+name+"\t"+count1+"\t"+count2+"\t"+enrichment1+"\t"+enrichment2+"\t"+p1+"\t"+p2+"\t"+windowP+"\t"+windowEnrichment+"\n");
			}		
		}
		writer.close();
	}
	
	
	
	
	
	private double getWindowNormPValue(int sampleReadCount, int inputReadCount, int sampleTotalCount, int inputTotalCount, int sampleWindows, int inputWindows) {
		int k=sampleReadCount;
		int n=sampleReadCount+inputReadCount;
		
		if(n>0){
			double elution=(double)sampleTotalCount/(double)sampleWindows;
			double input=(double) inputTotalCount/(double)inputWindows;
			double p=elution/(elution+input);
			Binomial b=new Binomial(n, p);
			return 1-b.cdf(k);
		}
		return 1.0;
	}



	public double getWindowNormEnrichment(int sampleCount, int maxInputCount, int sampleTotalCount, int inputTotalCount, int sampleWindows, int inputWindows){
		double numerator=(double)sampleCount/(double)maxInputCount;
		double elution=(double)sampleTotalCount/(double)sampleWindows;
		double input=(double)inputTotalCount/(double)inputWindows;
		double denominator=elution/input;
		return numerator/denominator;
	}


	private void differential(PeakCalling samplePeaks, PeakCalling inputPeaks, Map<Annotation, Double> geneScores) {
		this.differentialEnrich=new TreeMap<Annotation, Double>();
		this.differentialPVal=new TreeMap<Annotation, Double>();
		
		Map<Annotation, Integer> counts1=samplePeaks.getCounts();
		this.maxInputScore=getInputScores(inputPeaks, samplePeaks, geneScores);
		
		int total1=samplePeaks.getTotalCount();
		int total2=inputPeaks.getTotalCount();
		
		for(Annotation r: counts1.keySet()) {
			int window1=counts1.get(r);
			int window2=maxInputScore.get(r);
			if(window2==0) {window2=1;}
			
			double enrichment=getEnrichment(window1, window2, total1, total2);
			double p=getPValue(window1, window2, total1, total2);
			
			this.differentialEnrich.put(r, enrichment);
			this.differentialPVal.put(r, p);
		}
		
	}
	
	private Map<Annotation, Integer> getInputScores(PeakCalling inputPeaks, PeakCalling samplePeaks, Map<Annotation, Double> geneScore) {
		this.windowsToGene=new TreeMap<Annotation, String>();
		
		Map<Annotation, Integer> rtrn=new TreeMap<Annotation, Integer>();
		
		Collection<Annotation> bins=samplePeaks.getCounts().keySet();
		for(Annotation bin: bins) {
			int count1=0;
			if(inputPeaks.getCounts().containsKey(bin)) {
				count1=inputPeaks.getCounts().get(bin);
			}
			Annotation maxGene=inputPeaks.getMaxGene(bin);
			int count2=0;
			if(maxGene!=null) {
				this.windowsToGene.put(bin, maxGene.getName());
				if(geneScore.containsKey(maxGene)) {
					count2=(int)Math.ceil(geneScore.get(maxGene));
				}
			}
			else {
				this.windowsToGene.put(bin, "intergenic");
				double mean=inputPeaks.getScore(bin);
				count2=(int)Math.ceil(mean);
				/*if(mean>0) {
					double scaled=inputPeaks.getPerms(mean);
					count2=(int)Math.ceil(scaled);
				}*/
			}
			
			rtrn.put(bin, Math.max(count1, count2));
		}
		
		
		return rtrn;
	}


	private static double getEnrichment(double window1, double window2, int total1, int total2) {
		double num=((double)window1)/(double)total1;
		double denom=((double)window2)/(double)total2;
		return num/denom;
	}
	
	private static double getPValue(int window_sampleCount, int window_inputCount, int total_sampleCount, int total_inputCount) {
		//Return binomial p
		int n=window_inputCount+window_sampleCount;
		if(n>0){
			double p=(double)total_sampleCount/(double)(total_inputCount+total_sampleCount);
			Binomial b=new Binomial(n, p);
			return 1-b.cdf(window_sampleCount);
		}
		return 1.0;
	}
	
	private static int readCount(File bam) {
		SAMFileReader inputReader= new SAMFileReader(bam);
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
		reads.close();
		inputReader.close();
		return totalCount;
	}


	public static void main(String[] args) throws IOException {
		if(args.length>4) {
			File sampleBam=new File(args[0]);
			File inputBam=new File(args[1]);
			int binSize=Integer.parseInt(args[2]);
			
			String save=args[4];
			
			Map<String, IntervalTree<Annotation>> genes=BEDFileIO.loadTreePlusIntronsAnnotation(args[3]);
			
			System.err.println("V6");
			new CLAPAnalysis(sampleBam, inputBam, binSize, genes, save);
		}
		else {System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=sample bam \n args[1]=input bam \n args[2]=bin size \n args[3]=genes (bed) \n args[4]=save";
}

package guttmanlab.core.clap.old;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.math.ScanStat;
import jsc.distributions.Binomial;

public class CLAPScore {
	
	private int sampleCount;
		private int inputCount;
		private int sampleTotal;
		private int inputTotal;
		Annotation parentAnnotation;
		private int inputPercentile;
		boolean hasPercentile;
		private int numberInputWindows;
		private int numberSampleWindows;
		private int geneSampleTotal;
		private int geneInputTotal;
		private int geneLength;
		private Annotation window;
		private String geneName;
		private int totalGeneSizes;
		
		
		public static String getHeader(){
			return "Gene Name \t window \t Normalized Sample Score \t Normalized Input Score \t Global Enrichment \t global p-value \t global scan p-value \t sampleCount \t input Count \t  Max Input Count \t Window Normalized Enrichment \t Window P-Value \t Local P-Value \t Local Scan P-Value";	
			
		}

		public CLAPScore(Annotation window,int sampleCount, int inputCount, int sampleTotal, int inputTotal, int totalGeneSizes){
			this.sampleCount=sampleCount;
			this.sampleTotal=sampleTotal;
			this.inputCount=inputCount;
			this.inputTotal=inputTotal;
			hasPercentile=false;
			this.window=window;
			this.geneName=window.toUCSC();
			this.totalGeneSizes=totalGeneSizes;
		}

		public void setGeneInputTotal(int geneInputTotal) {
			this.geneInputTotal=geneInputTotal;
		}

		public void setGeneSampleTotal(int geneSampleTotal) {
			this.geneSampleTotal=geneSampleTotal;
		}

		public String toString(){
			String rtrn=geneName+"\t"+window.toUCSC(window.getOrientation())+"\t"+getNormalizedSampleScore()+"\t"+this.getPercentileNormalizedInputScore()+"\t"+getGlobalEnrichment()+"\t"+this.getPValue()+"\t"+this.getScanPValue(window.size(), totalGeneSizes)+"\t"+sampleCount+"\t"+inputCount+"\t"+getMaxInputCount()+"\t"+getWindowNormEnrichment()+"\t"+getWindowPVal()+"\t"+this.getLocalPVal()+"\t"+getLocalScanPValue(window.size());	
			//String rtrn=geneName+"\t"+window.toUCSC(window.getOrientation())+"\t"+sampleCount+"\t"+inputCount+"\t"+getMaxInputCount()+"\t"+getGlobalEnrichment()+"\t"+this.getPValue()+"\t"+this.getScanPValue(window.size(), totalGeneSizes);
			
			return rtrn;
		}
		
		

		public void setNumberSampleWindows(int clipWindows) {
			numberSampleWindows=clipWindows;
		}

		public void setNumberInputWindows(int inputWindows) {
			numberInputWindows=inputWindows;
		}
		
		public int getNumberSampleWindows(){
			return this.numberSampleWindows;
		}
		
		public int getNumberInputWindows(){
			return this.numberInputWindows;
		}
		
		public double getWindowNormEnrichment(){
			double numerator=(double)getSampleReadCount()/(double)getMaxInputCount();
			double elution=(double)this.getSampleTotalCount()/(double)this.getNumberSampleWindows();
			double input=(double)this.getInputTotalCount()/(double)this.getNumberInputWindows();
			double denominator=elution/input;
			return numerator/denominator;
		}
		

		public double getScanPValue(int windowSize, int totalSize) {
			return ScanStat.getPValue(getMaxInputCount(), getSampleReadCount(), inputTotal, sampleTotal, windowSize, totalSize);
		}
		
		
		public double getLocalScanPValue(int windowSize) {
			return ScanStat.getPValue(getMaxInputCount(), getSampleReadCount(), this.geneInputTotal, this.geneSampleTotal, windowSize, this.geneLength);
		}

		/**
		 * P values based on local "gene" level
		 * @return window within gene p-value
		 */
		public double getLocalPVal() {
			int k=this.getSampleReadCount();
			int n=getSampleReadCount()+getMaxInputCount();
			
			if(n>0){
				double elution=(double)this.geneSampleTotal;
				double input=(double) this.geneInputTotal;
				double p=elution/(elution+input);
				//System.err.println(elution+" "+input+" "+p);
				if(elution==0 && input==0){return 1.0;}
				if(p==1.0){return 1.0;}
				Binomial b=new Binomial(n, p);
				return 1-b.cdf(k);
			}
			return 1.0;
			
		}
		
		
		
		/**
		 * P values based on window normalization
		 * @return window normalized p-value
		 */
		public double getWindowPVal() {
			int k=this.getSampleReadCount();
			int n=getSampleReadCount()+getMaxInputCount();
			
			if(n>0){
				double elution=(double)this.getSampleTotalCount()/(double)this.getNumberSampleWindows();
				double input=(double) this.getInputTotalCount()/(double)this.getNumberInputWindows();
				double p=elution/(elution+input);
				if(elution ==0 && input ==0){return 1.0;}
				Binomial b=new Binomial(n, p);
				return 1-b.cdf(k);
			}
			return 1.0;
			
		}

		/**
		 * 
		 * @return Return the maximum of the input percentile of the value in this window
		 */
		public int getMaxInputCount() {
			int max=Math.max(this.inputPercentile, getInputReadCount());
			return max;
		}

		public double getGlobalEnrichment() {
			return getNormalizedSampleScore()/getPercentileNormalizedInputScore();
		}

		
		private double getPercentileNormalizedInputScore() {
			double num=(double)getMaxInputCount();
			return 1000000*(num/(double)inputTotal);
		}

		public double getNormalizedSampleScore() {
			int numerator=getSampleReadCount();
			return 1000000*((double)numerator/(double)sampleTotal);
		}

		public void setInputPercentile(int inputPercentile) {
			this.inputPercentile=inputPercentile;
			this.hasPercentile=true;
		}

		public double getPValue() {
			//Return binomial p
			
			int n=getSampleReadCount()+getInputReadCount();
			if(n>0){
				double p=(double)getSampleTotalCount()/((double)getInputTotalCount()+(double)getSampleTotalCount());
				Binomial b=new Binomial(n, p);
				return 1-b.cdf(getSampleReadCount());
			}
			return 1.0;
		}

		private int getInputTotalCount() {
			return this.inputTotal;
		}

		private int getSampleTotalCount() {
			return this.sampleTotal;
		}

		

		public int getSampleReadCount(){
			int count=this.sampleCount;
			if(count==0){count=1;}
			return count;
		}
		
		public int getInputReadCount(){
			int count= this.inputCount;
			if(count==0){count=1;}
			return count;
		}	

		public Annotation getParentAnnotation() {
			return this.parentAnnotation;
		}

		public void setParentAnnotation(Annotation region) {
			this.parentAnnotation=region;
		}

		public void setGeneName(String name) {
			this.geneName=name;
			
		}

		public void setGeneLength(int size) {
			this.geneLength=size;
			
		}
		
	}


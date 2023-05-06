package guttmanlab.core.splicing.speckle;

import guttmanlab.core.annotation.Gene;

public class JunctionScore {

	private Gene junction;
	private int exonJunctionCounts;
	private int exonCounts;
	private int intronCounts;
	
	public JunctionScore(Gene junction) {
		this.junction=junction;
	}
	
	public void setExonJunctionCounts(int count) {
		this.exonJunctionCounts=count;
	}
	
	public void setExonCounts(int count) {
		this.exonCounts=count;
	}
	
	public void setIntronCounts(int count) {
		this.intronCounts=count;
	}
	
	public double normalizedExonFraction() {
		double num=exonCounts+exonJunctionCounts;
		double denom=junction.size();
		return (num/denom);
	}
	
	public double normalizedIntronFraction() {
		double num=intronCounts;
		double denom=junction.getGenomicLength()-junction.size();
		return num/denom;
	}

	public void increment(int exonJunctionCounts, int exonCounts, int intronCounts) {
		this.exonJunctionCounts+=exonJunctionCounts;
		this.exonCounts+=exonCounts;
		this.intronCounts+=intronCounts;
	}

	public int getExonJunctionCounts() {
		return this.exonJunctionCounts;
	}

	public int getAllExonCounts() {
		return this.exonCounts+this.exonJunctionCounts;
	}

	public int getIntronCounts() {
		return this.intronCounts;
	}
	
}

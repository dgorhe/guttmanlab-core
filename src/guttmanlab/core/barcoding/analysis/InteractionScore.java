package guttmanlab.core.barcoding.analysis;

import java.util.Collection;
import java.util.Set;

import guttmanlab.core.annotation.Score;

public class InteractionScore implements Score {

	Collection<String> sharedBarcodes;
	Collection<String> barcodes1;
	Collection<String> barcodes2;
	double score;
	static final InteractionScore self=new InteractionScore(-1.0);
	
	public InteractionScore(Collection<String> sharedBarcodes, Collection<String> barcodes1, Collection<String> barcodes2) {
		this.sharedBarcodes=sharedBarcodes;
		this.barcodes1=barcodes1;
		this.barcodes2=barcodes2;
		this.score=sharedBarcodes.size();
	}
	
	public InteractionScore(double score) {
		this.score=score;
	}
	
	public Collection<String> getSharedBarcodes(){return this.sharedBarcodes;}
	
	public double normalizedScore(int numberOfBins1, int numberOfBins2){
		double prob1=(double)barcodes1.size()/(double)numberOfBins1;
		double prob2=(double)barcodes2.size()/(double)numberOfBins2;
		
		double expected=prob1*prob2;
		double observed=sharedBarcodes.size();
		return observed/expected;
	}
	
	@Override
	public double getScore() {
		return score;
	}

	@Override
	public boolean isSignificant() {
		// TODO Auto-generated method stub
		return false;
	}

}

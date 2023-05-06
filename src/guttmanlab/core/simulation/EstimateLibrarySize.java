package guttmanlab.core.simulation;

import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.math.Statistics;

public class EstimateLibrarySize {

	int numberORFs=250;
	int numberOfIterations=1;
	int numPerm=1000;
	
	public EstimateLibrarySize(){
		
		double[] perms=new double[numPerm];
		for(int j=0; j<numPerm; j++){
			Collection<Integer> nums=new TreeSet<Integer>();
			for(int i=0; i<numberOfIterations; i++){
				nums.addAll(sample(numberORFs));
			}
			double percent=(double)nums.size()/(double)numberORFs;
			perms[j]=percent;
			System.err.println(j+" "+nums.size()+" "+numberORFs+" "+percent);
		}
		
		System.err.println(numberOfIterations+" mean "+Statistics.mean(perms));
		
	}

	private Collection<Integer> sample(int n) {
		Collection<Integer> rtrn=new TreeSet<Integer>();
		for(int i=0; i<n; i++){
			Integer val=new Double((Math.random()*n)).intValue();
			rtrn.add(val);
		}
		return rtrn;
	}
	
	public static void main(String[] args){
		new EstimateLibrarySize();
	}
	
}

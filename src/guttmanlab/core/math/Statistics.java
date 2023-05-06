package guttmanlab.core.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.apache.commons.math.distribution.HypergeometricDistributionImpl;

public class Statistics {

	/**
	 * Uses the Weighted average method.
	 * @param list - Ordered list of numbers (i.e. short, float, double, etc)
	 * @param pct  - Desired quantile
	 * @return the estimated quantile requested: x s.t. P(X <= x) >= pct
	 */
	public static double quantile(List<Double> list, double pct) {
		if(list.size() == 1) { return list.get(0).doubleValue();}
		if(list.size() == 0) { return 0;}
		Collections.sort(list);
		if(pct==0.0){
			return (list.get(0).doubleValue());
		}
		if(pct==1.0){
			return list.get(list.size()-1).doubleValue();
		}
		int idx = (int)  Math.floor( pct * (list.size() - 1));
		double reminder = pct * (list.size() - 1) - idx;
		double idxthTerm = list.get(idx).doubleValue();
		double idxNextthTerm = list.get(idx + 1).doubleValue();
		//System.out.println("pct " + pct + " # " + list.size() + " reminder " + reminder + " idx " + idx + " idxthTerm " + idxthTerm + " idxNextthTerm " + idxNextthTerm);
		return  idxthTerm + reminder*(idxNextthTerm - idxthTerm);
	}
	
	
	public static int quantileInt(List<Integer> list, double pct) {
		if(list.size() == 1) { return list.get(0);}
		if(list.size() == 0) { return 0;}
		Collections.sort(list);
		if(pct==0.0){
			return (list.get(0));
		}
		if(pct==1.0){
			return list.get(list.size()-1);
		}
		
		
		int idx = (int)  Math.ceil( pct * (list.size() - 1));
		int idxthTerm = list.get(idx);
		return  idxthTerm;
	}

	public static double quantile(int [] vals, double pct) {
		if(vals.length == 1) { return vals[0];}
		if(vals.length == 0) { return 0;}
		Arrays.sort(vals);
		int idx = (int)  Math.floor( pct * (vals.length - 1));
		double reminder = pct * (vals.length - 1) - idx;
		double idxthTerm = vals[idx];
		double idxNextthTerm = vals[idx + 1];
		//System.out.println("pct " + pct + " # " + list.size() + " reminder " + reminder + " idx " + idx + " idxthTerm " + idxthTerm + " idxNextthTerm " + idxNextthTerm);
		return  idxthTerm + reminder*(idxNextthTerm - idxthTerm);
	}
	
	public static double quantile(double [] vals, double pct) {
		if(vals.length == 1) { return vals[0];}
		if(vals.length == 0) { return 0;}
		Arrays.sort(vals);
		int idx = (int)  Math.floor( pct * (vals.length - 1));
		double reminder = pct * (vals.length - 1) - idx;
		double idxthTerm = vals[idx];
		double idxNextthTerm = vals[idx + 1];
		//System.out.println("pct " + pct + " # " + list.size() + " reminder " + reminder + " idx " + idx + " idxthTerm " + idxthTerm + " idxNextthTerm " + idxNextthTerm);
		return  idxthTerm + reminder*(idxNextthTerm - idxthTerm);
	}

	public static double fisherExact(int a, int b, int c, int d){
		int populationSize=a+b+c+d;
		int numberOfSuccesses=a+b;
		int sampleSize=a+c; //num of draws
		
		if(numberOfSuccesses ==0 || sampleSize==0){return 1.0;}
		
		HypergeometricDistributionImpl dist=new HypergeometricDistributionImpl(populationSize, numberOfSuccesses, sampleSize);
		
		return 1-dist.cumulativeProbability(a);
	}
	
	
	public static double fisherExact(double a, double b, double c, double d){
		int populationSize=(int)(a+b+c+d);
		int numberOfSuccesses=(int)(a+b);
		int sampleSize=(int)(a+c); //num of draws
		
		if(numberOfSuccesses ==0 || sampleSize==0){return 1.0;}
		
		HypergeometricDistributionImpl dist=new HypergeometricDistributionImpl(populationSize, numberOfSuccesses, sampleSize);
		
		return 1-dist.cumulativeProbability((int)a);
	}
	
	public static double mean(double [] values) {
		double total = 0;
		int counter=0;
		for (int i = 0; i < values.length; i++) {
			if(!new Double(values[i]).equals(Double.NaN)){
				total = total + values[i];
				counter++;
			}
		}
		
		return total/counter;
	}
	
	public static double mean(double [] values, double val) {
		double total = val;
		int counter=1;
		for (int i = 0; i < values.length; i++) {
			if(!new Double(values[i]).equals(Double.NaN)){
				total = total + values[i];
				counter++;
			}
		}
		
		return total/counter;
	}
	
	public static double mean(Collection<? extends Number>  values) {
		double total = 0;
		int size = values.size();
		for (Number n: values) {
			total = total + n.doubleValue();
		}
		
		return total/(double)size;
	}
	
	
	public static double sum(Collection<? extends Number>  values) {
		double total = 0;
		int size = values.size();
		for (Number n: values) {
			total = total + n.doubleValue();
		}
		
		return total;
	}
	
	public static double min(Collection<? extends Number>  values) {
		double min=Double.MAX_VALUE;
		for (Number n: values) {
			min=Math.min(min, n.doubleValue());
		}
		
		return min;
	}
	
	public static double max(Collection<? extends Number>  values) {
		double max=-Double.MAX_VALUE;
		for (Number n: values) {
			max=Math.max(max, n.doubleValue());
		}
		
		return max;
	}

	public static double percentLessThan(double observed, double[] random) {
		int count=0;
		for(int i=0; i<random.length; i++){
			if(random[i]<=observed){count++;}
		}
		return (double)count/(double)random.length;
	}
	
	public static double percentLessThan(int observed, int[] random) {
		int count=0;
		for(int i=0; i<random.length; i++){
			if(random[i]<=observed){count++;}
		}
		return (double)count/((double)random.length);
	}
	
	public static double percentGreaterThan(int observed, int[] random) {
		int count=0;
		for(int i=0; i<random.length; i++){
			if(random[i]>=observed){count++;}
		}
		return (double)count/((double)random.length);
	}
	
	public static double percentLessThan(double observed, List<Double> random) {
		int count=0;
		for(int i=0; i<random.size(); i++){
			if(random.get(i)<=observed){count++;}
		}
		return (double)count/(double)random.size();
	}
	
	
	public static double percentGreaterThan(double observed, double[] random) {
		int count=0;
		for(int i=0; i<random.length; i++){
			if(random[i]>=observed){count++;}
		}
		return (double)count/(double)random.length;
	}

	public static double max(double[] scores) {
		double max=-Double.MAX_VALUE;
		for(int i=0; i<scores.length; i++){
			max=Math.max(max, scores[i]);
		}
		return max;
	}
	
	public static int max(int[] scores) {
		int max=-Integer.MAX_VALUE;
		for(int i=0; i<scores.length; i++){
			max=Math.max(max, scores[i]);
		}
		return max;
	}
	
	
	public static double max(double[] scores, double exclude) {
		double max=-Double.MAX_VALUE;
		for(int i=0; i<scores.length; i++){
			if(scores[i]!=exclude) {
				max=Math.max(max, scores[i]);
			}
		}
		return max;
	}
	
	public static double min(double[] scores) {
		double min=Double.MAX_VALUE;
		for(int i=0; i<scores.length; i++){
			min=Math.min(min, scores[i]);
		}
		return min;
	}

	public static double zScore(double actual, ArrayList<Double> otherScores) {
		double mean=mean(otherScores);
		double sem=sem(otherScores, mean);
		return (actual-mean)/sem;
	}

	
	public static double zScore(double actual, double[] otherScores) {
		double mean=mean(otherScores);
		double sem=sem(otherScores, mean);
		return (actual-mean)/sem;
	}
	
	
	public static double sem(double[] otherScores, double mean) {
		double sum=0;
		for(Double val: otherScores){
			sum+=Math.pow(val-mean,2);
		}
		return Math.sqrt(sum/otherScores.length);
	}
	
	public static double sem(List<Double> otherScores, double mean) {
		double sum=0;
		for(Double val: otherScores){
			sum+=Math.pow(val-mean,2);
		}
		return Math.sqrt(sum/otherScores.size());
	}

	public static double meanWithoutZero(double[] values) {
		double total = 0;
		int counter=0;
		for (int i = 0; i < values.length; i++) {
			if(values[i]>0){
				counter++;
				total=total+values[i];
			}
		
		}
		
		return total/counter;
	}

	public static double sum(double[] vals) {
		double sum=0;
		for(int i=0; i<vals.length; i++){
			sum+=vals[i];
		}
		return sum;
	}


	public static double max(List<Double> list) {
		double max=list.iterator().next();
		for(Double val: list){
			max=Math.max(val, max);
		}
		return max;
	}


	public static double maxAbs(double[] vals) {
		double max=0;
		for(Double val: vals){
			max=Math.max(Math.abs(val), max);
		}
		return max;
	}


	public static double mean(int[] vals) {
		double sum=0;
		
		for(int i=0; i<vals.length; i++) {sum+=vals[i];}
		
		return sum/(double)vals.length;
	}

	
	public static int sum(int[] vals) {
		
		int sum=0;
		for(int i=0; i<vals.length; i++){
			sum+=vals[i];
		}
		return sum;
		
		
	}

	public static int[] sum(int[] list, int[] list2) {
		int[] rtrn=new int[list.length];
		
		for(int i=0; i<list.length; i++) {
			rtrn[i]=list[i]+list2[i];
		}
		
		return rtrn;
	}
}

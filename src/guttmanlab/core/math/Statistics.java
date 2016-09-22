package guttmanlab.core.math;

import java.util.Collection;
import java.util.List;

public class Statistics {

	/**
	 * Uses the Weighted average method.
	 * @param list - Ordered list of numbers (i.e. short, float, double, etc)
	 * @param pct  - Desired quantile
	 * @return the estimated quantile requested: x s.t. P(X <= x) >= pct
	 */
	public static double quantile(List<? extends Number> list, double pct) {
		if(list.size() == 1) { return list.get(0).doubleValue();}
		if(list.size() == 0) { return 0;}
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
	
	/**
	 * Get a quantile of a sorted array
	 * @param vals Sorted values
	 * @param pct Quantile to get
	 * @return The value
	 */
	public static double quantile(double [] vals, double pct) {
		if(vals.length == 1) { return vals[0];}
		if(vals.length == 0) { return 0;}
		int idx = (int)  Math.floor( pct * (vals.length - 1));
		double reminder = pct * (vals.length - 1) - idx;
		double idxthTerm = vals[idx];
		double idxNextthTerm = vals[idx + 1];
		//System.out.println("pct " + pct + " # " + list.size() + " reminder " + reminder + " idx " + idx + " idxthTerm " + idxthTerm + " idxNextthTerm " + idxNextthTerm);
		return  idxthTerm + reminder*(idxNextthTerm - idxthTerm);
	}

	/**
	 * Get the mean of an array of values
	 * @param values Values
	 * @return Mean value
	 */
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
	
	/**
	 * Get the mean of a collection of values
	 * @param values Values
	 * @return Mean value
	 */
	public static double mean(Collection<? extends Number>  values) {
		double total = 0;
		int size = values.size();
		for (Number n: values) {
			total = total + n.doubleValue();
		}
		
		return total/(double)size;
	}
}

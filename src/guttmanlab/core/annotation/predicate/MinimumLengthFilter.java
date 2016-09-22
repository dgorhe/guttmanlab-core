package guttmanlab.core.annotation.predicate;

import guttmanlab.core.annotation.Annotation;

import org.apache.commons.collections15.Predicate;

public class MinimumLengthFilter<T extends Annotation> implements Predicate<T>{

	int minSize;
	
	/**
	 * @param minSize Minimum feature size
	 */
	public MinimumLengthFilter(int minSize){
		this.minSize=minSize;
	}
	
	@Override
	public boolean evaluate(T annot) {
		if(annot.size()>minSize){return true;}
		return false;
	}

}

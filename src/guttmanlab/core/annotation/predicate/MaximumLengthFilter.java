package guttmanlab.core.annotation.predicate;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;

import org.apache.commons.collections15.Predicate;

public class MaximumLengthFilter<T extends Annotation> implements Predicate<T>{

	int maxSize;
	
	/**
	 * @param maxSize Maximum feature size
	 */
	public MaximumLengthFilter(int maxSize){
		this.maxSize=maxSize;
	}
	
	@Override
	public boolean evaluate(T annot) {
		if(annot.getReferenceEndPosition()-annot.getReferenceStartPosition() < maxSize){return true;}
		return false;
	}

}

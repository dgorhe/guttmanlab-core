package guttmanlab.core.annotation.predicate;

import org.apache.commons.collections15.Predicate;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.MappedFragment;

/**
 * A predicate which evaluates whether a MappedFragment is a subset of a particular Annotation
 */
public class ContainedByFilter<T extends MappedFragment> implements Predicate<T> {

	private Annotation annot;
	
	public ContainedByFilter(Annotation annot) {
		this.annot = annot;
	}
	
	@Override
	public boolean evaluate(T fragment) {
		return annot.contains(fragment);
	}
}
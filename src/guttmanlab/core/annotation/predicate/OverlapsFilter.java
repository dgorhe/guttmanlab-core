package guttmanlab.core.annotation.predicate;

import org.apache.commons.collections15.Predicate;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.MappedFragment;

/**
 * A Predicate which evaluates whether a MappedFragment overlaps a given Annotation.
 */
public class OverlapsFilter<T extends MappedFragment> implements Predicate<T> {

	private Annotation annot;
	
	/**
	 * @param annot Annotation that fragments are required to overlap
	 */
	public OverlapsFilter(Annotation annot) {
		this.annot = annot;
	}
	
	@Override
	public boolean evaluate(T fragment) {
		return fragment.overlaps(annot);
	}
}
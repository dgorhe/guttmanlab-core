package guttmanlab.core.annotation.predicate;

import guttmanlab.core.annotation.MappedFragment;

import org.apache.commons.collections15.Predicate;


public class SAMFragmentNumHitsFilter implements Predicate<MappedFragment> {

	private int maxNumHits;
	
	/**
	 * @param maxHits Max number of hits allowed
	 */
	public SAMFragmentNumHitsFilter(int maxHits) {
		maxNumHits = maxHits;
	}
	
	@Override
	public boolean evaluate(MappedFragment fragment) {
		return fragment.getNumHits() <= maxNumHits;
	}

}

package guttmanlab.core.annotation.predicate;

import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;

import org.apache.commons.collections15.Predicate;

public class InsertSizeFilter<T extends PairedMappedFragment<SAMFragment>> implements Predicate<T>{

	private int maxSize;
	
	/**
	 * @param maxSize Max insert size between mates
	 */
	public InsertSizeFilter(int maxSize){
		this.maxSize=maxSize;
	}
	
	@Override
	public boolean evaluate(T fragment) {
		if(fragment.fragmentLength()<maxSize){return true;}
		return false;
	}

}

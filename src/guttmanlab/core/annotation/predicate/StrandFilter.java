package guttmanlab.core.annotation.predicate;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;

import org.apache.commons.collections15.Predicate;

public class StrandFilter<T extends Annotation> implements Predicate<T> {

	Strand feat_dir;
	
	/**
	 * @param f Required strand
	 */
	public StrandFilter(Strand f){
		this.feat_dir = f;
	}
	
	@Override
	public boolean evaluate(Annotation read) {
		//if read unknown or invalid, false
		Strand read_dir = read.getOrientation();
		if(read_dir.equals(Strand.BOTH)) {
			return true;
		}
		if(read_dir.equals(Strand.INVALID) || read_dir.equals(Strand.UNKNOWN)) {
			return false;
		}
		//if read positive or negative
		if(read_dir.equals(Strand.POSITIVE)) {
			return (feat_dir.equals(Strand.BOTH) || feat_dir.equals(Strand.POSITIVE));
		}
		if(read_dir.equals(Strand.NEGATIVE)) {
			return (feat_dir.equals(Strand.BOTH) || feat_dir.equals(Strand.NEGATIVE));
		}
		return false;
	}

}

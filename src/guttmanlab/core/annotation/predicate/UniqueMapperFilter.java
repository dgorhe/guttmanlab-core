package guttmanlab.core.annotation.predicate;

import org.apache.commons.collections15.Predicate;

import guttmanlab.core.annotation.SAMFragment;

public class UniqueMapperFilter implements Predicate<SAMFragment> {

	@Override
	public boolean evaluate(SAMFragment read) {
		if(read.getMappingQuality()==255){return true;}
		return false;
	}

}

package guttmanlab.core.annotation.predicate;

import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;

import org.apache.commons.collections15.Predicate;

public class SecondReadFilter<T extends PairedMappedFragment<SAMFragment>> implements Predicate<T> {
		
	public boolean evaluate(T frag) {
		return frag.getSamRecord(null).getSecondOfPairFlag();
	}

}

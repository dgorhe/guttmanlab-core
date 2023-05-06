package guttmanlab.core.barcoding.analysis;

import java.util.Collection;

import org.apache.commons.collections15.Predicate;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SingleInterval;

public class RegionPredicate implements Predicate<Cluster> {
	
	Collection<Annotation> regions;

	public RegionPredicate(Collection<Annotation> regions){}
	
	@Override
	public boolean evaluate(Cluster cluster) {
		for(SingleInterval interval: cluster.getAllIntervals()){
			for(Annotation region: regions){
				if(interval.overlaps(region)){return true;}
			}
		}
		return false;
	}

}

package guttmanlab.core.annotation;

import guttmanlab.core.annotation.predicate.ReadFlag;

import java.util.Collection;

/**
 * An Abstract Class representing functionality specific to mapped reads
 * @author mguttman
 *
 */
public interface MappedFragment extends Annotation{

	public Collection<? extends ReadFlag> getReadFlags();
	
	public int getNumHits();
	
	public int getMappingQuality();
	
}

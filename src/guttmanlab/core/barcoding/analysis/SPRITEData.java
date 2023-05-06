package guttmanlab.core.barcoding.analysis;

import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.Annotation;
import htsjdk.samtools.util.CloseableIterator;

public interface SPRITEData extends CloseableIterator<Cluster>{

	public Collection<Cluster> getClustersOverlappingRegion(Annotation region) throws IOException;
	
}

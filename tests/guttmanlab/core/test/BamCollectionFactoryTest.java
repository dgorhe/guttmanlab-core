package guttmanlab.core.test;

import static org.junit.Assert.*;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotationcollection.AbstractAnnotationCollection;
import guttmanlab.core.annotationcollection.BAMFragmentCollectionFactory;
import net.sf.samtools.util.CloseableIterator;

import org.junit.Test;

public class BamCollectionFactoryTest {

	private static String singleBam = "/Users/cburghard/Downloads/Ribosome.bam";
	private static String pairedBam = "/Users/cburghard/Downloads/input.bam";
	
	@Test
	public void testSingleDetection() {
		assertTrue(!BAMFragmentCollectionFactory.isPairedEnd(singleBam));
	}

	@Test
	public void testPairedDetection() {
		assertTrue(BAMFragmentCollectionFactory.isPairedEnd(pairedBam));
	}
	
	@Test
	public void testSingleImplementation() {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(singleBam);
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		iter.hasNext();
		assertEquals("annotation.SAMFragment", iter.next().getClass().getName());
	}
	
	@Test
	public void testForceSingleImplementation() {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam, true);
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		iter.hasNext();
		assertEquals("annotation.SAMFragment", iter.next().getClass().getName());
	}
	
	@Test
	public void testPairedImplementation() {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam);
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		iter.hasNext();
		assertEquals("annotation.PairedMappedFragment", iter.next().getClass().getName());
	}
	
}

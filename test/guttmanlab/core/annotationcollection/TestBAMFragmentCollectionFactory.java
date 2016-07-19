package guttmanlab.core.annotationcollection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import guttmanlab.core.annotation.MappedFragment;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.util.CloseableIterator;

import org.junit.Test;

public class TestBAMFragmentCollectionFactory {
	
	private static File singleBam = new File("test/resources/SingleCollectionTest.bam");
	private static File pairedBam = new File("test/resources/PairedCollectionTest.bam");	

	@Test
	public void testSingleDetection() {
		assertTrue(!BAMFragmentCollectionFactory.isPairedEnd(singleBam));
	}

	@Test
	public void testPairedDetection() {
		assertTrue(BAMFragmentCollectionFactory.isPairedEnd(pairedBam));
	}

	@Test
	public void testSingleImplementation() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(singleBam);
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		assertEquals("guttmanlab.core.annotation.SAMFragment", iter.next().getClass().getName());
		iter.close();
	}
	
	@Test
	public void testForceSingleImplementation() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam, true);
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		assertEquals("guttmanlab.core.annotation.SAMFragment", iter.next().getClass().getName());
		iter.close();
	}
	
	@Test
	public void testPairedImplementation() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam);
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		assertEquals("guttmanlab.core.annotation.PairedMappedFragment", iter.next().getClass().getName());
		iter.close();
	}

	
}

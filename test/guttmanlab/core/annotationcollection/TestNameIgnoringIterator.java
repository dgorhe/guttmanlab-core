package guttmanlab.core.annotationcollection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.util.CloseableIterator;

import org.junit.Before;
import org.junit.Test;

public class TestNameIgnoringIterator {
	
	private static final CoordinateSpace refSpace = CoordinateSpace.MM9;
	
	private static final File bedFile = new File("test/resources/iter_ignore_name_test.bed");
	private static AnnotationCollection<BEDFileRecord> genes;
	
	private static final File singleEndBamFileSorted = new File("test/resources/SingleCollectionTest_iter_ignore_name.sorted.bam");
	private static final File pairedEndBamFileSorted = new File("test/resources/PairedCollectionTest_iter_ignore_name.sorted.bam");
	private static final File singleEndBamFileUnsorted = new File("test/resources/SingleCollectionTest_iter_ignore_name.unsorted.bam");
	private static final File pairedEndBamFileUnsorted = new File("test/resources/PairedCollectionTest_iter_ignore_name.unsorted.bam");
	private AnnotationCollection<? extends MappedFragment> singleReadsSorted;
	private AnnotationCollection<? extends MappedFragment> pairedReadsSorted;
	private AnnotationCollection<? extends MappedFragment> singleReadsUnsorted;
	private AnnotationCollection<? extends MappedFragment> pairedReadsUnsorted;
	
	@Before
	public void setUp() {
		try {
			genes = BEDFileIO.loadFromFile(bedFile, refSpace);
			singleReadsSorted = BAMFragmentCollectionFactory.createFromBam(singleEndBamFileSorted);
			pairedReadsSorted = BAMFragmentCollectionFactory.createFromBam(pairedEndBamFileSorted);
			singleReadsUnsorted = BAMFragmentCollectionFactory.createFromBam(singleEndBamFileUnsorted);
			pairedReadsUnsorted = BAMFragmentCollectionFactory.createFromBam(pairedEndBamFileUnsorted);
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	private static int iterCount(CloseableIterator<?> iter) throws IllegalStateException {
		int num = 0;
		while(iter.hasNext()) {
			num++;
			iter.next();
		}
		iter.close();
		return num;
	}
	
	@Test
	public void testBed() {
		CloseableIterator<BEDFileRecord> iter = genes.sortedIteratorIgnoreName();
		assertEquals(40, iterCount(iter));
	}
	
	
	@Test
	public void testSingleEndBamSorted() {
		CloseableIterator<? extends MappedFragment> iter = singleReadsSorted.sortedIteratorIgnoreName();
		assertEquals(7, iterCount(iter));
	}
	
	@Test
	public void testPairedEndBamSorted() {
		/**
		 * The BAMPairedFragmentCollection iterator only returns fragments with both mates mapped
		 * To see why the expected number is 8, run
		 * samtools view -f 129 -F 12 PairedCollectionTest_iter_ignore_name.sorted.bam | awk '{print $3 "\t" $4 "\t" $6}' | sort -u
		 */
		CloseableIterator<? extends MappedFragment> iter = pairedReadsSorted.sortedIteratorIgnoreName();
		assertEquals(8, iterCount(iter));
	}
	
	/**
	 * Fails if the SAMRecordIterator fails to throw an appropriate IllegalStateException for unsorted data
	 * @param data
	 */
	private static void testUnsorted(AnnotationCollection<? extends MappedFragment> data) {
		CloseableIterator<? extends MappedFragment> iter = data.sortedIteratorIgnoreName();
		boolean correctExceptionThrown = false;
		try {
			@SuppressWarnings("unused")
			int c = iterCount(iter);
		} catch(IllegalStateException e) {
			String message = e.getMessage();
			if(message.contains("should come after") && message.contains("when sorting with")) {
				correctExceptionThrown = true;
			} else {
				throw e;
			}
		}
		assertTrue("An IllegalStateException with message \"Records ... should come after ... when sorting with ...\" "
				+ "should have been thrown for unsorted bam file.", correctExceptionThrown);
	}
	
	@Test
	public void testSingleEndBamUnsorted() {
		testUnsorted(singleReadsUnsorted);
	}
	
	@Test
	public void testPairedEndBamUnsorted() {
		testUnsorted(pairedReadsUnsorted);
	}
	

}

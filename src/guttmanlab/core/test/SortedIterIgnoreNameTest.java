package guttmanlab.core.test;

import static org.junit.Assert.*;

import java.io.IOException;

import net.sf.samtools.util.CloseableIterator;

import org.junit.Before;
import org.junit.Test;

import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMFragmentCollectionFactory;
import guttmanlab.core.coordinatespace.CoordinateSpace;

public class SortedIterIgnoreNameTest {
	
	private static final String refSpaceFile = "/storage/shared/CoreTestData/sizes";
	private static final CoordinateSpace refSpace = new CoordinateSpace(refSpaceFile);
	
	private static final String bedFile = "/storage/shared/CoreTestData/iter_ignore_name_test.bed";
	private AnnotationCollection<BEDFileRecord> genes;
	
	private static final String singleEndBamFileSorted = "/storage/shared/CoreTestData/SingleCollectionTest_iter_ignore_name.sorted.bam";
	private static final String pairedEndBamFileSorted = "/storage/shared/CoreTestData/PairedCollectionTest_iter_ignore_name.sorted.bam";
	private static final String singleEndBamFileUnsorted = "/storage/shared/CoreTestData/SingleCollectionTest_iter_ignore_name.unsorted.bam";
	private static final String pairedEndBamFileUnsorted = "/storage/shared/CoreTestData/PairedCollectionTest_iter_ignore_name.unsorted.bam";
	private AnnotationCollection<? extends MappedFragment> singleReadsSorted;
	private AnnotationCollection<? extends MappedFragment> pairedReadsSorted;
	@SuppressWarnings("unused")
	private AnnotationCollection<? extends MappedFragment> singleReadsUnsorted;
	@SuppressWarnings("unused")
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
	@SuppressWarnings("unused")
	private static void testUnsorted(AnnotationCollection<? extends MappedFragment> data) {
		CloseableIterator<? extends MappedFragment> iter = data.sortedIteratorIgnoreName();
		boolean correctExceptionThrown = false;
		try {
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
	
//	@Test
//	public void testSingleEndBamUnsorted() {
//		testUnsorted(singleReadsUnsorted);
//	}
//	
//	@Test
//	public void testPairedEndBamUnsorted() {
//		testUnsorted(pairedReadsUnsorted);
//	}
	
	
	
	
}

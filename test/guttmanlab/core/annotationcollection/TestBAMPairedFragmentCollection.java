package guttmanlab.core.annotationcollection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.predicate.MaximumLengthFilter;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.coordinatespace.GenomeSize;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;

import org.junit.BeforeClass;
import org.junit.Test;

public class TestBAMPairedFragmentCollection {
	
	private static File pairedBam = new File("test/resources/PairedCollectionTest.bam");	
	private static File chr19bam = new File("test/resources/chr19.clean.sorted.bam");
	private static SingleInterval malat1 = new SingleInterval("chr19", 5795689, 5802671, Strand.BOTH, "Malat1");
	private static AnnotationCollection<BEDFileRecord> refSeqFeatures;	

	@BeforeClass
	public static void setUp() throws IOException
	{
		BEDFileIO io =  new BEDFileIO(new CoordinateSpace(GenomeSize.MM9));
		refSeqFeatures = io.loadFromFile(new File("test/resources/RefSeq_mm9.bed"));
	}

	@Test
	public void testPairedIteratorCount() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam);
		assertEquals(8, data.getNumAnnotations());
	}

	@Test(expected = UnsupportedOperationException.class)
	public void testPairedIteratorRemove() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam);
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		iter.remove();
		iter.close();
	}

	@Test
	public void testPairedInitialNext() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam);
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		assertTrue(iter.next() != null);
		iter.close();
	}

	@Test
	public void testPairedIteration() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam);
		CloseableIterator<? extends MappedFragment> iter = data.sortedIterator();
		for (int i = 0; i < 8; i++) {
			assertTrue(iter.hasNext());
			iter.next();			
		}
		assertFalse(iter.hasNext());
		iter.close();
	}

	@Test
	public void testPairedNumOverlappersPartial() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam);
		assertEquals(6, data.numOverlappers(malat1, false));
	}
	
	@Test
	public void testPairedNumOverlappersFull() throws IOException {
		AbstractAnnotationCollection<? extends MappedFragment> data = BAMFragmentCollectionFactory.createFromBam(pairedBam);
		assertEquals(4, data.numOverlappers(malat1, true));
	}

	@Test  //Pass, paired ends take a long time
	//Verifies that the correct number of reads are returned for Ccdc87, a single exon gene on the positive strand, using paired end reads
	public void ccdc87SPairedEndOverlapReadCount() throws IOException{
		
		BAMPairedFragmentCollection bamPair = new BAMPairedFragmentCollection(chr19bam);
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_207268"))
				break;
		}
		iter.close();
		
		int count =0;
		bamPair.addFilter(new MaximumLengthFilter<PairedMappedFragment<SAMFragment>>(10000));
		CloseableIterator<PairedMappedFragment<SAMFragment>> f_iter = bamPair.sortedIterator(a, true);
		while(f_iter.hasNext())
		{
			PairedMappedFragment<SAMFragment> f = f_iter.next();
			count++;
		}
		
		f_iter.close();
		assertEquals("20 unconverted reads should overlap Ccdc87.",20,count);
		
	}

	
	@Test 
	//Verifies that the correct number of reads are returned for Ccdc87, a single exon gene on the positive strand, using paired end reads
	public void ccdc87SPairedEndConvertCoordinates() throws IOException{
		
		BAMPairedFragmentCollection bamPair = new BAMPairedFragmentCollection(chr19bam);
		
		CloseableIterator<BEDFileRecord> iter = refSeqFeatures.sortedIterator();

		Annotation a = null;
		while(iter.hasNext()) 
		{
			a = iter.next();
			if(a.getName().equals("NM_207268"))
				break;
		}
		iter.close();
		
		int count =0;
		Annotation b = new SingleInterval(a.getName(), 0, a.size()-1);
		b.setOrientation(Strand.BOTH);
		AnnotationCollection<DerivedAnnotation<PairedMappedFragment<SAMFragment>>> converted = refSeqFeatures.convertCoordinates(bamPair, CoordinateSpace.MM9, true);	
		CloseableIterator<DerivedAnnotation<PairedMappedFragment<SAMFragment>>> c_iter = converted.sortedIterator(b,true);
		
		while(c_iter.hasNext())
		{
			DerivedAnnotation<PairedMappedFragment<SAMFragment>> c = c_iter.next();
			count++;
		}
		
		assertEquals("20 unconverted reads should overlap Ccdc87.",20,count);
		
	}
	

	
}

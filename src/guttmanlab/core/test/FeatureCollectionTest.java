package guttmanlab.core.test;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Before;
import org.junit.Test;

import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;

public class FeatureCollectionTest {

	private static final String refSpaceFile = "/storage/shared/CoreTestData/sizes";
	private static final CoordinateSpace refSpace = new CoordinateSpace(refSpaceFile);
	private static final String bedFile = "/storage/shared/CoreTestData/iter_ignore_name_test.bed";
	private FeatureCollection<BEDFileRecord> genes;

	@Before
	public void setUp() {
		try {
			genes = (FeatureCollection<BEDFileRecord>) BEDFileIO.loadFromFile(bedFile, refSpace);
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		genes.add(new BEDFileRecord(new SingleInterval("chrM", 25, 30, Strand.NEGATIVE)));
		genes.add(new BEDFileRecord(new SingleInterval("chr19", 61342420, 61342430, Strand.POSITIVE)));
	}

	@Test
	public void testGeneCount() {
		assertEquals(46, genes.getNumAnnotations());
	}
	
	@Test
	public void testOverlaps() {
		
		SingleInterval interval1 = new SingleInterval("chr2", 130196216, 130196217);
		assertTrue("Strand not specified", genes.overlaps(interval1));

		SingleInterval interval2 = new SingleInterval("chr2", 130196216, 130196217, Strand.POSITIVE);
		assertTrue("Same strand", genes.overlaps(interval2));

		SingleInterval interval3 = new SingleInterval("chr2", 130196216, 130196217, Strand.BOTH);
		assertTrue("Both strands", genes.overlaps(interval3));

		SingleInterval interval4 = new SingleInterval("chr2", 130196216, 130196217, Strand.NEGATIVE);
		assertTrue("Wrong strand", !genes.overlaps(interval4));

		SingleInterval interval5 = new SingleInterval("chr2", 130196216, 130196217, Strand.UNKNOWN);
		assertTrue("Strand unknown", genes.overlaps(interval5));

		SingleInterval interval6 = new SingleInterval("chr2", 130196216, 130196217, Strand.INVALID);
		assertTrue("Invalid strand", !genes.overlaps(interval6));

		SingleInterval interval7 = new SingleInterval("chr2", 130196217, 130276013, Strand.POSITIVE);
		assertTrue("No overlap plus strand", !genes.overlaps(interval7));

		SingleInterval interval8 = new SingleInterval("chr2", 130196217, 130276013, Strand.NEGATIVE);
		assertTrue("No overlap minus strand", !genes.overlaps(interval8));
		
		SingleInterval interval9 = new SingleInterval("chrM", 20, 30, Strand.NEGATIVE);
		assertTrue("Overlap chrM", genes.overlaps(interval9));
		
		SingleInterval interval10 = new SingleInterval("chrM", 10, 20, Strand.POSITIVE);
		assertTrue("No overlap chrM", !genes.overlaps(interval10));
		
		SingleInterval interval11 = new SingleInterval("chr7", 0, 5, Strand.POSITIVE);
		assertTrue("No overlap at beginning of chr7", !genes.overlaps(interval11));
		
		SingleInterval interval12 = new SingleInterval("chr19", 61342410, 61342425, Strand.POSITIVE);
		assertTrue("Overlap at end of chr19", genes.overlaps(interval12));
		
	}
	
	
	
	
	
	
	
	
	
	
	
}

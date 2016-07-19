package guttmanlab.core.annotationcollection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import org.junit.Before;
import org.junit.Test;

public class TestFeatureCollection {
	
	private FeatureCollection<BEDFileRecord> fewGenes;

	@Before
	public void setUp() {
		try {
			fewGenes = (FeatureCollection<BEDFileRecord>) BEDFileIO.loadFromFile(new File("test/resources/iter_ignore_name_test.bed"), CoordinateSpace.MM9);
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		fewGenes.add(new BEDFileRecord(new SingleInterval("chrM", 25, 30, Strand.NEGATIVE)));
		fewGenes.add(new BEDFileRecord(new SingleInterval("chr19", 61342420, 61342430, Strand.POSITIVE)));
	}

	@Test
	public void testGeneCount() {
		assertEquals(46, fewGenes.getNumAnnotations());
	}
	
	@Test
	public void testOverlaps() {
		
		SingleInterval interval1 = new SingleInterval("chr2", 130196216, 130196217);
		assertTrue("Strand not specified", fewGenes.overlaps(interval1));

		SingleInterval interval2 = new SingleInterval("chr2", 130196216, 130196217, Strand.POSITIVE);
		assertTrue("Same strand", fewGenes.overlaps(interval2));

		SingleInterval interval3 = new SingleInterval("chr2", 130196216, 130196217, Strand.BOTH);
		assertTrue("Both strands", fewGenes.overlaps(interval3));

		SingleInterval interval4 = new SingleInterval("chr2", 130196216, 130196217, Strand.NEGATIVE);
		assertTrue("Wrong strand", !fewGenes.overlaps(interval4));

		SingleInterval interval5 = new SingleInterval("chr2", 130196216, 130196217, Strand.UNKNOWN);
		assertTrue("Strand unknown", fewGenes.overlaps(interval5));

		SingleInterval interval6 = new SingleInterval("chr2", 130196216, 130196217, Strand.INVALID);
		assertTrue("Invalid strand", !fewGenes.overlaps(interval6));

		SingleInterval interval7 = new SingleInterval("chr2", 130196217, 130276013, Strand.POSITIVE);
		assertTrue("No overlap plus strand", !fewGenes.overlaps(interval7));

		SingleInterval interval8 = new SingleInterval("chr2", 130196217, 130276013, Strand.NEGATIVE);
		assertTrue("No overlap minus strand", !fewGenes.overlaps(interval8));
		
		SingleInterval interval9 = new SingleInterval("chrM", 20, 30, Strand.NEGATIVE);
		assertTrue("Overlap chrM", fewGenes.overlaps(interval9));
		
		SingleInterval interval10 = new SingleInterval("chrM", 10, 20, Strand.POSITIVE);
		assertTrue("No overlap chrM", !fewGenes.overlaps(interval10));
		
		SingleInterval interval11 = new SingleInterval("chr7", 0, 5, Strand.POSITIVE);
		assertTrue("No overlap at beginning of chr7", !fewGenes.overlaps(interval11));
		
		SingleInterval interval12 = new SingleInterval("chr19", 61342410, 61342425, Strand.POSITIVE);
		assertTrue("Overlap at end of chr19", fewGenes.overlaps(interval12));
		
	}
	
	

	@Test
	public void featureCollectionMerge2to1() {
		Map<String,Integer> mapping = new TreeMap<String,Integer>();
		mapping.put("a1",1);
		mapping.put("a2",1000);
		CoordinateSpace fcspace = new CoordinateSpace(mapping);
		FeatureCollection<BlockedAnnotation> fc = new FeatureCollection<BlockedAnnotation>(fcspace);
		
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		a1.addBlocks(block1);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block2 = new SingleInterval("a1",200,500);		
		a2.addBlocks(block2);
		
		fc.add(a1);
		fc.add(a2);
		
		FeatureCollection<BlockedAnnotation> fc_merged = fc.merge();
		assertEquals("merged fc contains one annotation",1,fc_merged.getNumAnnotations());
	}
	
	@Test
	public void featureCollectionMerge3to1() {
		Map<String,Integer> mapping = new TreeMap<String,Integer>();
		mapping.put("a1",1);
		mapping.put("a2",1000);
		CoordinateSpace fcspace = new CoordinateSpace(mapping);
		FeatureCollection<BlockedAnnotation> fc = new FeatureCollection<BlockedAnnotation>(fcspace);
		
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		a1.addBlocks(block1);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block2 = new SingleInterval("a1",200,500);		
		a2.addBlocks(block2);
		
		BlockedAnnotation a3 = new BlockedAnnotation();
		SingleInterval block3 = new SingleInterval("a1",350,600);		
		a3.addBlocks(block3);
		
		fc.add(a1);
		fc.add(a2);
		fc.add(a3);
		
		FeatureCollection<BlockedAnnotation> fc_merged = fc.merge();
		assertEquals("merged fc contains one annotation",1,fc_merged.getNumAnnotations());
	}

	@Test
	public void featureCollectionMerge4to2() {
		Map<String,Integer> mapping = new TreeMap<String,Integer>();
		mapping.put("a1",1);
		mapping.put("a2",1000);
		CoordinateSpace fcspace = new CoordinateSpace(mapping);
		FeatureCollection<BlockedAnnotation> fc = new FeatureCollection<BlockedAnnotation>(fcspace);
		
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		a1.addBlocks(block1);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block2 = new SingleInterval("a1",200,500);		
		a2.addBlocks(block2);
		
		BlockedAnnotation a3 = new BlockedAnnotation();
		SingleInterval block3 = new SingleInterval("a1",600,700);		
		a3.addBlocks(block3);
		
		BlockedAnnotation a4 = new BlockedAnnotation();
		SingleInterval block4 = new SingleInterval("a1",650,750);		
		a4.addBlocks(block4);
		
		fc.add(a1);
		fc.add(a2);
		fc.add(a3);
		fc.add(a4);
		
		FeatureCollection<BlockedAnnotation> fc_merged = fc.merge();
		assertEquals("merged fc contains two annotations",2,fc_merged.getNumAnnotations());
	}

	
}

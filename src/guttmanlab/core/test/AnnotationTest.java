package guttmanlab.core.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.Iterator;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;

import org.junit.Before;
import org.junit.Test;

public class AnnotationTest {
		private SingleInterval block1;
		private SingleInterval block2;
		private SingleInterval block3;
		
	@Before
	public void setUp() {
		block1 = new SingleInterval("a1", 100, 300);
		block2 = new SingleInterval("a1", 350, 500);
		block3 = new SingleInterval("a1", 600, 700);
	}

	@Test
	public void testStrand() {
		Annotation pos1 = new BlockedAnnotation();
		Annotation pos2 = new BlockedAnnotation();
		Annotation neg = new BlockedAnnotation();
		Annotation both = new BlockedAnnotation();
		Annotation unknown = new BlockedAnnotation();
		
		pos1.setOrientation(Strand.POSITIVE);
		pos2.setOrientation(Strand.POSITIVE);
		neg.setOrientation(Strand.NEGATIVE);
		both.setOrientation(Strand.BOTH);
		
		Strand posStrand1 = pos1.getOrientation();
		Strand posStrand2 = pos2.getOrientation();
		Strand negStrand = neg.getOrientation();
		Strand unknownStrand = unknown.getOrientation();
		
		assertEquals("Two positive strands should be equal.", posStrand1, posStrand2);
		assertNotEquals("Two different orientations should be not equal.", posStrand2, negStrand);
		assertEquals("An annotation with an unset orientation should return Strand.UNKNOWN", unknownStrand, Strand.UNKNOWN);
	}
	
	@Test
	public void testSetOrientation() {
		BlockedAnnotation blocked = new BlockedAnnotation();
		blocked.addBlocks(block1);
		blocked.addBlocks(block2);
		blocked.addBlocks(block3);
		
		blocked.setOrientation(Strand.POSITIVE);
		Iterator<SingleInterval> iter = blocked.getBlocks();
		assertEquals("First block orientation not equal to blocked orientation", iter.next().getOrientation(), blocked.getOrientation());
		assertEquals("Middle block orientation not equal to blocked orientation", iter.next().getOrientation(), blocked.getOrientation());
		assertEquals("Last block orientation not equal to blocked orientation", iter.next().getOrientation(), blocked.getOrientation());
	}
	
	@Test 
	public void testTrimAnnotationPos() {
		BlockedAnnotation blocked = new BlockedAnnotation();
		blocked.addBlocks(block1);
		blocked.addBlocks(block2);
		blocked.addBlocks(block3);
		
		blocked.setOrientation(Strand.POSITIVE);
		Annotation blocked1 = blocked.trim(0, 800);    // no trimming
		Annotation blocked2 = blocked.trim(150, 650);  // trim first and last exon
		Annotation blocked3 = blocked.trim(330, 550);  // between exons
		
		Iterator<SingleInterval> iter = blocked1.getBlocks();
		SingleInterval trimmedBlock1 = new SingleInterval(iter.next(), block1.getName());
		SingleInterval trimmedBlock2 = new SingleInterval(iter.next(), block2.getName());
		SingleInterval trimmedBlock3 = new SingleInterval(iter.next(), block3.getName());
		
		assertEquals("Trimming should not have removed any blocks.", blocked1.getNumberOfBlocks(), 3);
		assertEquals("Trimming should not affect first block.", block1, trimmedBlock1);
		assertEquals("Trimming should not affect middle block.", block2, trimmedBlock2);
		assertEquals("Trimming should not affect last block.", block3, trimmedBlock3);
		
		iter = blocked2.getBlocks();
		trimmedBlock1 = new SingleInterval(iter.next(), block1.getName());
		trimmedBlock2 = new SingleInterval(iter.next(), block2.getName());
		trimmedBlock3 = new SingleInterval(iter.next(), block3.getName());
		
		assertEquals("Trimming should not have removed any blocks.", blocked2.getNumberOfBlocks(), 3);
		assertEquals("First block is trimmed.", new SingleInterval("a1", 150, 300, Strand.POSITIVE), trimmedBlock1);
		assertEquals("Trimming should not affect middle block.", block2, trimmedBlock2);
		assertEquals("Last block is trimmed.", new SingleInterval("a1", 600, 650, Strand.POSITIVE), trimmedBlock3);
		
		iter = blocked3.getBlocks();
		trimmedBlock2 = new SingleInterval(iter.next(), block2.getName());
		
		assertEquals("Trimming should have removed the first and last block.", blocked3.getNumberOfBlocks(), 1);
		assertEquals("Trimming should not affect middle block.", block2, trimmedBlock2);
	}
	
	@Test 
	public void testTrimAnnotationNeg() {
		BlockedAnnotation blocked = new BlockedAnnotation();
		blocked.addBlocks(block1);
		blocked.addBlocks(block2);
		blocked.addBlocks(block3);
		
		blocked.setOrientation(Strand.NEGATIVE);
		Annotation blocked1 = blocked.trim(0, 800);    // no trimming
		Annotation blocked2 = blocked.trim(150, 650);  // trim first and last exon
		Annotation blocked3 = blocked.trim(330, 550);  // between exons
		
		Iterator<SingleInterval> iter = blocked1.getBlocks();
		SingleInterval trimmedBlock1 = new SingleInterval(iter.next(), block1.getName());
		SingleInterval trimmedBlock2 = new SingleInterval(iter.next(), block2.getName());
		SingleInterval trimmedBlock3 = new SingleInterval(iter.next(), block3.getName());
		
		assertEquals("Trimming should not have removed any blocks.", blocked1.getNumberOfBlocks(), 3);
		assertEquals("Trimming should not affect first block.", block1, trimmedBlock1);
		assertEquals("Trimming should not affect middle block.", block2, trimmedBlock2);
		assertEquals("Trimming should not affect last block.", block3, trimmedBlock3);
		
		iter = blocked2.getBlocks();
		trimmedBlock1 = new SingleInterval(iter.next(), block1.getName());
		trimmedBlock2 = new SingleInterval(iter.next(), block2.getName());
		trimmedBlock3 = new SingleInterval(iter.next(), block3.getName());
		
		assertEquals("Trimming should not have removed any blocks.", blocked2.getNumberOfBlocks(), 3);
		assertEquals("First block is trimmed.", new SingleInterval("a1", 150, 300, Strand.NEGATIVE), trimmedBlock1);
		assertEquals("Trimming should not affect middle block.", block2, trimmedBlock2);
		assertEquals("Last block is trimmed.", new SingleInterval("a1", 600, 650, Strand.NEGATIVE), trimmedBlock3);
		
		iter = blocked3.getBlocks();
		trimmedBlock2 = new SingleInterval(iter.next(), block2.getName());
		
		assertEquals("Trimming should have removed the first and last block.", blocked3.getNumberOfBlocks(), 1);
		assertEquals("Trimming should not affect middle block.", block2, trimmedBlock2);
	}
	
	@Test 
	public void testHashCodeEquals() {
		BlockedAnnotation blocked1 = new BlockedAnnotation();
		BlockedAnnotation blocked2 = new BlockedAnnotation();
		BlockedAnnotation blocked3 = new BlockedAnnotation();
	
		blocked1.addBlocks(block1);
		blocked1.addBlocks(block2);
		blocked1.addBlocks(block3);
		
		blocked2.addBlocks(block1);
		blocked2.addBlocks(block2);
		blocked2.addBlocks(block3);
		
		blocked3.addBlocks(block3);
				
		assertNotEquals("blocked1 and blocked3 should not be equal.", blocked1, blocked3);
		assertEquals("blocked1 and blocked2 should be equal.", blocked1, blocked2);
		assertEquals("blocked3 and block3 should be equal.", block3, blocked3);
		assertEquals("blocked3 and block3 should have same hash code.", blocked3.hashCode(), block3.hashCode());
		assertNotEquals("blocked2 and block3 should have different hash codes.", blocked3.hashCode(), blocked2.hashCode());
	}
}

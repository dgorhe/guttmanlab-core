package guttmanlab.core.annotation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.Iterator;

import guttmanlab.core.annotation.Annotation.Strand;

import org.junit.Before;
import org.junit.Test;

public class TestBlockedAnnotation {
	
	private BlockedAnnotation blocked1;
	private BlockedAnnotation blocked2;
	private SingleInterval block1;
	private SingleInterval block2;
	private SingleInterval block3;

	
	@Before
	public void setUp() {
		blocked1 = new BlockedAnnotation("chr1");
		blocked1.addBlocks(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		blocked1.addBlocks(new SingleInterval("chr1", 400, 500, Strand.POSITIVE));
		blocked1.addBlocks(new SingleInterval("chr1", 550, 600, Strand.POSITIVE));

		blocked2 = new BlockedAnnotation("chr1");
		blocked2.addBlocks(new SingleInterval("chr1", 0, 150, Strand.POSITIVE));
		blocked2.addBlocks(new SingleInterval("chr1", 300, 600, Strand.POSITIVE));
		
		block1 = new SingleInterval("a1", 100, 300);
		block2 = new SingleInterval("a1", 350, 500);
		block3 = new SingleInterval("a1", 600, 700);

	}
	
	//@Test //appears orientation is not used by overlaps() 
	public void annotationsMustHaveCompatibleOrientations() {
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		a1.addBlocks(block1);
		a1.setOrientation(Strand.POSITIVE);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block2 = new SingleInterval("a1",200,500);		
		a2.addBlocks(block2);
		a2.setOrientation(Strand.NEGATIVE);

		Annotation a3 = a1.merge(a2);
		assertEquals("merged should have 0 blocks.",0,a3.getNumberOfBlocks());
		
		a2.setOrientation(Strand.BOTH);
		a3 = a1.merge(a2);
		assertEquals("merged should have 1 block.",1,a3.getNumberOfBlocks());
	}

	
	@Test
	public void mergeTwoSingleBlockAnnotations() {
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		a1.addBlocks(block1);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block2 = new SingleInterval("a1",200,500);		
		a2.addBlocks(block2);
		
		Annotation a3 = a1.merge(a2);
		assertEquals("merged should have 1 block.",1,a3.getNumberOfBlocks());
		assertEquals("merged start = 100",100,a3.getReferenceStartPosition());
		assertEquals("merged end = 500",500,a3.getReferenceEndPosition());
	}
	
	@Test
	public void mergeTwoBlockedAnnotations() {
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		SingleInterval block2 = new SingleInterval("a1",350,500);
		a1.addBlocks(block1);
		a1.addBlocks(block2);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block3 = new SingleInterval("a1",50,200);
		SingleInterval block4 = new SingleInterval("a1",325,400);
		a2.addBlocks(block3);
		a2.addBlocks(block4);
		
		Annotation a3 = a1.merge(a2);
		assertEquals("merged should have 2 blocks.",2,a3.getNumberOfBlocks());
		Iterator<SingleInterval> iter = a3.getBlocks();
		Annotation a31 = iter.next();
		Annotation a32 = iter.next();
		
		assertEquals("block1 start = 50",50,a31.getReferenceStartPosition());
		assertEquals("block1 end = 300",300,a31.getReferenceEndPosition());
		assertEquals("block2 start = 325",325,a32.getReferenceStartPosition());
		assertEquals("block2 end = 500",500,a32.getReferenceEndPosition());
	}
	
	@Test
	public void collapseToSingleInterval() {
		BlockedAnnotation a1 = new BlockedAnnotation();
		SingleInterval block1 = new SingleInterval("a1",100,300);
		SingleInterval block2 = new SingleInterval("a1",350,500);
		a1.addBlocks(block1);
		a1.addBlocks(block2);
		
		BlockedAnnotation a2 = new BlockedAnnotation();
		SingleInterval block3 = new SingleInterval("a1",50,200);
		SingleInterval block4 = new SingleInterval("a1",275,400);
		a1.addBlocks(block3);
		a2.addBlocks(block4);
		
		Annotation a3 = a1.merge(a2);
		assertEquals("merged should have 1 blocks.",1,a3.getNumberOfBlocks());
		assertEquals("block1 start = 50",50,a3.getReferenceStartPosition());
		assertEquals("block2 end = 500",500,a3.getReferenceEndPosition());
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
	public void testFlatten() {
		int[] flat = blocked1.flatten();
		int[] cmp = new int[] {100, 200, 400, 500, 550, 600};
		assertEquals("Flatten failed: check length", flat.length, cmp.length);
		for (int i = 0; i < flat.length; i++) {
			assertEquals("Flatten failed.", flat[i], cmp[i]);
		}
	}

	@Test
	public void testMergeTwoBlockedAnnotations() {
		Annotation intersection = blocked1.intersect(blocked2);
		BlockedAnnotation cmpIntersection = new BlockedAnnotation("chr1");
		cmpIntersection.addBlocks(new SingleInterval("chr1", 100, 150, Strand.POSITIVE));
		cmpIntersection.addBlocks(new SingleInterval("chr1", 400, 500, Strand.POSITIVE));
		cmpIntersection.addBlocks(new SingleInterval("chr1", 550, 600, Strand.POSITIVE));
		assertEquals("Intersection between two blocked annotations failed.", intersection, cmpIntersection);
		assertEquals("Intersection is not reflexive.", intersection, blocked2.intersect(blocked1));
		
		Annotation union = blocked1.union(blocked2);
		BlockedAnnotation cmpUnion = new BlockedAnnotation("chr1");
		cmpUnion.addBlocks(new SingleInterval("chr1", 0, 200, Strand.POSITIVE));
		cmpUnion.addBlocks(new SingleInterval("chr1", 300, 600, Strand.POSITIVE));
		assertEquals("Union between two blocked annotations failed.", union, cmpUnion);
		assertEquals("Union is not reflexive.", union, blocked2.union(blocked1));
		
		Annotation diff12 = blocked1.minus(blocked2);
		assertEquals("Difference between two blocked annotations failed (1).", diff12, new SingleInterval("chr1", 150, 200, Strand.POSITIVE, diff12.getName()));
		
		Annotation diff21 = blocked2.minus(blocked1);
		BlockedAnnotation cmpDiff21 = new BlockedAnnotation("chr1");
		cmpDiff21.addBlocks(new SingleInterval("chr1", 0, 100, Strand.POSITIVE));
		cmpDiff21.addBlocks(new SingleInterval("chr1", 300, 400, Strand.POSITIVE));
		cmpDiff21.addBlocks(new SingleInterval("chr1", 500, 550, Strand.POSITIVE));
		assertEquals("Difference between two blocked annotations failed (2).", diff21, cmpDiff21);
		
		Annotation xor = blocked2.xor(blocked1);
		BlockedAnnotation cmpXor = new BlockedAnnotation("chr1");
		cmpXor.addBlocks(new SingleInterval("chr1", 0, 100, Strand.POSITIVE));
		cmpXor.addBlocks(new SingleInterval("chr1", 150, 200, Strand.POSITIVE));
		cmpXor.addBlocks(new SingleInterval("chr1", 300, 400, Strand.POSITIVE));
		cmpXor.addBlocks(new SingleInterval("chr1", 500, 550, Strand.POSITIVE));
		assertEquals("Xor between two blocked annotations failed.", xor, cmpXor);
		assertEquals("Xor is not reflexive.", xor, blocked2.xor(blocked1));
	}

	
}

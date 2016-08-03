package guttmanlab.core.annotation;

import static org.junit.Assert.*;

import java.util.Iterator;

import guttmanlab.core.annotation.Annotation.Strand;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class TestBlockedAnnotation {
	
	@Rule
	public ExpectedException thrown = ExpectedException.none();

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
	
	private void checkHashCodeEquals(BlockedAnnotation b1, BlockedAnnotation b2, boolean shouldBeEqual) {
		boolean eq = b1.equals(b2);
		String s = (shouldBeEqual ? "Should be equal:" : "Should not be equal:")  + "\n" + b1.toBED() + "\n" + b2.toBED();
		assertTrue(s, shouldBeEqual == eq);
		int h1 = b1.hashCode();
		int h2 = b2.hashCode();
		if(eq) {
			assertEquals(h1, h2);
		} else {
			assert(h1 != h2);
		}
	}
	
	@Test
	public void hashCodeEquals1() {
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
	public void hashCodeEquals2() {
		
		BlockedAnnotation b1 = new BlockedAnnotation("b1");
		BlockedAnnotation b2 = new BlockedAnnotation("b2");
		BlockedAnnotation b3 = new BlockedAnnotation("b1");
		checkHashCodeEquals(b1, b2, false);
		checkHashCodeEquals(b1, b3, true);
		checkHashCodeEquals(b2, b3, false);
		b1.setOrientation(Strand.POSITIVE);
		checkHashCodeEquals(b1, b2, false);
		checkHashCodeEquals(b1, b3, false);
	}
	
	@Test
	public void hashCodeEquals3() {
		
		BlockedAnnotation b4 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.NEGATIVE), "b4");
		BlockedAnnotation b5 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.NEGATIVE), "b5");
		BlockedAnnotation b6 = new BlockedAnnotation(b4, "b6");
		BlockedAnnotation b7 = new BlockedAnnotation(b6, "b4");
		checkHashCodeEquals(b4, b5, false);
		checkHashCodeEquals(b4, b6, false);
		checkHashCodeEquals(b4, b7, true);
	}
	
	@Test
	public void hashCodeEquals4() {
		
		BlockedAnnotation b8 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.NEGATIVE));
		b8.addBlocks(new SingleInterval("chr1", 300, 400, Strand.NEGATIVE));
		BlockedAnnotation b9 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.NEGATIVE));
		b9.addBlocks(new SingleInterval("chr1", 300, 400, Strand.NEGATIVE));
		BlockedAnnotation b10 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.UNKNOWN));
		b10.addBlocks(new SingleInterval("chr1", 300, 400, Strand.UNKNOWN));
		BlockedAnnotation b11 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.UNKNOWN));
		b11.addBlocks(new SingleInterval("chr1", 300, 400, Strand.UNKNOWN, "si"));
		checkHashCodeEquals(b8, b9, true);
		checkHashCodeEquals(b9, b10, false);
		checkHashCodeEquals(b10, b11, true);
	}
	
	@Test
	public void refName() {
		assertEquals(blocked1.getReferenceName(), "chr1");
		assertEquals(block3.getReferenceName(), "a1");
		BlockedAnnotation empty = new BlockedAnnotation();
		assertNull(empty.getReferenceName());
		BlockedAnnotation nameOnly = new BlockedAnnotation("name");
		assertNull(nameOnly.getReferenceName());
		nameOnly.addBlocks(new SingleInterval("chr1", 100, 200));
		assertEquals(nameOnly.getReferenceName(), "chr1");
	}
	
	@Test
	public void numBlocks() {
		assertEquals(3, blocked1.getNumberOfBlocks());
		assertEquals(2, blocked2.getNumberOfBlocks());
		BlockedAnnotation b = new BlockedAnnotation();
		assertEquals(0, b.getNumberOfBlocks());
		b.addBlocks(new SingleInterval("chr1", 100, 200));
		assertEquals(1, b.getNumberOfBlocks());
		b.addBlocks(new SingleInterval("chr1", 150, 250));
		assertEquals(1, b.getNumberOfBlocks());
		b.addBlocks(new SingleInterval("chr1", 300, 400));
		assertEquals(2, b.getNumberOfBlocks());
		Annotation b2 = b.intersect(new SingleInterval("chr1", 50, 275));
		assertEquals(1, b2.getNumberOfBlocks());
	}
	
	@Test
	public void endpoints() {
		assertEquals(100, blocked1.getReferenceStartPosition());
		assertEquals(600, blocked1.getReferenceEndPosition());
		BlockedAnnotation b = new BlockedAnnotation();
		b.setOrientation(Strand.NEGATIVE);
		b.addBlocks(new SingleInterval("chr1", 100, 200));
		assertEquals(100, b.getReferenceStartPosition());
		b.addBlocks(new SingleInterval("chr1", 150, 250));
		assertEquals(250, b.getReferenceEndPosition());
		b.addBlocks(new SingleInterval("chr1", 300, 400));
		Annotation b2 = b.intersect(new SingleInterval("chr1", 50, 275));
		assertEquals(100, b2.getReferenceStartPosition());
		assertEquals(250, b2.getReferenceEndPosition());
		Annotation b3 = b.intersect(new SingleInterval("chr1", 125, 350));
		assertEquals(125, b3.getReferenceStartPosition());
		assertEquals(350, b3.getReferenceEndPosition());
	}
	
	@Test
	public void emptyEndpoints() {
		thrown.expect(NullPointerException.class);
		BlockedAnnotation b = new BlockedAnnotation();
		@SuppressWarnings("unused")
		int s = b.getReferenceStartPosition();
	}
	
	@Test
	public void illegalAddBlocks1() {
		thrown.expect(IllegalArgumentException.class);
		BlockedAnnotation b = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		b.addBlocks(new SingleInterval("chr1", 300, 400, Strand.NEGATIVE));
	}
	
	@Test
	public void illegalAddBlocks2() {
		thrown.expect(IllegalArgumentException.class);
		BlockedAnnotation b = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		b.addBlocks(new SingleInterval("chr1", 300, 400));
	}
	
	@Test
	public void illegalAddBlocks3() {
		thrown.expect(IllegalArgumentException.class);
		BlockedAnnotation b = new BlockedAnnotation(new SingleInterval("chr1", 100, 200));
		b.addBlocks(new SingleInterval("chr1", 300, 400, Strand.POSITIVE));
	}
	
	@Test
	public void illegalAddBlocks4() {
		thrown.expect(IllegalArgumentException.class);
		BlockedAnnotation b = new BlockedAnnotation(new SingleInterval("chr1", 100, 200));
		b.addBlocks(new SingleInterval("chr2", 300, 400));
	}
	
	@Test
	public void addBlocks() {
		BlockedAnnotation b = new BlockedAnnotation();
		assertTrue(b.size() == 0);
		b.addBlocks(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		assertTrue(1 == b.getNumberOfBlocks());
		assertTrue(b.size() == 100);
		b.addBlocks(new BlockedAnnotation(new SingleInterval("chr1", 150, 250, Strand.POSITIVE)));
		assertEquals(Strand.POSITIVE, b.getOrientation());
		assertTrue(1 == b.getNumberOfBlocks());
		assertTrue(b.size() == 150);
		b.addBlocks(new BlockedAnnotation(new SingleInterval("chr1", 350, 450, Strand.POSITIVE)));
		assertTrue(2 == b.getNumberOfBlocks());
		assertTrue(b.size() == 250);
		BlockedAnnotation b2 = new BlockedAnnotation(new SingleInterval("chr1", 400, 500, Strand.POSITIVE));
		b2.addBlocks(new SingleInterval("chr1", 600, 700, Strand.POSITIVE));
		b.addBlocks(b2);
		assertTrue(3 == b.getNumberOfBlocks());
		assertTrue(b.size() == 400);
	}
	
	@Test
	public void getOrientation() {
		BlockedAnnotation b = new BlockedAnnotation();
		assertEquals(Strand.UNKNOWN, b.getOrientation());
		b.setOrientation(Strand.NEGATIVE);
		assertEquals(Strand.NEGATIVE, b.getOrientation());
		BlockedAnnotation b2 = new BlockedAnnotation();
		b2.addBlocks(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		assertEquals(Strand.POSITIVE, b2.getOrientation());
	}
	
	@Test
	public void emptyConstructor() {
		BlockedAnnotation b = new BlockedAnnotation();
		assertNull(b.getName());
		assertNull(b.getReferenceName());
		assertEquals(Strand.UNKNOWN, b.getOrientation());
		assertTrue(0 == b.getNumberOfBlocks());
		assertEquals(b, new BlockedAnnotation());
		assertTrue(0 == b.size());
	}
	
	@Test
	public void emptyConstructorEndpoints1() {
		thrown.expect(IllegalStateException.class);
		BlockedAnnotation b = new BlockedAnnotation();
		@SuppressWarnings("unused")
		int s = b.getReferenceStartPosition();
	}
	
	@Test
	public void emptyConstructorEndpoints2() {
		thrown.expect(IllegalStateException.class);
		BlockedAnnotation b = new BlockedAnnotation();
		@SuppressWarnings("unused")
		int e = b.getReferenceEndPosition();
	}
	
	@Test
	public void nameOnlyConstructor() {
		BlockedAnnotation b = new BlockedAnnotation("name");
		assertEquals("name", b.getName());
		assertNull(b.getReferenceName());
		assertEquals(Strand.UNKNOWN, b.getOrientation());
		assertTrue(0 == b.getNumberOfBlocks());
		assertEquals(b, new BlockedAnnotation());
		assertTrue(0 == b.size());
	}
	
	@Test
	public void nameOnlyConstructorEndpoints1() {
		thrown.expect(IllegalStateException.class);
		BlockedAnnotation b = new BlockedAnnotation("name");
		@SuppressWarnings("unused")
		int s = b.getReferenceStartPosition();
	}
	
	@Test
	public void nameOnlyConstructorEndpoints2() {
		thrown.expect(IllegalStateException.class);
		BlockedAnnotation b = new BlockedAnnotation("name");
		@SuppressWarnings("unused")
		int e = b.getReferenceEndPosition();
	}
	
	@Test
	public void name() {
		BlockedAnnotation b = new BlockedAnnotation();
		assertNull(b.getName());
		b.addBlocks(new SingleInterval("chr1", 5, 10));
		assertNull(b.getName());
		BlockedAnnotation b2 = new BlockedAnnotation("b2");
		assertEquals("b2", b2.getName());
		b2.addBlocks(new SingleInterval("chr1", 20, 30));
		BlockedAnnotation b3 = new BlockedAnnotation(b2);
		assertEquals("b2", b3.getName());
		BlockedAnnotation b4 = new BlockedAnnotation(b3, "b4");
		assertEquals("b4", b4.getName());
	}

	@Test
	public void nameIntersect() {
		BlockedAnnotation b2 = new BlockedAnnotation("b2");
		assertEquals("b2", b2.getName());
		b2.addBlocks(new SingleInterval("chr1", 20, 30));
		Annotation b5 = b2.intersect(new SingleInterval("chr1", 7, 25));
		assertNull("Name of intersection should be null but is: " + b5.getName() + ".", b5.getName());
	}

	@Test
	public void nameTrim() {
		BlockedAnnotation b2 = new BlockedAnnotation("b2");
		assertEquals("b2", b2.getName());
		b2.addBlocks(new SingleInterval("chr1", 20, 30));
		Annotation b6 = b2.trim(22, 25);
		assertNull("Name of trimmed annotation should be null but is: " + b6.getName() + ".", b6.getName());
	}

	@Test //appears orientation is not used by overlaps() 
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
	public void getBlocksSeveral() {
		BlockedAnnotation b = new BlockedAnnotation();
		b.addBlocks(new SingleInterval("chr1", 100, 200));
		b.addBlocks(new SingleInterval("chr1", 500, 600));
		b.addBlocks(new SingleInterval("chr1", 300, 400));
		b.addBlocks(new SingleInterval("chr1", 450, 550));
		Iterator<SingleInterval> iter = b.getBlocks();
		assertEquals(new SingleInterval("chr1", 100, 200), iter.next());
		assertEquals(new SingleInterval("chr1", 300, 400), iter.next());
		assertEquals(new SingleInterval("chr1", 450, 600), iter.next());
		assertTrue(!iter.hasNext());
	}
	
	@Test
	public void getBlocksName() {
		BlockedAnnotation b1 = new BlockedAnnotation();
		Iterator<SingleInterval> i1 = b1.getBlocks();
		assertTrue(!i1.hasNext());
		SingleInterval s1 = new SingleInterval("chr1", 100, 200);
		b1.addBlocks(new SingleInterval(s1, "name"));
		Iterator<SingleInterval> i2 = b1.getBlocks();
		SingleInterval in1 = i2.next();
		assertTrue("Iterator should be empty", !i2.hasNext());
		assertEquals("Block should not have a name", new SingleInterval("chr1", 100, 200), in1);
		assertNull("Block should not have a name", in1.getName());
	}
	
	@Test
	public void getBlocksStrand() {
		BlockedAnnotation b1 = new BlockedAnnotation();
		b1.addBlocks(new SingleInterval("chr1", 100, 200, Strand.BOTH));
		Iterator<SingleInterval> i1 = b1.getBlocks();
		SingleInterval in1 = i1.next();
		assertEquals(Strand.BOTH, in1.getOrientation());
		BlockedAnnotation b2 = new BlockedAnnotation();
		b2.addBlocks(new SingleInterval("chr1", 100, 200));
		Iterator<SingleInterval> i2 = b2.getBlocks();
		SingleInterval in2 = i2.next();
		assertEquals(Strand.UNKNOWN, in2.getOrientation());
	}
	
	@Test
	public void size() {
		BlockedAnnotation b1 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		BlockedAnnotation b2 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.NEGATIVE));
		BlockedAnnotation b3 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.UNKNOWN));
		BlockedAnnotation b4 = new BlockedAnnotation(new SingleInterval("chr1", 100, 200, Strand.BOTH));
		assertTrue(b1.size() == 100);
		assertTrue(b2.size() == 100);
		assertTrue(b3.size() == 100);
		assertTrue(b4.size() == 100);
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
		
		assertEquals(posStrand1, Strand.POSITIVE);
		assertEquals(negStrand, Strand.NEGATIVE);
		assertEquals(unknownStrand, Strand.UNKNOWN);
		assertEquals("Two positive strands should be equal.", posStrand1, posStrand2);
		assertNotEquals("Two different orientations should be not equal.", posStrand2, negStrand);
		assertEquals("An annotation with an unset orientation should return Strand.UNKNOWN", unknownStrand, Strand.UNKNOWN);
		
		pos1.setOrientation(Strand.BOTH);
		assertEquals(pos1.getOrientation(), Strand.BOTH);
		
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

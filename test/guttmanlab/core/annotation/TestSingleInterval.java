package guttmanlab.core.annotation;

import static org.junit.Assert.*;

import java.util.Collections;
import java.util.Iterator;

import guttmanlab.core.annotation.Annotation.Strand;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class TestSingleInterval {
	
	@Rule
	public ExpectedException thrown = ExpectedException.none();
	
	
	@Test
	public void name() {
		SingleInterval s1 = new SingleInterval("chr1", 100, 200);
		assertEquals("", s1.getName());
		SingleInterval s2 = new SingleInterval(s1, "s2");
		assertEquals("s2", s2.getName());
		SingleInterval s3 = new SingleInterval("chr1", 100, 200, Strand.POSITIVE, "s3");
		assertEquals("s3", s3.getName());
		SingleInterval s4 = new SingleInterval(s3);
		assertEquals("s3", s4.getName());
	}
	
	@Test
	public void endpointsBackwards() {
		thrown.expect(IllegalArgumentException.class);
		@SuppressWarnings("unused")
		SingleInterval s1 = new SingleInterval("chr1", 200, 100);
	}
	
	@Test
	public void sizeZero() {
		thrown.expect(IllegalArgumentException.class);
		@SuppressWarnings("unused")
		SingleInterval s1 = new SingleInterval("chr1", 100, 100);
	}
	
	@Test
	public void negativeEndpoint() {
		thrown.expect(IllegalArgumentException.class);
		@SuppressWarnings("unused")
		SingleInterval s1 = new SingleInterval("chr1", -10, 100);
	}
	
	@Test
	public void invalidStrand() {
		thrown.expect(IllegalArgumentException.class);
		@SuppressWarnings("unused")
		SingleInterval s1 = new SingleInterval("chr1", 10, 100, Strand.INVALID);
	}
	
	private static void checkSingletonIterator(SingleInterval s, Iterator<SingleInterval> iter) {
		assertTrue(iter.hasNext());
		SingleInterval next = iter.next();
		assertEquals(s, next);
		assertTrue(!iter.hasNext());
	}
	
	@Test
	public void blocks() {
		SingleInterval s1 = new SingleInterval("chr1", 100, 200);
		checkSingletonIterator(s1, s1.getBlocks());
		SingleInterval s2 = new SingleInterval(s1, "s2");
		checkSingletonIterator(s2, s2.getBlocks());
		SingleInterval s3 = new SingleInterval("chr1", 100, 200, Strand.POSITIVE, "s3");
		checkSingletonIterator(s3, s3.getBlocks());
		SingleInterval s4 = new SingleInterval(s3);
		checkSingletonIterator(s4, s4.getBlocks());
	}
	
	@Test
	public void orientation() {
		SingleInterval s1 = new SingleInterval("chr1", 100, 200);
		assertEquals(Strand.UNKNOWN, s1.getOrientation());
		SingleInterval s2 = new SingleInterval(s1, "s2");
		assertEquals(Strand.UNKNOWN, s2.getOrientation());
		s2.setOrientation(Strand.NEGATIVE);
		assertEquals(Strand.UNKNOWN, s1.getOrientation());
		assertEquals(Strand.NEGATIVE, s2.getOrientation());
		SingleInterval s3 = new SingleInterval("chr1", 100, 200, Strand.POSITIVE, "s3");
		assertEquals(Strand.POSITIVE, s3.getOrientation());
	}
	
	@Test
	public void numBlocks() {
		SingleInterval pos = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		SingleInterval neg = new SingleInterval("chr1", 100, 200, Strand.NEGATIVE);
		SingleInterval both = new SingleInterval("chr1", 100, 200, Strand.BOTH);
		SingleInterval unknown = new SingleInterval("chr1", 100, 200, Strand.UNKNOWN);
		assertEquals(1, pos.getNumberOfBlocks());
		assertEquals(1, neg.getNumberOfBlocks());
		assertEquals(1, both.getNumberOfBlocks());
		assertEquals(1, unknown.getNumberOfBlocks());
	}
	
	@Test
	public void contains() {
		
		SingleInterval pos = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		SingleInterval neg = new SingleInterval("chr1", 100, 200, Strand.NEGATIVE);
		SingleInterval both = new SingleInterval("chr1", 100, 200, Strand.BOTH);
		SingleInterval unknown = new SingleInterval("chr1", 100, 200, Strand.UNKNOWN);
		SingleInterval posSmall = new SingleInterval("chr1", 125, 175, Strand.POSITIVE);
		SingleInterval negSmall = new SingleInterval("chr1", 125, 175, Strand.NEGATIVE);
		SingleInterval bothSmall = new SingleInterval("chr1", 125, 175, Strand.BOTH);
		SingleInterval unknownSmall = new SingleInterval("chr1", 125, 175, Strand.UNKNOWN);
		
		assertTrue(!pos.contains(new SingleInterval("chr1", 150, 250, Strand.POSITIVE)));
		assertTrue(!pos.contains(new SingleInterval("chr1", 225, 250, Strand.POSITIVE)));
		
		assertTrue(pos.contains(pos));
		assertTrue(!pos.contains(neg));
		assertTrue(pos.contains(both));
		assertTrue(pos.contains(unknown));
		assertTrue(pos.contains(posSmall));
		assertTrue(!pos.contains(negSmall));
		assertTrue(pos.contains(bothSmall));
		assertTrue(pos.contains(unknownSmall));
		assertTrue(!posSmall.contains(pos));
		assertTrue(!posSmall.contains(neg));
		assertTrue(!posSmall.contains(both));
		assertTrue(!posSmall.contains(unknown));
		
		assertTrue(!neg.contains(pos));
		assertTrue(neg.contains(neg));
		assertTrue(neg.contains(both));
		assertTrue(neg.contains(unknown));
		assertTrue(!neg.contains(posSmall));
		assertTrue(neg.contains(negSmall));
		assertTrue(neg.contains(bothSmall));
		assertTrue(neg.contains(unknownSmall));
		assertTrue(!negSmall.contains(pos));
		assertTrue(!negSmall.contains(neg));
		assertTrue(!negSmall.contains(both));
		assertTrue(!negSmall.contains(unknown));
		
		assertTrue(both.contains(pos));
		assertTrue(both.contains(neg));
		assertTrue(both.contains(both));
		assertTrue(both.contains(unknown));
		assertTrue(both.contains(posSmall));
		assertTrue(both.contains(negSmall));
		assertTrue(both.contains(bothSmall));
		assertTrue(both.contains(unknownSmall));
		assertTrue(!bothSmall.contains(pos));
		assertTrue(!bothSmall.contains(neg));
		assertTrue(!bothSmall.contains(both));
		assertTrue(!bothSmall.contains(unknown));
		
		assertTrue(unknown.contains(pos));
		assertTrue(unknown.contains(neg));
		assertTrue(unknown.contains(both));
		assertTrue(unknown.contains(unknown));
		assertTrue(unknown.contains(posSmall));
		assertTrue(unknown.contains(negSmall));
		assertTrue(unknown.contains(bothSmall));
		assertTrue(unknown.contains(unknownSmall));
		assertTrue(!unknownSmall.contains(pos));
		assertTrue(!unknownSmall.contains(neg));
		assertTrue(!unknownSmall.contains(both));
		assertTrue(!unknownSmall.contains(unknown));
				
	}
	
	@Test
	public void equalsAndHashCode() {
		SingleInterval s1 = new SingleInterval("chr1", 100, 200);
		SingleInterval s2 = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		SingleInterval s3 = new SingleInterval("chr1", 100, 200, Strand.BOTH);
		SingleInterval s4 = new SingleInterval("chr1", 100, 200, Strand.UNKNOWN);
		SingleInterval s5 = new SingleInterval(s1);
		SingleInterval s6 = new SingleInterval(s2, "name");
		BlockedAnnotation b1 = new BlockedAnnotation(Collections.singleton(s1), "");
		BlockedAnnotation b2 = new BlockedAnnotation(Collections.singleton(s2), "name");
		
		assertEquals(s1, s5);
		assertTrue(!s1.equals(s2));
		assertTrue(!s1.equals(s3));
		assertEquals(s1, s4);
		assertTrue(!s2.equals(s6));
		assertEquals(s1, b1);
		assertEquals(s6, b2);
		
		assertEquals(s1.hashCode(), s5.hashCode());
		assertTrue(s1.hashCode() != s2.hashCode());
		assertTrue(s1.hashCode() != s3.hashCode());
		assertEquals(s1.hashCode(), s4.hashCode());
		assertTrue(s2.hashCode() != s6.hashCode());
		assertEquals(s1.hashCode(), b1.hashCode());
		assertEquals(s6.hashCode(), b2.hashCode());
	}
	
	@Test
	public void overlaps() {
		
		SingleInterval pos = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		SingleInterval neg = new SingleInterval("chr1", 100, 200, Strand.NEGATIVE);
		SingleInterval both = new SingleInterval("chr1", 100, 200, Strand.BOTH);
		SingleInterval unknown = new SingleInterval("chr1", 100, 200, Strand.UNKNOWN);
		SingleInterval posSmall = new SingleInterval("chr1", 125, 175, Strand.POSITIVE);
		SingleInterval negSmall = new SingleInterval("chr1", 125, 175, Strand.NEGATIVE);
		SingleInterval bothSmall = new SingleInterval("chr1", 125, 175, Strand.BOTH);
		SingleInterval unknownSmall = new SingleInterval("chr1", 125, 175, Strand.UNKNOWN);
		SingleInterval posOverlap = new SingleInterval("chr1", 150, 250, Strand.POSITIVE);
		SingleInterval negOverlap = new SingleInterval("chr1", 150, 250, Strand.NEGATIVE);
		SingleInterval bothOverlap = new SingleInterval("chr1", 150, 250, Strand.BOTH);
		SingleInterval unknownOverlap = new SingleInterval("chr1", 150, 250, Strand.UNKNOWN);
		
		assertTrue(pos.overlaps(posOverlap));
		assertTrue(!pos.overlaps(negOverlap));
		assertTrue(pos.overlaps(bothOverlap));
		assertTrue(pos.overlaps(unknownOverlap));
		
		assertTrue(!pos.overlaps(new SingleInterval("chr1", 200, 201, Strand.POSITIVE)));
		assertTrue(!pos.overlaps(new SingleInterval("chr1", 99, 100, Strand.POSITIVE)));
		
		assertTrue(pos.overlaps(pos));
		assertTrue(!pos.overlaps(neg));
		assertTrue(pos.overlaps(both));
		assertTrue(pos.overlaps(unknown));
		assertTrue(pos.overlaps(posSmall));
		assertTrue(!pos.overlaps(negSmall));
		assertTrue(pos.overlaps(bothSmall));
		assertTrue(pos.overlaps(unknownSmall));
		assertTrue(posSmall.overlaps(pos));
		assertTrue(!posSmall.overlaps(neg));
		assertTrue(posSmall.overlaps(both));
		assertTrue(posSmall.overlaps(unknown));
		
		assertTrue(!neg.overlaps(pos));
		assertTrue(neg.overlaps(neg));
		assertTrue(neg.overlaps(both));
		assertTrue(neg.overlaps(unknown));
		assertTrue(!neg.overlaps(posSmall));
		assertTrue(neg.overlaps(negSmall));
		assertTrue(neg.overlaps(bothSmall));
		assertTrue(neg.overlaps(unknownSmall));
		assertTrue(!negSmall.overlaps(pos));
		assertTrue(negSmall.overlaps(neg));
		assertTrue(negSmall.overlaps(both));
		assertTrue(negSmall.overlaps(unknown));
		
		assertTrue(both.overlaps(pos));
		assertTrue(both.overlaps(neg));
		assertTrue(both.overlaps(both));
		assertTrue(both.overlaps(unknown));
		assertTrue(both.overlaps(posSmall));
		assertTrue(both.overlaps(negSmall));
		assertTrue(both.overlaps(bothSmall));
		assertTrue(both.overlaps(unknownSmall));
		assertTrue(bothSmall.overlaps(pos));
		assertTrue(bothSmall.overlaps(neg));
		assertTrue(bothSmall.overlaps(both));
		assertTrue(bothSmall.overlaps(unknown));
		
		assertTrue(unknown.overlaps(pos));
		assertTrue(unknown.overlaps(neg));
		assertTrue(unknown.overlaps(both));
		assertTrue(unknown.overlaps(unknown));
		assertTrue(unknown.overlaps(posSmall));
		assertTrue(unknown.overlaps(negSmall));
		assertTrue(unknown.overlaps(bothSmall));
		assertTrue(unknown.overlaps(unknownSmall));
		assertTrue(unknownSmall.overlaps(pos));
		assertTrue(unknownSmall.overlaps(neg));
		assertTrue(unknownSmall.overlaps(both));
		assertTrue(unknownSmall.overlaps(unknown));

	}
	
	@Test
	public void trimRelative() {
		SingleInterval pos = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		SingleInterval neg = new SingleInterval("chr1", 100, 200, Strand.NEGATIVE);
		SingleInterval both = new SingleInterval("chr1", 100, 200, Strand.BOTH);
		SingleInterval unknown = new SingleInterval("chr1", 100, 200, Strand.UNKNOWN);
		assertEquals(pos, pos.trimRelative(0, 100));
		assertEquals(neg, neg.trimRelative(0, 100));
		assertEquals(new SingleInterval("chr1", 105, 190, Strand.POSITIVE), pos.trimRelative(5, 90));
		assertEquals(new SingleInterval("chr1", 110, 195, Strand.NEGATIVE), neg.trimRelative(5, 90));
		assertEquals(new SingleInterval("chr1", 105, 190, Strand.BOTH), both.trimRelative(5, 90));
		assertEquals(new SingleInterval("chr1", 105, 190, Strand.UNKNOWN), unknown.trimRelative(5, 90));
	}
	
	@Test
	public void trimNegative() {
		thrown.expect(IllegalArgumentException.class);
		SingleInterval s1 = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		@SuppressWarnings("unused")
		SingleInterval s2 = s1.trimRelative(-1, 10);
	}
	
	@Test
	public void trimReversed() {
		thrown.expect(IllegalArgumentException.class);
		SingleInterval s1 = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		@SuppressWarnings("unused")
		SingleInterval s2 = s1.trimRelative(6, 5);
	}
	
	@Test
	public void trimOutOfBounds() {
		thrown.expect(IllegalArgumentException.class);
		SingleInterval s1 = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		@SuppressWarnings("unused")
		SingleInterval s2 = s1.trimRelative(50, 101);
	}
	
	@Test
	public void size() {
		SingleInterval pos = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		SingleInterval neg = new SingleInterval("chr1", 100, 200, Strand.NEGATIVE);
		SingleInterval both = new SingleInterval("chr1", 100, 200, Strand.BOTH);
		SingleInterval unknown = new SingleInterval("chr1", 100, 200, Strand.UNKNOWN);
		assertEquals(100, pos.size());
		assertEquals(100, neg.size());
		assertEquals(100, both.size());
		assertEquals(100, unknown.size());
	}
	
	@Test
	public void relPos() {
		SingleInterval pos = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		SingleInterval neg = new SingleInterval("chr1", 100, 200, Strand.NEGATIVE);
		SingleInterval both = new SingleInterval("chr1", 100, 200, Strand.BOTH);
		SingleInterval unknown = new SingleInterval("chr1", 100, 200, Strand.UNKNOWN);
		
		assertEquals(0, pos.getRelativePositionFrom5PrimeOfFeature(100));
		assertEquals(99, neg.getRelativePositionFrom5PrimeOfFeature(100));
		assertEquals(0, both.getRelativePositionFrom5PrimeOfFeature(100));
		assertEquals(0, unknown.getRelativePositionFrom5PrimeOfFeature(100));
		assertEquals(99, pos.getRelativePositionFrom5PrimeOfFeature(199));
		assertEquals(0, neg.getRelativePositionFrom5PrimeOfFeature(199));
		assertEquals(99, both.getRelativePositionFrom5PrimeOfFeature(199));
		assertEquals(99, unknown.getRelativePositionFrom5PrimeOfFeature(199));

		assertEquals(1, pos.getRelativePositionFrom5PrimeOfFeature(101));
		assertEquals(98, neg.getRelativePositionFrom5PrimeOfFeature(101));
		assertEquals(1, both.getRelativePositionFrom5PrimeOfFeature(101));
		assertEquals(1, unknown.getRelativePositionFrom5PrimeOfFeature(101));
		assertEquals(98, pos.getRelativePositionFrom5PrimeOfFeature(198));
		assertEquals(1, neg.getRelativePositionFrom5PrimeOfFeature(198));
		assertEquals(98, both.getRelativePositionFrom5PrimeOfFeature(198));
		assertEquals(98, unknown.getRelativePositionFrom5PrimeOfFeature(198));
}
	
	
	@Test
	public void relPosOutOfBounds1() {
		thrown.expect(IllegalArgumentException.class);
		SingleInterval s1 = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		@SuppressWarnings("unused")
		int p = s1.getRelativePositionFrom5PrimeOfFeature(99);
	}
	
	@Test
	public void relPosOutOfBounds2() {
		thrown.expect(IllegalArgumentException.class);
		SingleInterval s1 = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
		@SuppressWarnings("unused")
		int p = s1.getRelativePositionFrom5PrimeOfFeature(200);
	}
	
	@Test
	public void relPosOutOfBounds3() {
		thrown.expect(IllegalArgumentException.class);
		SingleInterval s1 = new SingleInterval("chr1", 100, 200, Strand.NEGATIVE);
		@SuppressWarnings("unused")
		int p = s1.getRelativePositionFrom5PrimeOfFeature(99);
	}
	
	@Test
	public void relPosOutOfBounds4() {
		thrown.expect(IllegalArgumentException.class);
		SingleInterval s1 = new SingleInterval("chr1", 100, 200, Strand.NEGATIVE);
		@SuppressWarnings("unused")
		int p = s1.getRelativePositionFrom5PrimeOfFeature(200);
	}
	

	@Test
	public void mergeTwoSingleBlockAnnotations() {		
		SingleInterval single1 = new SingleInterval("chr1", 100, 300, Strand.POSITIVE);
		SingleInterval single2 = new SingleInterval("chr1", 200, 400, Strand.POSITIVE);
		Annotation intersection = single1.intersect(single2);
		Annotation union = single1.union(single2);
		Annotation diff = single1.minus(single2);
		Annotation xor = single1.xor(single2);
		BlockedAnnotation cmpXor = new BlockedAnnotation("chr1");
		cmpXor.addBlocks(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		cmpXor.addBlocks(new SingleInterval("chr1", 300, 400, Strand.POSITIVE));

		assertEquals("Intersection between two single intervals failed.", intersection, new SingleInterval("chr1", 200, 300, Strand.POSITIVE, intersection.getName()));
		assertEquals("Intersection not reflexive", intersection, single2.intersect(single1));
		assertEquals("Union between two single intervals failed.", union, new SingleInterval("chr1", 100, 400, Strand.POSITIVE, union.getName()));
		assertEquals("Union not reflexive", union, single2.union(single1));
		assertEquals("Difference between two single intervals failed.", diff, new SingleInterval("chr1", 100, 200, Strand.POSITIVE, diff.getName()));
		assertEquals("Symmetric difference between two single intervals failed.", xor, cmpXor);
		assertEquals("Symmetric difference not reflexive", xor, cmpXor);
	}
	

}

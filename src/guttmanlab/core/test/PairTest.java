package guttmanlab.core.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import guttmanlab.core.datastructures.Pair;

public class PairTest {
	
	@Test
	public void testEmptyPairConstructor() {
		assertTrue(new Pair<Integer>().isEmpty());
	}
	
	@Test
	public void testPairConstructor() {
		assertTrue(new Pair<Integer>(6, 7).isComplete());
	}

	@Test
	public void testEmptyPairIsEmpty() {
		assertTrue(Pair.of(null, null).isEmpty());
	}
	
	@Test
	public void testValue1PairIsNotEmpty() {
		assertFalse(Pair.of(1, null).isEmpty());
	}
	
	@Test
	public void testValue2PairIsNotEmpty() {
		assertFalse(Pair.of(null, 2).isEmpty());
	}
	
	@Test
	public void testCompletePairIsNotEmpty() {
		assertFalse(Pair.of(1, 2).isEmpty());
	}
	
	@Test
	public void testEmptyPairIsNotComplete() {
		assertFalse(Pair.of(null, null).isComplete());
	}
	
	@Test
	public void testValue1PairIsNotComplete() {
		assertFalse(Pair.of(1, null).isComplete());
	}
	
	@Test
	public void testValue2PairIsNotComplete() {
		assertFalse(Pair.of(null, 2).isComplete());
	}
	
	@Test
	public void testCompletePairIsComplete() {
		assertTrue(Pair.of(1, 2).isComplete());
	}
	
	@Test
	public void testEmptyPairToString() {
		assertEquals("(null, null)", Pair.of(null, null).toString());
	}
	
	@Test
	public void testValue1PairToString() {
		assertEquals("(1, null)", Pair.of(1, null).toString());
	}

	@Test
	public void testValue2PairToString() {
		assertEquals("(null, 1)", Pair.of(null, 1).toString());
	}
	
	@Test
	public void testCompletePairToString() {
		assertEquals("(2, 1)", Pair.of(2, 1).toString());
	}
	
	@Test
	public void testSameReferenceEquals() {
		Pair<Integer> a = Pair.of(1, 2);
		Pair<Integer> b = a;
		assertTrue(a == b && a.equals(b));
	}
	
	@Test
	public void testDifferentTypeNotEquals() {
		Pair<Integer> a = Pair.of(1, 2);
		Pair<Double> b = Pair.of(1.0, 2.0);
		assertNotEquals(a, b);
	}
	
	@Test
	public void testSameValuesDifferentObjectEquals() {
		Pair<Integer> a = Pair.of(1, 2);
		Pair<Integer> b = Pair.of(1, 2);
		assertTrue(a != b && a.equals(b));
	}
	
	@Test
	public void testBothValuesDifferentNotEquals() {
		Pair<Integer> a = Pair.of(1, 2);
		Pair<Integer> b = Pair.of(3, 4);
		assertNotEquals(a, b);
	}
	
	@Test
	public void testEqualsSymmetry() {
		Pair<Integer> a = Pair.of(1, 2);
		Pair<Integer> b = Pair.of(1, 2);
		assertTrue(a.equals(b) && b.equals(a));
	}
}
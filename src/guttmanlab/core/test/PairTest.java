package guttmanlab.core.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import org.junit.Test;

import guttmanlab.core.datastructures.Pair;

public class PairTest {
	
	@Test
	public void testPairSimpleEqualsBehavior() {
		Pair<Integer> a = Pair.of(1, 2);
		Pair<Integer> b = Pair.of(1, 2);
		Pair<Integer> leftNull1 = Pair.of(null, 1);
		Pair<Integer> leftNull2 = Pair.of(null, 1);
		Pair<Integer> rightNull1 = Pair.of(1, null);
		Pair<Integer> rightNull2 = Pair.of(1, null);
		Pair<Integer> allNull1 = Pair.of(null, null);
		Pair<Integer> allNull2 = Pair.of(null, null);

		assertEquals(a, b);
		assertEquals(leftNull1, leftNull2);
		assertEquals(rightNull1, rightNull2);
		assertEquals(allNull1, allNull2);
	}
	
	@Test
	public void testPairSimpleNotEqualsBehavior() {
		Pair<Integer> noNull1 = Pair.of(1, 2);
		Pair<Integer> noNull2 = Pair.of(1, 3);
		Pair<Integer> noNull3 = Pair.of(2, 2);
		Pair<Integer> leftNull1 = Pair.of(null, 1);
		Pair<Integer> leftNull2 = Pair.of(null, 2);
		Pair<Integer> rightNull1 = Pair.of(1, null);
		Pair<Integer> rightNull2 = Pair.of(2, null);
		Pair<Integer> allNull = Pair.of(null, null);

		assertNotEquals(noNull1, noNull2);
		assertNotEquals(noNull2, noNull3);
		assertNotEquals(noNull1, noNull3);
		assertNotEquals(leftNull1, leftNull2);
		assertNotEquals(rightNull1, rightNull2);
		assertNotEquals(noNull1, leftNull1);
		assertNotEquals(noNull2, rightNull2);
		assertNotEquals(allNull, leftNull1);
		assertNotEquals(allNull, rightNull1);
	}
}
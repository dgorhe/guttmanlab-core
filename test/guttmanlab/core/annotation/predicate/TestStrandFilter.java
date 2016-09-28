package guttmanlab.core.annotation.predicate;

import static org.junit.Assert.*;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SingleInterval;

import org.junit.Test;

public class TestStrandFilter {
	
	private static SingleInterval plus = new SingleInterval("chr1", 100, 200, Strand.POSITIVE);
	private static SingleInterval minus = new SingleInterval("chr1", 100, 200, Strand.NEGATIVE);
	private static SingleInterval both = new SingleInterval("chr1", 100, 200, Strand.BOTH);
	private static SingleInterval unknown = new SingleInterval("chr1", 100, 200, Strand.UNKNOWN);
	
	@Test
	public void testStrandFilterPositive() {
		StrandFilter<SingleInterval> sf = new StrandFilter<SingleInterval>(Strand.POSITIVE);
		assertTrue(sf.evaluate(plus));
		assertTrue(!sf.evaluate(minus));
		assertTrue(sf.evaluate(both));
		assertTrue(!sf.evaluate(unknown));
	}
	
	@Test
	public void testStrandFilterNegative() {
		StrandFilter<SingleInterval> sf = new StrandFilter<SingleInterval>(Strand.NEGATIVE);
		assertTrue(!sf.evaluate(plus));
		assertTrue(sf.evaluate(minus));
		assertTrue(sf.evaluate(both));
		assertTrue(!sf.evaluate(unknown));
	}
	
	@Test
	public void testStrandFilterBoth() {
		StrandFilter<SingleInterval> sf = new StrandFilter<SingleInterval>(Strand.BOTH);
		assertTrue(sf.evaluate(plus));
		assertTrue(sf.evaluate(minus));
		assertTrue(sf.evaluate(both));
		assertTrue(!sf.evaluate(unknown));
	}
	
	@Test(expected=IllegalArgumentException.class)
	public void testStrandFilterUnknown() {
		StrandFilter<SingleInterval> sf = new StrandFilter<SingleInterval>(Strand.UNKNOWN);
	}
	
	@Test(expected=IllegalArgumentException.class)
	public void testStrandFilterInvalid() {
		StrandFilter<SingleInterval> sf = new StrandFilter<SingleInterval>(Strand.INVALID);
	}
	
}

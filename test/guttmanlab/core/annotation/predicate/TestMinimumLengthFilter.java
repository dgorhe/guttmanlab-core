package guttmanlab.core.annotation.predicate;

import static org.junit.Assert.assertTrue;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;

import org.junit.Test;

public class TestMinimumLengthFilter {
	
	private static BlockedAnnotation blockedPos = new BlockedAnnotation();
	private static BlockedAnnotation blockedNeg = new BlockedAnnotation();
	private static SingleInterval interval = new SingleInterval("chr1", 100, 300, Strand.POSITIVE);
	
	static {
		blockedPos.addBlocks(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		blockedPos.addBlocks(new SingleInterval("chr1", 300, 400, Strand.POSITIVE));
		blockedNeg.addBlocks(new SingleInterval("chr1", 100, 200, Strand.NEGATIVE));
		blockedNeg.addBlocks(new SingleInterval("chr1", 300, 400, Strand.NEGATIVE));
	}
	
	@Test
	public void testMinimumLengthFilterBlocked199() {
		MinimumLengthFilter<BlockedAnnotation> f199 = new MinimumLengthFilter<BlockedAnnotation>(199);
		assertTrue("Min length filter 199 on blocked annotation of length 200 should evaluate to true", f199.evaluate(blockedPos));
		assertTrue("Min length filter 199 on blocked annotation of length 200 should evaluate to true", f199.evaluate(blockedNeg));
	}

	@Test
	public void testMinimumLengthFilterBlocked200() {
		MinimumLengthFilter<BlockedAnnotation> f200 = new MinimumLengthFilter<BlockedAnnotation>(200);
		assertTrue("Min length filter 200 on blocked annotation of length 200 should evaluate to true", f200.evaluate(blockedPos));
		assertTrue("Min length filter 200 on blocked annotation of length 200 should evaluate to true", f200.evaluate(blockedNeg));
	}

	@Test
	public void testMinimumLengthFilterBlocked201() {
		MinimumLengthFilter<BlockedAnnotation> f201 = new MinimumLengthFilter<BlockedAnnotation>(201);
		assertTrue("Min length filter 201 on blocked annotation of length 200 should evaluate to false", !f201.evaluate(blockedPos));
		assertTrue("Min length filter 201 on blocked annotation of length 200 should evaluate to false", !f201.evaluate(blockedNeg));
	}

	@Test
	public void testMinimumLengthFilterBlocked250() {
		MinimumLengthFilter<BlockedAnnotation> f250 = new MinimumLengthFilter<BlockedAnnotation>(250);
		assertTrue("Min length filter 250 on blocked annotation of length 200 should evaluate to false", !f250.evaluate(blockedPos));
		assertTrue("Min length filter 250 on blocked annotation of length 200 should evaluate to false", !f250.evaluate(blockedNeg));
	}

	@Test
	public void testMinimumLengthFilterSingleInterval() {
		MinimumLengthFilter<SingleInterval> f199 = new MinimumLengthFilter<SingleInterval>(199);
		MinimumLengthFilter<SingleInterval> f200 = new MinimumLengthFilter<SingleInterval>(200);
		MinimumLengthFilter<SingleInterval> f201 = new MinimumLengthFilter<SingleInterval>(201);
		assertTrue("Min length filter 199 on single interval of length 200 should evaluate to true", f199.evaluate(interval));
		assertTrue("Min length filter 200 on single interval of length 200 should evaluate to true", f200.evaluate(interval));
		assertTrue("Min length filter 201 on single interval of length 200 should evaluate to false", !f201.evaluate(interval));
	}

}

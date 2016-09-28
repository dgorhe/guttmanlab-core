package guttmanlab.core.annotation.predicate;

import static org.junit.Assert.*;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;

import org.junit.Test;

public class TestMaximumLengthFilter {
	
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
	public void testMaximumLengthFilterBlocked199() {
		MaximumLengthFilter<BlockedAnnotation> f199 = new MaximumLengthFilter<BlockedAnnotation>(199);
		assertTrue("Max length filter 199 on blocked annotation of length 200 should evaluate to false", !f199.evaluate(blockedPos));
		assertTrue("Max length filter 199 on blocked annotation of length 200 should evaluate to false", !f199.evaluate(blockedNeg));
	}

	@Test
	public void testMaximumLengthFilterBlocked200() {
		MaximumLengthFilter<BlockedAnnotation> f200 = new MaximumLengthFilter<BlockedAnnotation>(200);
		assertTrue("Max length filter 200 on blocked annotation of length 200 should evaluate to true", f200.evaluate(blockedPos));
		assertTrue("Max length filter 200 on blocked annotation of length 200 should evaluate to true", f200.evaluate(blockedNeg));
	}

	@Test
	public void testMaximumLengthFilterBlocked201() {
		MaximumLengthFilter<BlockedAnnotation> f201 = new MaximumLengthFilter<BlockedAnnotation>(201);
		assertTrue("Max length filter 201 on blocked annotation of length 200 should evaluate to true", f201.evaluate(blockedPos));
		assertTrue("Max length filter 201 on blocked annotation of length 200 should evaluate to true", f201.evaluate(blockedNeg));
	}

	@Test
	public void testMaximumLengthFilterBlocked250() {
		MaximumLengthFilter<BlockedAnnotation> f250 = new MaximumLengthFilter<BlockedAnnotation>(250);
		assertTrue("Max length filter 250 on blocked annotation of length 200 should evaluate to true", f250.evaluate(blockedPos));
		assertTrue("Max length filter 250 on blocked annotation of length 200 should evaluate to true", f250.evaluate(blockedNeg));
	}

	@Test
	public void testMaximumLengthFilterSingleInterval() {
		MaximumLengthFilter<SingleInterval> f199 = new MaximumLengthFilter<SingleInterval>(199);
		MaximumLengthFilter<SingleInterval> f200 = new MaximumLengthFilter<SingleInterval>(200);
		MaximumLengthFilter<SingleInterval> f201 = new MaximumLengthFilter<SingleInterval>(201);
		assertTrue("Max length filter 199 on single interval of length 200 should evaluate to false", !f199.evaluate(interval));
		assertTrue("Max length filter 200 on single interval of length 200 should evaluate to true", f200.evaluate(interval));
		assertTrue("Max length filter 201 on single interval of length 200 should evaluate to true", f201.evaluate(interval));
	}

}

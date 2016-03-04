package guttmanlab.core.test;
import static org.junit.Assert.assertEquals;

import org.junit.Before;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;

import org.junit.Test;

public class AnnotationMergeTest {
	
	private SingleInterval single1;
	private SingleInterval single2;
	private BlockedAnnotation blocked1;
	private BlockedAnnotation blocked2;
	
	@Before
	public void setUp() {
		single1 = new SingleInterval("chr1", 100, 300, Strand.POSITIVE);
		single2 = new SingleInterval("chr1", 200, 400, Strand.POSITIVE);

		blocked1 = new BlockedAnnotation("chr1");
		blocked1.addBlocks(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		blocked1.addBlocks(new SingleInterval("chr1", 400, 500, Strand.POSITIVE));
		blocked1.addBlocks(new SingleInterval("chr1", 550, 600, Strand.POSITIVE));

		blocked2 = new BlockedAnnotation("chr1");
		blocked2.addBlocks(new SingleInterval("chr1", 0, 150, Strand.POSITIVE));
		blocked2.addBlocks(new SingleInterval("chr1", 300, 600, Strand.POSITIVE));
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
	public void testMergeTwoSingleBlockAnnotations() {		
		Annotation intersection = single1.intersect(single2);
		Annotation union = single1.union(single2);
		Annotation diff = single1.minus(single2);
		Annotation xor = single1.xor(single2);
		BlockedAnnotation cmpXor = new BlockedAnnotation("chr1");
		cmpXor.addBlocks(new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		cmpXor.addBlocks(new SingleInterval("chr1", 300, 400, Strand.POSITIVE));

		assertEquals("Intersection between two single intervals failed.", intersection, new SingleInterval("chr1", 200, 300, Strand.POSITIVE));
		assertEquals("Intersection not reflexive", intersection, single2.intersect(single1));
		assertEquals("Union between two single intervals failed.", union, new SingleInterval("chr1", 100, 400, Strand.POSITIVE));
		assertEquals("Union not reflexive", union, single2.union(single1));
		assertEquals("Difference between two single intervals failed.", diff, new SingleInterval("chr1", 100, 200, Strand.POSITIVE));
		assertEquals("Symmetric difference between two single intervals failed.", xor, cmpXor);
		assertEquals("Symmetric difference not reflexive", xor, cmpXor);
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
		assertEquals("Difference between two blocked annotations failed (1).", diff12, new SingleInterval("chr1", 150, 200, Strand.POSITIVE));
		
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

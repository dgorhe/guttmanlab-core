package guttmanlab.core.annotation;

import static org.junit.Assert.assertEquals;
import guttmanlab.core.annotation.Annotation.Strand;

import org.junit.Before;
import org.junit.Test;

public class TestSingleInterval {
	
	private SingleInterval single1;
	private SingleInterval single2;

	@Before
	public void setUp() {
		single1 = new SingleInterval("chr1", 100, 300, Strand.POSITIVE);
		single2 = new SingleInterval("chr1", 200, 400, Strand.POSITIVE);
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

		assertEquals("Intersection between two single intervals failed.", intersection, new SingleInterval("chr1", 200, 300, Strand.POSITIVE, intersection.getName()));
		assertEquals("Intersection not reflexive", intersection, single2.intersect(single1));
		assertEquals("Union between two single intervals failed.", union, new SingleInterval("chr1", 100, 400, Strand.POSITIVE, union.getName()));
		assertEquals("Union not reflexive", union, single2.union(single1));
		assertEquals("Difference between two single intervals failed.", diff, new SingleInterval("chr1", 100, 200, Strand.POSITIVE, diff.getName()));
		assertEquals("Symmetric difference between two single intervals failed.", xor, cmpXor);
		assertEquals("Symmetric difference not reflexive", xor, cmpXor);
	}
	

}

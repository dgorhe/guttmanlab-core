package guttmanlab.core.util;

import static org.junit.Assert.*;

import org.junit.Test;

public class TestSAMFlagDecoder {
	
	@Test
	public void allFlagsTrue() {
		SAMFlagDecoder d = new SAMFlagDecoder(4095);
		assertTrue(d.readPaired());
		assertTrue(d.properPair());
		assertTrue(d.readUnmapped());
		assertTrue(d.mateUnmapped());
		assertTrue(d.readReverseStrand());
		assertTrue(d.mateReverseStrand());
		assertTrue(d.firstInPair());
		assertTrue(d.secondInPair());
		assertTrue(d.notPrimaryAlignment());
		assertTrue(d.failsQualityChecks());
		assertTrue(d.pcrOrOpticalDuplicate());
		assertTrue(d.supplementaryAlignment());
	}
	
	@Test
	public void allFlagsFalse() {
		SAMFlagDecoder d = new SAMFlagDecoder(0);
		assertTrue(!d.readPaired());
		assertTrue(!d.properPair());
		assertTrue(!d.readUnmapped());
		assertTrue(!d.mateUnmapped());
		assertTrue(!d.readReverseStrand());
		assertTrue(!d.mateReverseStrand());
		assertTrue(!d.firstInPair());
		assertTrue(!d.secondInPair());
		assertTrue(!d.notPrimaryAlignment());
		assertTrue(!d.failsQualityChecks());
		assertTrue(!d.pcrOrOpticalDuplicate());
		assertTrue(!d.supplementaryAlignment());
	}
	
	@Test
	public void alternatingFlagsTrue() {
		SAMFlagDecoder d = new SAMFlagDecoder(1365);
		assertTrue(d.readPaired());
		assertTrue(!d.properPair());
		assertTrue(d.readUnmapped());
		assertTrue(!d.mateUnmapped());
		assertTrue(d.readReverseStrand());
		assertTrue(!d.mateReverseStrand());
		assertTrue(d.firstInPair());
		assertTrue(!d.secondInPair());
		assertTrue(d.notPrimaryAlignment());
		assertTrue(!d.failsQualityChecks());
		assertTrue(d.pcrOrOpticalDuplicate());
		assertTrue(!d.supplementaryAlignment());
	}
	
	
	
}

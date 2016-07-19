package guttmanlab.core.annotation;

import static org.junit.Assert.assertEquals;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.BEDFileRecord.BEDBuilder;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class TestBEDFileRecord {
	
	private final static double DEFAULT_SCORE = 0;
	private final static Color DEFAULT_COLOR = Color.BLACK;
	
	private static final String CHR_1 = "chr1";
	private static final int POS_1 = 100;
	private static final int POS_2 = 200;
	private static final int POS_3 = 300;
	private static final int POS_4 = 400;
	private static final String NAME_1 ="name1";
	private static final double SCORE_1 = 25;
	private static final Color COLOR_1 = new Color(150, 150, 150);
	private static final String COLOR_STRING_1 = "150,150,150";
	private static final String BLOCK_STARTS_1 = "0,200,";
	private static final String BLOCK_LENS_1 = "100,100,";
	
	private static final double EPSILON = 1e-15; // for testing equality of doubles
	
	private static Annotation simplestAnnot = new SingleInterval(CHR_1, POS_1, POS_2);
	private static BlockedAnnotation complexAnnot; // Build in setUp()
	
	private static FeatureCollection<BEDFileRecord> refSeqChr1; // All the RefSeq annotations on chr1
	private static BEDFileRecord geneMinus; // Some gene on the minus strand
	private static BEDFileRecord genePlus; // Some gene on the plus strand

	
	@Rule
	public final ExpectedException parseException = ExpectedException.none();
	
	@BeforeClass
	public static void setUp() {

		// Build complexAnnot
		complexAnnot = new BlockedAnnotation(NAME_1);
		complexAnnot.addBlocks(new SingleInterval(CHR_1, POS_1, POS_2, Strand.POSITIVE));
		complexAnnot.addBlocks(new SingleInterval(CHR_1, POS_3, POS_4, Strand.POSITIVE));
		
		// Read RefSeq chr1 annotation
		File bedFile = new File("resources/RefSeq_mm9.chr1.bed");
		Map<String, FeatureCollection<BEDFileRecord>> genesByReference = null;
		try {
			genesByReference = BEDFileIO.loadFromFileByReferenceName(bedFile, CoordinateSpace.MM9);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		refSeqChr1 = genesByReference.get("chr1");
		
		// Create geneMinus and genePlus
		Iterator<BEDFileRecord> iter = refSeqChr1.iterator();
		iter.hasNext();
		geneMinus = iter.next();
		while(iter.hasNext()) {
			genePlus = iter.next();
			if(genePlus.getOrientation().equals(Strand.POSITIVE)) break;
		}

	}

	@Test
	public void testSimplestBEDBuilder() {
		BEDFileRecord simplest = (new BEDBuilder(simplestAnnot)).build();
		assertEquals(simplest.getReferenceName(), CHR_1);
		assertEquals(simplest.getReferenceStartPosition(), POS_1);
		assertEquals(simplest.thickStart(), POS_1);
		assertEquals(simplest.thickEnd(), POS_1);
		assertEquals(simplest.getReferenceEndPosition(), POS_2);		
		assertEquals(simplest.getName(), "");
		assertEquals(simplest.color(), DEFAULT_COLOR);
		assertEquals(simplest.score(), DEFAULT_SCORE, EPSILON);
		assertEquals(simplest.getNumberOfBlocks(), 1);
	}
	
	@Test
	public void testComplexBEDBuilder() {
		BEDFileRecord complex = (new BEDBuilder(complexAnnot)).score(SCORE_1)
						   									  .color(COLOR_1)
						   									  .thickStart(POS_3)
						   									  .thickEnd(POS_4)
						   									  .build();
		assertEquals(complex.getReferenceName(), CHR_1);
		assertEquals(complex.getReferenceStartPosition(), POS_1);
		assertEquals(complex.thickStart(), POS_3);
		assertEquals(complex.thickEnd(), POS_4);
		assertEquals(complex.getReferenceEndPosition(), POS_4);		
		assertEquals(complex.getName(), NAME_1);
		assertEquals(complex.color(), COLOR_1);
		assertEquals(complex.score(), SCORE_1, EPSILON);
		assertEquals(complex.getNumberOfBlocks(), 2);
	}
	
	/************************************************
	 * Test String-to-BED and BED-to-String parsing *
	 ************************************************/
	
	@Test
	public void testStringParse1() {
		String str1 = CHR_1;
		parseException.expect(IllegalArgumentException.class);
		BEDFileRecord.fromFormattedString(str1);
	}
	
	@Test
	public void testStringParse2() {
		String str2 = CHR_1 + "\t" + POS_1;
		parseException.expect(IllegalArgumentException.class);
		BEDFileRecord.fromFormattedString(str2);
	}
	
	@Test
	public void testStringParse3() {
		String str3 = CHR_1 + "\t" + POS_1 + "\t" + POS_2;
		BEDFileRecord threeField = BEDFileRecord.fromFormattedString(str3);
		assertEquals(threeField.getReferenceName(), CHR_1);
		assertEquals(threeField.getReferenceStartPosition(), POS_1);
		assertEquals(threeField.thickStart(), POS_1);
		assertEquals(threeField.thickEnd(), POS_1);
		assertEquals(threeField.getReferenceEndPosition(), POS_2);		
		assertEquals(threeField.getName(), "");
		assertEquals(threeField.getOrientation(), Strand.UNKNOWN);
		assertEquals(threeField.color(), DEFAULT_COLOR);
		assertEquals(threeField.score(), DEFAULT_SCORE, EPSILON);
		assertEquals(threeField.getNumberOfBlocks(), 1);
		assertEquals(str3, threeField.toFormattedString(3));
	}
	
	@Test
	public void testStringParse4() {
		String str4 = CHR_1 + "\t" + POS_1 + "\t" + POS_2 + "\t" + NAME_1;
		BEDFileRecord fourField = BEDFileRecord.fromFormattedString(str4);
		assertEquals(fourField.getReferenceName(), CHR_1);
		assertEquals(fourField.getReferenceStartPosition(), POS_1);
		assertEquals(fourField.thickStart(), POS_1);
		assertEquals(fourField.thickEnd(), POS_1);
		assertEquals(fourField.getReferenceEndPosition(), POS_2);		
		assertEquals(fourField.getName(), NAME_1);
		assertEquals(fourField.getOrientation(), Strand.UNKNOWN);
		assertEquals(fourField.color(), DEFAULT_COLOR);
		assertEquals(fourField.score(), DEFAULT_SCORE, EPSILON);
		assertEquals(fourField.getNumberOfBlocks(), 1);
		assertEquals(str4, fourField.toFormattedString(4));
	}
	
	@Test
	public void testStringParse5() {
		String str5 = CHR_1 + "\t" + POS_1 + "\t" + POS_2 + "\t" + NAME_1 + "\t" + SCORE_1;
		BEDFileRecord fiveField = BEDFileRecord.fromFormattedString(str5);
		assertEquals(fiveField.getReferenceName(), CHR_1);
		assertEquals(fiveField.getReferenceStartPosition(), POS_1);
		assertEquals(fiveField.thickStart(), POS_1);
		assertEquals(fiveField.thickEnd(), POS_1);
		assertEquals(fiveField.getReferenceEndPosition(), POS_2);		
		assertEquals(fiveField.getName(), NAME_1);
		assertEquals(fiveField.getOrientation(), Strand.UNKNOWN);
		assertEquals(fiveField.color(), DEFAULT_COLOR);
		assertEquals(fiveField.score(), SCORE_1, EPSILON);
		assertEquals(fiveField.getNumberOfBlocks(), 1);
		assertEquals(str5, fiveField.toFormattedString(5));
	}
	
	@Test
	public void testStringParse6() {
		String str6 = CHR_1 + "\t" + POS_1 + "\t" + POS_2 + "\t" + NAME_1 + "\t" + SCORE_1 + "\t+";
		BEDFileRecord sixField = BEDFileRecord.fromFormattedString(str6);
		assertEquals(sixField.getReferenceName(), CHR_1);
		assertEquals(sixField.getReferenceStartPosition(), POS_1);
		assertEquals(sixField.thickStart(), POS_1);
		assertEquals(sixField.thickEnd(), POS_1);
		assertEquals(sixField.getReferenceEndPosition(), POS_2);		
		assertEquals(sixField.getName(), NAME_1);
		assertEquals(sixField.getOrientation(), Strand.POSITIVE);
		assertEquals(sixField.color(), DEFAULT_COLOR);
		assertEquals(sixField.score(), SCORE_1, EPSILON);
		assertEquals(sixField.getNumberOfBlocks(), 1);
		assertEquals(str6, sixField.toFormattedString(6));
	}
	
	@Test
	public void testStringParse7() {
		String str7 = CHR_1 + "\t" + POS_1 + "\t" + POS_4 + "\t" + NAME_1 + "\t" + SCORE_1 + "\t+\t" +
				POS_3;
		parseException.expect(IllegalArgumentException.class);
		BEDFileRecord.fromFormattedString(str7);
	}
	
	@Test
	public void testStringParse8() {
		String str8 = CHR_1 + "\t" + POS_1 + "\t" + POS_4 + "\t" + NAME_1 + "\t" + SCORE_1 + "\t+\t" +
				POS_3 + "\t" + POS_4;
		BEDFileRecord eightField = BEDFileRecord.fromFormattedString(str8);
		assertEquals(eightField.getReferenceName(), CHR_1);
		assertEquals(eightField.getReferenceStartPosition(), POS_1);
		assertEquals(eightField.thickStart(), POS_3);
		assertEquals(eightField.thickEnd(), POS_4);
		assertEquals(eightField.getReferenceEndPosition(), POS_4);		
		assertEquals(eightField.getName(), NAME_1);
		assertEquals(eightField.getOrientation(), Strand.POSITIVE);
		assertEquals(eightField.color(), DEFAULT_COLOR);
		assertEquals(eightField.score(), SCORE_1, EPSILON);
		assertEquals(eightField.getNumberOfBlocks(), 1);
		assertEquals(str8, eightField.toFormattedString(8));
	}
	
	@Test
	public void testStringParse9() {
		String str9 = CHR_1 + "\t" + POS_1 + "\t" + POS_4 + "\t" + NAME_1 + "\t" + SCORE_1 + "\t+\t" +
				POS_3 + "\t" + POS_4 + "\t" + COLOR_STRING_1; 
		BEDFileRecord nineField = BEDFileRecord.fromFormattedString(str9);
		assertEquals(nineField.getReferenceName(), CHR_1);
		assertEquals(nineField.getReferenceStartPosition(), POS_1);
		assertEquals(nineField.thickStart(), POS_3);
		assertEquals(nineField.thickEnd(), POS_4);
		assertEquals(nineField.getReferenceEndPosition(), POS_4);		
		assertEquals(nineField.getName(), NAME_1);
		assertEquals(nineField.getOrientation(), Strand.POSITIVE);
		assertEquals(nineField.color(), COLOR_1);
		assertEquals(nineField.score(), SCORE_1, EPSILON);
		assertEquals(nineField.getNumberOfBlocks(), 1);
		assertEquals(str9, nineField.toFormattedString(9));
	}
	
	@Test
	public void testStringParse10() {
		String str10 = CHR_1 + "\t" + POS_1 + "\t" + POS_4 + "\t" + NAME_1 + "\t" + SCORE_1 + "\t+\t" +
				POS_3 + "\t" + POS_4 + "\t" + COLOR_STRING_1 + "\t2";
		parseException.expect(IllegalArgumentException.class);
		BEDFileRecord.fromFormattedString(str10);
	}
	
	@Test
	public void testStringParse11() {
		String str11 = CHR_1 + "\t" + POS_1 + "\t" + POS_4 + "\t" + NAME_1 + "\t" + SCORE_1 + "\t+\t" +
				POS_3 + "\t" + POS_4 + "\t" + COLOR_STRING_1 + "\t2\t" + BLOCK_LENS_1;
		parseException.expect(IllegalArgumentException.class);
		BEDFileRecord.fromFormattedString(str11);
	}
	
	@Test
	public void testStringParse12() {
		String str12 = CHR_1 + "\t" + POS_1 + "\t" + POS_4 + "\t" + NAME_1 + "\t" + SCORE_1 + "\t+\t" +
				POS_3 + "\t" + POS_4 + "\t" + COLOR_STRING_1 + "\t2\t" + BLOCK_LENS_1 + 
				"\t" + BLOCK_STARTS_1; 
		BEDFileRecord twelveField = BEDFileRecord.fromFormattedString(str12);
		assertEquals(twelveField.getReferenceName(), CHR_1);
		assertEquals(twelveField.getReferenceStartPosition(), POS_1);
		assertEquals(twelveField.thickStart(), POS_3);
		assertEquals(twelveField.thickEnd(), POS_4);
		assertEquals(twelveField.getReferenceEndPosition(), POS_4);		
		assertEquals(twelveField.getName(), NAME_1);
		assertEquals(twelveField.getOrientation(), Strand.POSITIVE);
		assertEquals(twelveField.color(), COLOR_1);
		assertEquals(twelveField.score(), SCORE_1, EPSILON);
		assertEquals(twelveField.getNumberOfBlocks(), 2);
		assertEquals(str12, twelveField.toFormattedString(12));
	}
	
	@Test
	public void testWindowsMinusStrand() {
		AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows = geneMinus.getWindows(100, 50);
		assertEquals(71, windows.getNumAnnotations());
		Iterator<DerivedAnnotation<? extends Annotation>> iter = windows.sortedIterator();
		iter.hasNext();
		assertEquals("chr1:3204562-3204662", iter.next().toUCSC());
		iter.hasNext();
		assertEquals("chr1:3204612-3204712", iter.next().toUCSC());
		iter.hasNext();
		assertEquals("chr1:3204662-3204762", iter.next().toUCSC());
		iter.hasNext();
		assertEquals("chr1:3204712-3204812", iter.next().toUCSC());
		iter.hasNext();
		assertEquals("chr1:3204762-3204862", iter.next().toUCSC());
		iter.hasNext();
		assertEquals(Strand.NEGATIVE, iter.next().getOrientation());
		DerivedAnnotation<? extends Annotation> lastWindow = null;
		while(iter.hasNext()) {
			lastWindow = iter.next();
		}
		assertEquals("chr1:3661445-3661545", lastWindow.toUCSC());
	}
	
	@Test
	public void testWindowsPlusStrand() {
		AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows = genePlus.getWindows(200, 100);
		assertEquals(23, windows.getNumAnnotations());
		Iterator<DerivedAnnotation<? extends Annotation>> iter = windows.sortedIterator();
		iter.hasNext();
		assertEquals("chr1:4797973-4820360", iter.next().toUCSC());
		iter.hasNext();
		assertEquals("chr1:4798545-4822455", iter.next().toUCSC());
		iter.hasNext();
		assertEquals("chr1:4820360-4829486", iter.next().toUCSC());
		iter.hasNext();
		assertEquals("chr1:4822455-4831053", iter.next().toUCSC());
		iter.hasNext();
		assertEquals("chr1:4829486-4831153", iter.next().toUCSC());
		iter.hasNext();
		assertEquals(Strand.POSITIVE, iter.next().getOrientation());
		DerivedAnnotation<? extends Annotation> lastWindow = null;
		while(iter.hasNext()) {
			lastWindow = iter.next();
		}
		assertEquals("chr1:4836583-4836783", lastWindow.toUCSC());
	}
	


}

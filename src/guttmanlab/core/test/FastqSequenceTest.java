package guttmanlab.core.test;

import org.junit.Before;
import org.junit.Test;

import guttmanlab.core.sequence.FastqSequence;
import guttmanlab.core.sequence.PhredEncoding;

import static org.junit.Assert.assertTrue;

import java.lang.reflect.Field;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThat;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.not;

public class FastqSequenceTest {
	private String name1 = "name1";
    private String bases1 = "ATCGTAGTAGCGTAG";
    private String bases1_rc = "CTACGCTACTACGAT";
    private String quals1Sanger = "FGDDGHGIII#####";
    private String quals1SangerRev = "#####IIIGHGDDGF";
    private byte[] quals1 = {37, 38, 35, 35, 38, 39, 38, 40, 40, 40, 2, 2, 2, 2, 2};
    private String seq1String = "@" + name1 + "\n" + bases1 + "\n+\n" + quals1Sanger;
    private String seq1FastaString = ">" + name1 + "\n" + bases1;

    private String name2 = "name2";
    private String bases2 = "GCTAGATGGCTGATA";
    private String quals2Sanger = "FGD!GHGIII##!##";
    
    private String shortQuals = "IIII";
    
    private String nameAllA = "allA";
    private String basesAllA = "aAAaaAAAAAAAAaA";
    
    private String nameAllT = "allT";
    private String basesAllT = "TtTTTtTttTTTTTT";
    
    private String nameAAndT = "aAndT";
    private String basesAAndT = "TtTTAaAttTTTTTT";
    
    private String perfectQuals = "IIIIIIIIIIIIIII";
    
    private FastqSequence seq1;
    private FastqSequence seq1_rc;
    private FastqSequence seq2_1;
    private FastqSequence seq2_2;
    private FastqSequence allA;
    private FastqSequence allT;
    private FastqSequence aAndT;
    
    @Before
    public void setUp() {
        seq1 = new FastqSequence(name1, bases1, quals1Sanger);
        seq1_rc = new FastqSequence(name1, bases1_rc, quals1SangerRev);
        seq2_1 = new FastqSequence(name2, bases2, quals2Sanger);
        seq2_2 = new FastqSequence(name2, bases2, quals2Sanger);
        allA = new FastqSequence(nameAllA, basesAllA, perfectQuals);
        allT = new FastqSequence(nameAllT, basesAllT, perfectQuals);
        aAndT = new FastqSequence(nameAAndT, basesAAndT, perfectQuals);
    }
    
    @Test(expected = IllegalArgumentException.class)
    public void nullNameTest() {
        new FastqSequence(null, bases1, quals1Sanger, PhredEncoding.SANGER);
    }
    
    @Test(expected = IllegalArgumentException.class)
    public void nullSequenceTest() {
    	new FastqSequence(name1, null, quals1Sanger, PhredEncoding.SANGER);
    }

    @Test(expected = IllegalArgumentException.class)
    public void nullQualsTest() {
    	new FastqSequence(name1, bases1, null, PhredEncoding.SANGER);
    }
    
    @Test(expected = IllegalArgumentException.class)
    public void nullEncodingTest() {
    	new FastqSequence(name1, bases1, quals1Sanger, null);
    }
    
    @Test(expected = IllegalArgumentException.class)
    public void baseQualsLengthMismatchTest() {
    	new FastqSequence(name1, bases1, shortQuals);
    }

    public void sangerToByteStringTest() throws NoSuchFieldException, SecurityException, IllegalArgumentException, IllegalAccessException {
        Field field = seq1.getClass().getDeclaredField("quality");
        field.setAccessible(true);
        assertThat(quals1, is(equalTo(field.get(seq1))));
    }
    
    @Test
    public void reverseComplementPosTest() {
        assertThat(seq1.reverseComplement(), is(equalTo(seq1_rc)));
    }

    @Test
    public void reverseComplementNegTest() {
        assertThat(seq1.reverseComplement(), is(not(equalTo(seq1))));
    }

    @Test
    public void isPolyANegTest() {
        assertFalse(seq1.isPolyA());
    }

    @Test
    public void isPolyAPosTest1() {
        assertTrue(allA.isPolyA());
    }
    
    @Test
    public void isPolyaPosTest2() {
        assertTrue(allT.isPolyA());
    }
    
    @Test
    public void isPolyaMixedTest() {
        assertFalse(aAndT.isPolyA());
    }
    
    @Test
    public void equalsSameObjectTest() {
        assertThat(seq1, is(equalTo(seq1)));
    }
    
    @Test
    public void equalsDifferentObjectTest() {
        assertThat(seq2_1, is(equalTo(seq2_2)));
    }
    
    @Test
    public void equalsNegTest() {
        assertThat(seq2_1, is(not(equalTo(seq1))));
    }
    
    @Test
    public void lengthTest() {
        assertThat(bases1.length(), is(equalTo(seq1.length())));
    }
    
    @Test
    public void toFastaTest() {
        assertThat(seq1.toFasta(), is(equalTo(seq1FastaString)));
    }
    
    @Test
    public void toFormattedStringTest() {
        assertThat(seq1.toFormattedString(), is(equalTo(seq1String)));
    }
}
package guttmanlab.core.test;

import org.junit.Before;
import org.junit.Test;

import guttmanlab.core.sequence.Sequence;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThat;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.equalTo;
import static org.hamcrest.CoreMatchers.not;

public class SequenceTest {

	private String name1 = "name1";
    private String bases1 = "ATCGTAGTAGCGTAG";
    private String bases1_rc = "CTACGCTACTACGAT";    
    private String seq1_string = ">" + name1 + "\n" + bases1;

    private String name2 = "name2";
    private String bases2 = "GCTAGATGGCTGATA";

    private String nameAllA = "allA";
    private String basesAllA = "aAAaaAAAAAAAAaA";
    
    private String nameAllT = "allT";
    private String basesAllT = "TtTTTtTttTTTTTT";
    
    private String nameAAndT = "aAndT";
    private String basesAAndT = "TtTTAaAttTTTTTT";
    
    private Sequence seq1;
    private Sequence seq1_rc;
    private Sequence seq2_1;
    private Sequence seq2_2;
    private Sequence allA;
    private Sequence allT;
    private Sequence aAndT;
    
    @Before
    public void setUp() {
        seq1 = new Sequence(name1, bases1);
        seq1_rc = new Sequence(name1, bases1_rc);
        seq2_1 = new Sequence(name2, bases2);
        seq2_2 = new Sequence(name2, bases2);
        allA = new Sequence(nameAllA, basesAllA);
        allT = new Sequence(nameAllT, basesAllT);
        aAndT = new Sequence(nameAAndT, basesAAndT);
    }
    
    @Test(expected = IllegalArgumentException.class)
    public void nullNameTest() {
        new Sequence(null, bases1);
    }
    
    @Test(expected = IllegalArgumentException.class)
    public void nullSequenceTest() {
        new Sequence(name1, null);
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
    public void toFormattedStringTest() {
        assertThat(seq1.toFasta(), is(equalTo(seq1_string)));
    }
}
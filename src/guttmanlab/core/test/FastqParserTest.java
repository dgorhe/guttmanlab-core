package guttmanlab.core.test;

import org.junit.Test;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThat;
import static org.hamcrest.CoreMatchers.is;
import static org.hamcrest.CoreMatchers.equalTo;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Collection;

import guttmanlab.core.sequence.FastqParser;
import guttmanlab.core.sequence.FastqSequence;
import guttmanlab.core.sequence.PhredEncoding;
import guttmanlab.core.sequence.Sequence;

public class FastqParserTest {
    
    private final URL completeFastq = this.getClass().getResource("/guttmanlab/core/test/CompleteFastq.fastq");
    private final static int NUM_COMPLETE_RECORDS = 3;
    
    private final URL incompleteFastq = this.getClass().getResource("/guttmanlab/core/test/IncompleteFastq.fastq");
    private final static int NUM_INCOMPLETE_RECORDS = 2;

    private final String name1 = "READ1: Lorem ipsum dolor sit amet, consectetur adipiscing elit.";
    private final String seq1 = "ATCGGATTAGGGCTATGGATAGGGCTA";
    private final String quals1 = "IIIIIIIIIIIIIIIIIIIIIIIIIII";
    
    @Test(expected = IllegalArgumentException.class)
    public void nullPathTest() throws IOException {
        try (FastqParser p = new FastqParser(null, PhredEncoding.SANGER)) {};
    }
    
    @Test(expected = IllegalArgumentException.class)
    public void nullEncodingTest() throws IOException {
        try (FastqParser p = new FastqParser((new File(completeFastq.getPath())).toPath(), null)) {};
    }
    
    @Test
    public void completeIterateTest() throws IOException {
        try (FastqParser p = new FastqParser((new File(completeFastq.getPath())).toPath(), PhredEncoding.SANGER)) {
            for (int i = 0; i < NUM_COMPLETE_RECORDS; i++) {
                assertTrue(p.hasNext());
                p.next();
            }
            assertFalse(p.hasNext());
        };
    }
    
    @Test
    public void completeIterateFirstRecordTest() throws IOException {
        try (FastqParser p = new FastqParser((new File(completeFastq.getPath())).toPath(), PhredEncoding.SANGER)) {
            assertThat(p.next(), is(equalTo(new FastqSequence(name1, seq1, quals1, PhredEncoding.SANGER))));
        };
    }
    
    @Test
    public void incompleteIterateTest() throws IOException {
        try (FastqParser p = new FastqParser((new File(incompleteFastq.getPath())).toPath(), PhredEncoding.SANGER)) {
            for (int i = 0; i < NUM_INCOMPLETE_RECORDS; i++) {
                assertTrue(p.hasNext());
                p.next();
            }
            assertFalse(p.hasNext());
        };
    }
    
    @Test
    public void incompleteIterateFirstRecordTest() throws IOException {
        try (FastqParser p = new FastqParser((new File(incompleteFastq.getPath())).toPath(), PhredEncoding.SANGER)) {
            assertThat(p.next(), is(equalTo(new FastqSequence(name1, seq1, quals1, PhredEncoding.SANGER))));
        };
    }
    
    @Test
    public void completeSlurpSizeTest() throws IOException {
        Collection<Sequence> c = FastqParser.slurp(new File(completeFastq.getPath()).toPath(), PhredEncoding.SANGER);
        assertThat(c.size(), is(equalTo(NUM_COMPLETE_RECORDS)));
    }
    
    @Test
    public void completeSlurpFirstRecordTest() throws IOException {
        Collection<Sequence> c = FastqParser.slurp(new File(completeFastq.getPath()).toPath(), PhredEncoding.SANGER);
        Sequence s = c.iterator().next(); // Depends on the underlying Collection being an ArrayList
        assertThat(s, is(equalTo(new FastqSequence(name1, seq1, quals1, PhredEncoding.SANGER))));
    }
    
    @Test
    public void incompleteSlurpSizeTest() throws IOException {
        Collection<Sequence> c = FastqParser.slurp(new File(incompleteFastq.getPath()).toPath(), PhredEncoding.SANGER);
        assertThat(c.size(), is(equalTo(NUM_INCOMPLETE_RECORDS)));
    }
    
    @Test
    public void incompleteSlurpFirstRecordTest() throws IOException {
        Collection<Sequence> c = FastqParser.slurp(new File(incompleteFastq.getPath()).toPath(), PhredEncoding.SANGER);
        Sequence s = c.iterator().next(); // Depends on the underlying Collection being an ArrayList
        
        assertThat(s, is(equalTo(new FastqSequence(name1, seq1, quals1, PhredEncoding.SANGER))));
    }
}
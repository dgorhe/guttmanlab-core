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

import guttmanlab.core.sequence.FastaParser;
import guttmanlab.core.sequence.Sequence;

public class FastaParserTest {
    
    private final URL completeFasta = this.getClass().getResource("/guttmanlab/core/test/CompleteFasta.fasta");
    private final static int NUM_COMPLETE_RECORDS = 3;
    
    private final URL incompleteFasta = this.getClass().getResource("/guttmanlab/core/test/IncompleteFasta.fasta");
    private final static int NUM_INCOMPLETE_RECORDS = 2;

    private final String name1 = "READ1: Lorem ipsum dolor sit amet, consectetur adipiscing elit.";
    private final String seq1 = "ATCGGATTAGGGCTATGGATAGGGCTA";
    
    @Test(expected = IllegalArgumentException.class)
    public void nullPathTest() throws IOException {
        try (FastaParser p = new FastaParser(null)) {};
    }
    
    @Test
    public void completeIterateTest() throws IOException {
        try (FastaParser p = new FastaParser((new File(completeFasta.getPath())).toPath())) {
            for (int i = 0; i < NUM_COMPLETE_RECORDS; i++) {
                assertTrue(p.hasNext());
                p.next();
            }
            assertFalse(p.hasNext());
        };
    }
    
    @Test
    public void completeIterateFirstRecordTest() throws IOException {
        try (FastaParser p = new FastaParser((new File(completeFasta.getPath())).toPath())) {
            assertThat(p.next(), is(equalTo(new Sequence(name1, seq1))));
        };
    }
    
    @Test
    public void incompleteIterateTest() throws IOException {
        try (FastaParser p = new FastaParser((new File(incompleteFasta.getPath())).toPath())) {
            for (int i = 0; i < NUM_INCOMPLETE_RECORDS; i++) {
                assertTrue(p.hasNext());
                p.next();
            }
            assertFalse(p.hasNext());
        };
    }
    
    @Test
    public void incompleteIterateFirstRecordTest() throws IOException {
        try (FastaParser p = new FastaParser((new File(incompleteFasta.getPath())).toPath())) {
            assertThat(p.next(), is(equalTo(new Sequence(name1, seq1))));
        };
    }
    
    @Test
    public void completeSlurpSizeTest() throws IOException {
        Collection<Sequence> c = FastaParser.slurp(new File(completeFasta.getPath()).toPath());
        assertThat(c.size(), is(equalTo(NUM_COMPLETE_RECORDS)));
    }
    
    @Test
    public void completeSlurpFirstRecordTest() throws IOException {
        Collection<Sequence> c = FastaParser.slurp(new File(completeFasta.getPath()).toPath());
        Sequence s = c.iterator().next(); // Depends on the underlying Collection being an ArrayList
        assertThat(s, is(equalTo(new Sequence(name1, seq1))));
    }
    
    @Test
    public void incompleteSlurpSizeTest() throws IOException {
        Collection<Sequence> c = FastaParser.slurp(new File(incompleteFasta.getPath()).toPath());
        assertThat(c.size(), is(equalTo(NUM_INCOMPLETE_RECORDS)));
    }
    
    @Test
    public void incompleteSlurpFirstRecordTest() throws IOException {
        Collection<Sequence> c = FastaParser.slurp(new File(incompleteFasta.getPath()).toPath());
        Sequence s = c.iterator().next(); // Depends on the underlying Collection being an ArrayList
        
        assertThat(s, is(equalTo(new Sequence(name1, seq1))));
    }
}
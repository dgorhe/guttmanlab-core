package guttmanlab.core.sequence;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import org.apache.log4j.Logger;

import htsjdk.samtools.util.CloseableIterator;

/**
 * Abstract class for parsing sequence files such as FASTA and FASTQ files.
 * <p><p>
 * Contains most of the shared iterator logic.
 * @param <T>  the type of records this parser understands (e.g., FastqSequence)
 */
public abstract class SequenceParser<T extends Sequence> implements CloseableIterator<T> {

    protected final Path p;
    protected final BufferedReader br;
    protected T next;
    
    public SequenceParser(Path p) throws IOException {
        if (p == null) {
            throw new IllegalArgumentException("SequenceParser constructed with"
                    + " null path.");
        }
        
        this.p = p;
        br = Files.newBufferedReader(p, StandardCharsets.US_ASCII);
    }
    
    /**
     * Gets the logger of this SequenceParser.
     * <p><p>
     * A workaround method to allow the use of the static logger of derived
     * classes within code of the SequenceParser superclass (for example,
     * in close()).
     */
    protected abstract Logger getLogger();
    
    @Override
    public void close() {
        try {
            br.close();
        } catch (IOException e) {
            getLogger().error("Exception caught when closing reader.", e);
        };
    }
    
    @Override
    public boolean hasNext() {
        return next != null;
    }
    
    @Override
    public T next() {
        T rtrn = next;
        findNext();
        return rtrn;
    }
    
    /**
     * Gets the next Sequence in the file.
     * <p><p>
     * A helper method which contains the file-format-specific parsing logic.
     * Advances the iteration in the file and stores the Sequence in a field
     * to be retrieved when next() is called.
     */
    protected abstract void findNext();
}
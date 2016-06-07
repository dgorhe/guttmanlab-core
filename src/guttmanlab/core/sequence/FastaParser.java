package guttmanlab.core.sequence;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import org.apache.log4j.Logger;

/**
 * A class for parsing FASTA files.
 * <p><p>
 * This class implements AutoCloseable, and is meant to be used in a try-with-
 * resources statement when used as an iterator:
 * <pre>
 * {@code
 * try (FastaParser p = new FastaParser(path)) {
 *     while (p.hasNext()) {
 *         doStuff();
 *     }
 * }
 * </pre>
 */
public final class FastaParser implements Iterator<Sequence>, AutoCloseable {

	private final Path p;
    private final BufferedReader br;
    private Sequence next;
    private final static int NUM_FASTA_LINES = 2;
    
    private static final Logger logger = Logger.getLogger(FastaParser.class);
    
    public FastaParser(Path p) throws IOException {
        if (p == null) {
            throw new IllegalArgumentException("FastqParser constructed with"
                    + " null path.");
        }

        this.p = p;
        br = Files.newBufferedReader(p, StandardCharsets.US_ASCII);
        findNext();
    }
    
    @Override
    public void close() throws IOException {
        br.close();
    }

    @Override
    public boolean hasNext() {
        return next != null;
    }

    private void findNext() {
        String[] s = new String[NUM_FASTA_LINES];
        String line = null;
        int i = 0;
        try {
            while (i < NUM_FASTA_LINES && (line = br.readLine()) != null) {
                s[i] = line;
                i++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        if (s[1] != null) {
            next = new Sequence(s[0].substring(1), s[1]);
        } else {
            next = null;
            logger.warn("FASTA file " + p.toString() + " has an incomplete final record.");
        }
    }
    
    @Override
    public Sequence next() {
        Sequence rtrn = next;
        findNext();
        return rtrn;
    }
    
    /**
     * Reads an entire FASTA file into memory.
     * <p><p>
     * This method will not throw an exception if the final record is truncated.
     * Rather, it will log a warning.
     * @param p  the path of the FASTA file
     * @return the FASTA records within the file
     * @throws IOException if an I/O error occurs opening the file
     */
    public static Collection<Sequence> slurp(Path p) throws IOException {

    	logger.info("Reading sequences from FASTA file " + p.toString());

    	Collection<Sequence> rtrn = new ArrayList<Sequence>();
		try (BufferedReader br = Files.newBufferedReader(p, StandardCharsets.US_ASCII)) {
			boolean firstLine = true;
			String line1 = "";
			String line2 = "";
			String line = null;
			while ((line = br.readLine()) != null) {
				if (firstLine) {
					line1 = line.substring(1); // Remove the initial '>'
				} else {
					line2 = line;
					rtrn.add(new Sequence(line1, line2));
				}
				firstLine = !firstLine;
			}
			if (!firstLine) {
				logger.warn("FASTA file " + p.toString() + " has an incomplete final record.");
			}
			
		}
		return rtrn;
    }
}
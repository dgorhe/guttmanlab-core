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

public final class FastqParser implements Iterator<FastqSequence>, AutoCloseable {

	private final Path p;
    private final BufferedReader br;
    private FastqSequence nextFastq;
    private final PhredEncoding pe;
    private final static int NUM_FASTQ_LINES = 4;
    
    private static final Logger logger = Logger.getLogger(FastqParser.class);

    /**
     * A class for parsing FASTQ files.
     * <p><p>
     * This class implements AutoCloseable, and is meant to be used in a try-with-
     * resources statement when used as an iterator:
     * <pre>
     * {@code
     * try (FastqParser p = new FastqParser(path, PhredEncoding.SANGER)) {
     *     while (p.hasNext()) {
     *         doStuff();
     *     }
     * }
     * </pre>
     */
    public FastqParser(Path p, PhredEncoding pe) throws IOException {
        if (p == null) {
            throw new IllegalArgumentException("FastqParser constructed with"
                    + " null path.");
        }
        
        if (pe == null) {
            throw new IllegalArgumentException("FastqParser constructed with"
                    + " null phred encoding.");
        }
        this.p = p;
        br = Files.newBufferedReader(p, StandardCharsets.US_ASCII);
        this.pe = pe;
        findNext();
    }
    
    @Override
    public void close() throws IOException {
        br.close();
    }

    @Override
    public boolean hasNext() {
        return nextFastq != null;
    }

    private void findNext() {
        String[] s = new String[NUM_FASTQ_LINES];
        String line = null;
        int i = 0;
        try {
            while (i < NUM_FASTQ_LINES && (line = br.readLine()) != null) {
                s[i] = line;
                i++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        if (s[3] != null) {
            nextFastq = new FastqSequence(s[0].substring(1), s[1], s[3], pe);
        } else {
            nextFastq = null;
            logger.warn("FASTQ file " + p.toString() + " has an incomplete final record.");
        }
    }
    
    @Override
    public FastqSequence next() {
        FastqSequence rtrn = nextFastq;
        findNext();
        return rtrn;
    }

    /**
     * Reads an entire FASTQ file into memory.
     * <p><p>
     * This method will not throw an exception if the final record is truncated.
     * Rather, it will log a warning. However, if the final record is truncated
     * in such a way that the sequence length differs from the Phred quality
     * length, the FastqSequence constructor will throw an exception.
     * @param p  the path of the FASTQ file
     * @param pe the Phred encoding scheme of the FASTQ file
     * @return the FASTQ records within the file
     * @throws IOException if an I/O error occurs opening the file
     */
    public static Collection<Sequence> slurp(Path p, PhredEncoding pe) throws IOException {

    	logger.info("Reading sequences from FASTQ file " + p.toString());

    	Collection<Sequence> rtrn = new ArrayList<Sequence>();
		try (BufferedReader br = Files.newBufferedReader(p, StandardCharsets.US_ASCII)) {
			int lineNum = 1;
			String line1 = "";
			String line2 = "";
			String line4 = "";
			String line = null;
			while ((line = br.readLine()) != null) {
				switch (lineNum) {
				case 1:
					line1 = line.substring(1); // Remove the initial '@'
					break;
				case 2:
					line2 = line;
					break;
				case 3:
					break;
				case 4:
					line4 = line;
					rtrn.add(new FastqSequence(line1, line2, line4, pe));
					break;
				default:
					logger.fatal("Unrecoverable error in FASTQ parsing.", new IllegalStateException());
					System.exit(1);
				}
				lineNum = lineNum == 4 ? 1 : lineNum + 1;
			}
			if (lineNum != 1) {
				logger.warn("FASTQ file " + p.toString() + " has an incomplete final record.");
			}
			
		}
		return rtrn;
    }
}
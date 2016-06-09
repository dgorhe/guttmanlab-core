package guttmanlab.core.sequence;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;

/**
 * A class for parsing FASTA files.
 * <p><p>
 * This class implements CloseableIterator, and is meant to be used in a try-
 * with-resources statement:
 * <pre>
 * {@code
 * try (FastaParser p = new FastaParser(path)) {
 *     while (p.hasNext()) {
 *         doStuff();
 *     }
 * }
 * </pre>
 */
public final class FastaParser extends SequenceParser<Sequence> {

    private static final int NUM_FASTA_LINES = 2;
    private static final Logger logger = Logger.getLogger(FastaParser.class);
    
    public FastaParser(Path p) throws IOException {
        super(p);
        findNext();
    }

    @Override
    protected void findNext() {
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

        checkHeapSize(p);
        
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
    
    /**
     * Checks if there's enough heap memory allocated to handle reading in the
     * file at the given Path.
     * <p><p>
     * Throws an error if the maximum heap size is greater than the size of the file.
     * @param p  path of the file
     * @throws OutOfMemoryError if there isn't enough heap memory.
     */
    private static void checkHeapSize(Path p) {
        long maxHeapSize = Runtime.getRuntime().maxMemory();
        long fileSize = 0;

        try {
            fileSize = Files.size(p);
        } catch (IOException e) {
            logger.error("Unable to check file size of " + p.toString(), e);
        }

        if (fileSize > maxHeapSize) {
            String msg = "File " + p.toString() + " is too large. File"
                    + "size: " + fileSize + " bytes. Max heap size: "
                    + maxHeapSize + " bytes.";
            logger.fatal(msg);
            throw new OutOfMemoryError(msg);  // What's the best way to deal with this?
        }
    }

    @Override
    protected Logger getLogger() {
        return logger;
    }
}
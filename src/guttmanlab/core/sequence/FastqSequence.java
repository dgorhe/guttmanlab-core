package guttmanlab.core.sequence;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import org.apache.commons.lang.builder.HashCodeBuilder;

/**
 * A representation of a FASTQ record.
 * <p><p>
 * A FastqSequence is simply a sequence of nucleotides bundled with Phred
 * quality scores and a name.
 */
public final class FastqSequence extends Sequence {

    // byte[] instead of String to handle different Phred encodings
    private final byte[] quality;
    
    /**
     * FastqSequence constructor.
     * <p><p>
     * This constructor does not check the seq parameter. Any String can be
     * passed as a sequence of bases.
     * @param name  the name of this FastqSequence
     * @param seq  the bases of this FastqSequence
     * @param quality  the Phred quality scores of this FastqSequence
     * @param pe  the Phred encoding of the quality scores
     * @throws IllegalArgumentException if null is passed as an argument
     * @throws IllegalArgumentException if the lengths of the bases and the
     * quality scores are not equal
     */
    public FastqSequence(String name, String seq, String quality,
            PhredEncoding pe) {
        super(name, seq);

        if (quality == null) {
            throw new IllegalArgumentException("Attemped to create"
                    + " FastqSequence with null quality String.");
        }
        
        if (pe == null) {
            throw new IllegalArgumentException("Attemped to create"
                    + " FastqSequence with null PhredEncoding.");
        }

        if (seq.length() != quality.length()) {
            throw new IllegalArgumentException("Sequence bases and quality"
                    + " scores do not agree.");
        }
        
        byte[] tmp = quality.getBytes(StandardCharsets.US_ASCII);
        for (int i = 0; i < tmp.length; i++) {
            tmp[i] -= pe.offset();
        }
        
        this.quality = tmp;
    }
    
    /**
     * FastqSequence constructor.
     * <p><p>
     * This constructor defaults to a Sanger Phred encoding. It also does not
     * check the seq parameter. Any String can be passed as a sequence of
     * bases.
     * @param name  the name of this FastqSequence
     * @param seq  the bases of this FastqSequence
     * @param quality  the Phred quality scores of this FastqSequence
     * @throws IllegalArgumentException if null is passed as an argument
     * @throws IllegalArgumentException if the lengths of the bases and the
     * quality scores are not equal
     */
    public FastqSequence(String name, String seq, String quality) {
    	this(name, seq, quality, PhredEncoding.SANGER);
    }
    
    /*
     * Constructor passing the quality as raw byte scores instead of a string.
     * Used internally.
     */
    private FastqSequence(String name, String seq, byte[] quality) {
        super(name, seq);
        this.quality = quality; // does not copy the array!
    }

    @Override
    public FastqSequence changeName(String name) {
        return new FastqSequence(name, sequence, Arrays.copyOf(quality, quality.length));
    }
    
    @Override
    public FastqSequence reverseComplement() {
        return reverseComplement(name);
    }

    @Override
    public FastqSequence reverseComplement(String name) {
        byte[] reverseQuals = new byte[quality.length];
        for (int i = 0; i < reverseQuals.length; i++) {
            reverseQuals[i] = quality[quality.length - i - 1];
        }
        return new FastqSequence(name, reverse(complement(sequence)),
                reverseQuals);
    }
    
    @Override
    public FastqSequence subsequence(int start, int end) {
        return subsequence(name, start, end);
    }
    
    @Override
    public FastqSequence subsequence(String name, int start, int end) {
        String subseq = sequence.substring(start, end);
        byte[] subqual = Arrays.copyOfRange(quality, start, end);
        return new FastqSequence(name, subseq, subqual);
    }
    
    /**
     * Converts this Sequence into a String representation suitable for
     * outputting to a FASTQ file. This String is not newline terminated.
     * @param pe  the Phred encoding to apply to the quality scores of
     * this FastqSequence
     */
    public String toFormattedString(PhredEncoding pe) {
        byte[] tmp = new byte[quality.length];
        for (int i = 0; i < quality.length; i++) {
            tmp[i] = (byte) (quality[i] + pe.offset());
        }
        return "@" + name + nl + sequence + nl + "+" + nl + new String(tmp);
    }

    /**
     * Converts this Sequence into a String representation suitable for
     * outputting to a FASTQ file. Assumes a Sanger Phred encoding. This
     * String is not newline terminated.
     */
    @Override
    public String toFormattedString() {
        return toFormattedString(PhredEncoding.SANGER);
    }
    
    @Override
    public boolean equals(Object other) {

    	if (this == other) {
    		return true;
    	}
    	
        // No need for null check. The instanceof operator returns false if
        // (other == null).
        if (!(other instanceof FastqSequence)) {
            return false;
        }

        // OK to cast this. Class was explicitly checked above
        FastqSequence o = (FastqSequence) other;
        
        if (quality.length != o.quality.length) {
            return false;
        }
        
        boolean qualityEquals = true;
        for (int i = 0; i < quality.length; i++) {
            if (quality[i] != o.quality[i]) {
                qualityEquals = false;
            }
        }
        return qualityEquals && name.equals(o.name)
                && sequence.equals(o.sequence);
    }
    
    @Override
    public int hashCode() {
        return new HashCodeBuilder(17, 37)
                .append(name)
                .append(sequence)
                .append(quality)
                .toHashCode();
    }
}
package guttmanlab.core.sequence;

import org.apache.commons.lang.builder.HashCodeBuilder;

/**
 * A representation of a FASTQ record.
 * <p><p>
 * A FastqSequence is simply a sequence of nucleotides bundled with Phred
 * quality scores and a name.
 */
public final class FastqSequence extends Sequence {

    private final PhredEncoding pe;
    private final String quality;
    
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
        
        this.quality = quality;
        this.pe = pe;
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

    @Override
    public FastqSequence changeName(String name) {
        return new FastqSequence(name, sequence, quality, pe);
    }
    
    @Override
    public FastqSequence reverseComplement() {
        return reverseComplement(name);
    }

    @Override
    public FastqSequence reverseComplement(String name) {
        return new FastqSequence(name, reverse(complement(sequence)),
                reverse(quality));
    }
    
    @Override
    public FastqSequence subsequence(int start, int end) {
        return subsequence(name, start, end);
    }
    
    @Override
    public FastqSequence subsequence(String name, int start, int end) {
        String subseq = sequence.substring(start, end);
        String subqual = quality.substring(start, end);
        return new FastqSequence(name, subseq, subqual, pe);
    }
    
    /**
     * Converts this Sequence into a String representation suitable for
     * outputting to a FASTQ file. This String is not newline terminated.
     */
    @Override
    public String toFormattedString() {
        return "@" + name + nl + sequence + nl + "+" + nl + quality;
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
        
        return quality.equals(o.quality) && name.equals(o.name)
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
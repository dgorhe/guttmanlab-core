package guttmanlab.core.sequence;

import org.apache.commons.lang.builder.HashCodeBuilder;

/**
 * A representation of a sequence of nucleotides.
 * <p><p>
 * A Sequence is simply a sequence of nucleotides bundled with a name.
 */
public class Sequence {
	
	protected final String sequence;
	protected String name;
	protected static final String nl = "\n";
	
	/**
	 * {@code Sequence} constructor.
	 * <p>
	 * This constructor does not check the seq parameter. Any String can be
	 * passed as a sequence of bases.
	 * @param name  the name of this Sequence
	 * @param seq   the bases of this Sequence
	 * @throws IllegalArgumentException if null is passed as an argument
	 */
	public Sequence(String name, String seq) {
		if (name == null) {
			throw new IllegalArgumentException("Attempted to construct a Sequence"
					+ " object with a null name");
		}
		if (seq == null) {
			throw new IllegalArgumentException("Attempted to construct a Sequence"
					+ " object with a null sequence");
		}
		this.name = name;
		this.sequence = seq;
	}

	/**
	 * @return the bases of this Sequence
	 */
	public String bases() {
		return sequence;
	}
	
	/**
	 * @return the name of this Sequence
	 */
	public String name() {
		return name;
	}
	
	/**
	 * Constructs a new Sequence which is identical to this Sequence, except
	 * with a different name. This method creates a new Sequence rather than
	 * modifying this Sequence in order to maintain immutability.
	 * @param name  the new name
	 * @return a Sequence with a different name but otherwise identical
	 */
    public Sequence changeName(String name) {
        return new Sequence(name, sequence);
    }
	
	/**
	 * Creates a new Sequence object in which the bases have been reverse-
	 * complemented. The name of the new Sequence remains the same.
	 * @return the reverse-complement of this Sequence
	 */
	public Sequence reverseComplement() {
		return reverseComplement(name);
	}


	/**
	 * Creates a new Sequence object in which the bases have been reverse-
     * complemented. 
	 * @param name  the name of the new Sequence
	 * @return the reverse-complement of this Sequence
	 */
	public Sequence reverseComplement(String name) {
		return new Sequence(name, reverse(complement(sequence)));
	}
	
	protected String reverse(String s) {
		return (new StringBuilder(s).reverse().toString());
	}
	
	protected String complement(String s) {
		char[] cs = s.toCharArray();
		char[] rtrn = new char[cs.length];
		for (int i = 0; i < cs.length; i++) {
			rtrn[i] = complement(cs[i]);
		}
		return String.valueOf(rtrn);
	}
	
	/**
	 * Complements a base.
	 * <p>
	 * Case is retained (for example, {@code complement('a') == 't'}).
	 * Only characters in [ACGTNacgtn] are supported.
	 * @param c  the base to complement
	 * @return the complement of this base
	 * @throws IllegalArgumentException if passed an unsupported base
	 */
	protected char complement(char c) {
		switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'T': return 'A';
		case 'N': return 'N';
		case 'a': return 't';
		case 'c': return 'g';
		case 'g': return 'c';
		case 't': return 'a';
		case 'n': return 'n';
		default: throw new IllegalArgumentException("Unsupported base: " + c);
		}
	}
	
	/**
     * Creates a subsequence of this Sequence as a new Sequence object.
     * The name of the new Sequence remains the same. The coordinates
     * are zero-based. The new Sequence begins at the specified start
     * coordinate and extends to the base at the specified end
     * coordinate minus 1. 
     * @param start  the start coordinate
     * @param end    the end coordinate
     * @return a subsequence of this Sequence
     */
	public Sequence subsequence(int start, int end) {
		return subsequence(name, start, end);
	}

	/**
     * Creates a subsequence of this Sequence as a new Sequence object
     * with a different name. The coordinates are zero-based. The new
     * Sequence begins at the specified start coordinate and extends to
     * the base at the specified end coordinate minus 1.
     * @param name   the new name 
     * @param start  the start coordinate
     * @param end    the end coordinate
     * @return a subsequence of this Sequence
     */
	public Sequence subsequence(String name, int start, int end) {
		String subseq = sequence.substring(Math.max(start, 0),
				Math.min(end, sequence.length()));
		return new Sequence(name, subseq);
	}
	
	/**
	 * @return the length of the sequence of bases in this Sequence
	 */
	public int length() {
		return sequence.length();
	}
	
	/**
	 * @return if this Sequence consists entirely of A's or entirely of
	 * T's, case-insensitive
	 */
	public boolean isPolyA() {
		return sequence.chars().allMatch(c -> c == 'a' || c == 'A') ||
				sequence.chars().allMatch(c -> c == 't' || c == 'T');
	}
	
	/**
	 * Converts this Sequence into a String representation suitable for
	 * outputting to a FASTA file. This String is not newline terminated.
	 */
	public String toFasta() {
		return ">" + name + nl + sequence;
	}

    /**
     * Converts this Sequence into a String representation suitable for
     * outputting to a FASTA file. This String is not newline terminated.
     */
	public String toFormattedString() {
	    return toFasta();
	}
	
	/**
	 *  Returns a String representation of this Sequence. The exact details of
	 *  this representation are unspecified and subject to change, but the
	 *  following may be regarded as typical: "name:sequence"
	 */
	@Override
	public String toString() {
		return name + ": " + sequence;
	}
	
	@Override
	public boolean equals(Object o) {
		
		if (this == o) {
			return true;
		}
		
		if (!(o instanceof Sequence)) {
			return false;
		}
		
		Sequence other = (Sequence) o;

		if (!name().equals(other.name()))	{
			return false;
		}
		if (!bases().equals(other.bases()))	{
			return false;
		}
		
		return true;
	}
	
	@Override
	public int hashCode() {
		return new HashCodeBuilder(17, 31)
			.append(name)
			.append(sequence)
			.hashCode();
	}	
}
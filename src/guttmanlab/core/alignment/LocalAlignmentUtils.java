package guttmanlab.core.alignment;

public class LocalAlignmentUtils {

	/**
	 * Whether there is a full length ungapped match of sequence 2 to sequence 1, starting at specified position of sequence 1, allowing mismatches
	 * @param seq1 Sequence 1
	 * @param seq2 Sequence 2
	 * @param seq1start Required start position of match on sequence 1
	 * @param ignoreCase Ignore case
	 * @return Whether the best ungapped alignment is full length and has at most the max number of mismatches
	 */
	public static boolean containsFullLengthUngappedMatchAtPosition(String seq1, String seq2, int seq1start, int maxMismatches, boolean ignoreCase) {
		int s2len = seq2.length();
		if(seq1.length() - seq1start > s2len) {
			throw new IllegalArgumentException("Seq1 not long enough");
		}
		return hammingDist(seq1.substring(seq1start, s2len), seq2, ignoreCase) <= maxMismatches;
	}

	/**
	 * Hamming distance between two strings
	 * @param seq1 String 1
	 * @param seq2 String 2
	 * @param ignoreCase Ignore case
	 * @return Hamming distance
	 */
	public static int hammingDist(String seq1, String seq2, boolean ignoreCase) {
		if(seq1.length() != seq2.length()) {
			throw new IllegalArgumentException("Sequences must have same length");
		}
		char[] s1 = ignoreCase ? seq1.toUpperCase().toCharArray() : seq1.toCharArray();
		char[] s2 = ignoreCase ? seq2.toUpperCase().toCharArray() : seq2.toCharArray();
		int rtrn = 0;
		for(int i = 0; i < s1.length; i++) {
			if(s1[i] != s2[i]) rtrn++;
		}
		return rtrn;
	}

}

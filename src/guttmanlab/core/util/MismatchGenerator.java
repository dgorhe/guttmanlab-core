package guttmanlab.core.util;

import java.util.Collection;
import java.util.HashSet;

/**
 * Introduce mismatches to a sequence
 * @author prussell
 *
 */
public final class MismatchGenerator {

	private static final char[] alphabet = {'A', 'C', 'G', 'T'};

	/**
	 * Prohibit instantiation
	 */
	private MismatchGenerator(){}
	
	/**
	 * Get the chars in the alphabet except this one
	 * @param c Any char. Throws exception if not in alphabet.
	 * @return Chars in the alphabet not equal to this one
	 */
	private static char[] otherChars(char c) {
		char[] rtrn = new char[alphabet.length - 1];
		boolean found = false;
		int curr = 0;
		int i = 0;
		while(i < alphabet.length) {
			if(c == alphabet[i]) {
				found = true;
				i++;
				continue;
			} else {
				rtrn[curr] = alphabet[i];
				curr++;
				i++;
			}
		}
		if(!found) {
			throw new IllegalArgumentException("Char " + c + " not in alphabet");
		}
		return rtrn;
	}

	/**
	 * Systematically replace each position with each other char in alphabet
	 * @param seq Sequence
	 * @return Set of sequences with hamming distance one from original sequence
	 */
	private static Collection<String> introduceOneMismatch(String seq) {
		char[] charSeq = seq.toUpperCase().toCharArray();
		Collection<String> rtrn = new HashSet<String>();
		for(int i = 0; i < charSeq.length; i++) {
			char realChar = charSeq[i];
			char[] otherChars = otherChars(realChar);
			for(int j = 0; j < otherChars.length; j++) {
				charSeq[i] = otherChars[j];
				rtrn.add(new String(charSeq));
			}
			charSeq[i] = realChar;
		}
		return rtrn;
	}

	/**
	 * Get a set containing all possible strings with hamming distance one from a string in the input set
	 * No strings in return set will have hamming distance > 1 from any string in the input set
	 * Some strings in return set may have hamming distance 0 from a string in input set
	 * @param seqs Input strings
	 * @return Set of all strings of hamming distance 1 from a string in input set, plus maybe some with hamming distance 0
	 */
	private static Collection<String> introduceAllPossibleSingleMismatches(Collection<String> seqs) {
		Collection<String> rtrn = new HashSet<String>();
		for(String seq : seqs) {
			rtrn.addAll(introduceOneMismatch(seq));
		}
		return rtrn;
	}
	
	/**
	 * Get all possible versions of the sequence with this many mismatches, plus maybe some with fewer mismatches
	 * @param seq Sequence
	 * @param mismatches Number of mismatches
	 * @return A set containing all versions of the sequence with the requested number of mismatches, plus maybe some with fewer mismatches
	 */
	public static Collection<String> getRepresentatives(String seq, int mismatches) {
		if(mismatches == 0) {
			Collection<String> rtrn = new HashSet<String>();
			rtrn.add(seq);
			return rtrn;
		}
		Collection<String> prev = new HashSet<String>();
		prev.add(seq);
		int mismatchesDone = 0;
		while(true) {
			Collection<String> mutated = introduceAllPossibleSingleMismatches(prev);
			mismatchesDone++;
			if(mismatchesDone == mismatches) {
				return mutated;
			}
			prev.clear();
			prev.addAll(mutated);
		}
	}
	
	/**
	 * Get all possible versions of the sequences with this many mismatches, plus maybe some with fewer mismatches
	 * @param seq Sequences
	 * @param mismatches Number of mismatches
	 * @return A set containing all versions of the sequences with the requested number of mismatches, plus maybe some with fewer mismatches
	 */
	public static Collection<String> getRepresentatives(Collection<String> seqs, int mismatches) {
		Collection<String> rtrn = new HashSet<String>();
		for(String seq : seqs) {
			rtrn.addAll(getRepresentatives(seq, mismatches));
		}
		return rtrn;
	}

}

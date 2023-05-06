package guttmanlab.core.probegeneration;


import java.util.List;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.lang.StringUtils;






/**
 * @author engreitz
 * Filter out low-complexity sequences.  Simple algorithm:  search probe sequences for perfect matches of
 * length n to all possible low-complexity repeats of length 1 to k by creating a hash of sequences of length n
 * from low-complexity repeats and comparing each n-mer in the probe sequence to this hash.
 */
public class LowComplexityFilter implements ProbeFilter  {

	public String name = "LowComplexity";
	
	// Setting K = 6 will check for sequencing containing repeats of up to k = 6 bases (e.g., (GGCACA GGCACA ...))
	private int maxK = 6;
	
	// Sets the length of the repeat that must be present to call a match
	private int[] minLength = new int[] { 10, 10, 12, 12, 18, 18 };
	
	private Map<Integer, HashSet<String>> hashes = null;
	
	
	
	/**
	 * Constructor with default settings
	 */
	public LowComplexityFilter() {
		hashes = getLowComplexityHashes(maxK, minLength);
	}
	
	
	/**
	 * @param k
	 * @param minLengths
	 * @return
	 * Set up hash sets containing all strings of minLengths (n) bases from k-mer repeats
	 */
	private Map<Integer, HashSet<String>> getLowComplexityHashes(int k, int[] minLengths) {
		Map<Integer, HashSet<String>> hashes = new TreeMap<Integer, HashSet<String>>();
		for (int repeatK = 0; repeatK < k; repeatK++) {
			HashSet<String> currHash = hashes.containsKey(minLengths[repeatK]) ?
				currHash = hashes.get(minLengths[repeatK]) : new HashSet<String>();
	
			List<String> allRepeats = new ArrayList<String>();
			addAllKmerRepeats(allRepeats, "", repeatK+1);
			for (String repeatElement : allRepeats) {
				addRepeatToHash(currHash, repeatElement, minLengths[repeatK]);
			}
			
			hashes.put(minLengths[repeatK], currHash);
		}
		return hashes;
	}
	
	
	final char chars[] = new char[] {'A','C','T','G'};
	/**
	 * @param allRepeats
	 * @param base
	 * @param k
	 * Recursive function to generate all possible k-mer sequences
	 */
	private void addAllKmerRepeats(List<String> allRepeats, String base, int k) {
		if (base.length() >= k) {
			//logger.info("Adding " + base + " to list.");
			allRepeats.add(base);
		} else {
			for (int i = 0; i < chars.length; i++)
				addAllKmerRepeats(allRepeats, base + chars[i], k);
		}
	}
	
	
	/**
	 * @param hash
	 * @param repeatElement
	 * @param hashLength
	 * Take a repeat element (e.g., GGA) and add all strings of hashLength bases to the hash
	 */
	private void addRepeatToHash(HashSet<String> hash, String repeatElement, int hashLength) {
		repeatElement = StringUtils.repeat(repeatElement, hashLength*2);
		//logger.info("Adding " + repeatElement);
		for (int i = 0; i < hashLength; i++) {
			hash.add(repeatElement.substring(i,hashLength));
		}
	}

	


	@Override
	public boolean rejectProbe(Probe probe) {
		return rejectSequence(probe.getProbeSequence().getSequenceBases());
	}
	

	public boolean rejectSequence(String s) {
		for (int period = 0; period < maxK; period++) {
			if (sequenceContainsNmerInHash(s, minLength[period])) return true;
		}
		return false;
	}


	private boolean sequenceContainsNmerInHash(String seq, int n) {
		for (int charIndex = 0; charIndex < seq.length() - n + 1; charIndex++) {
			String substr = seq.substring(charIndex, charIndex + n);
			if (hashes.get(n).contains(substr)) return true;
		}
		return false;
	}
	
	
	
	
	
}

package guttmanlab.core.probegeneration;



/**
 * @author prussell
 * Filter probes with too many repeat masked bases
 */
public class RepeatFilter implements ProbeFilter {
	
	private boolean includesLower;
	private boolean includesN;
	private double maxPct;
	
	/**
	 * 
	 */
	public RepeatFilter() {}
	
	/**
	 * @param maxRepeatPct Maximum allowable percentage of repeat masked bases (0 to 1)
	 * @param includeLowerCase Lower case bases count as repeats
	 * @param includeN Ns count as repeats
	 */
	public RepeatFilter(double maxRepeatPct, boolean includeLowerCase, boolean includeN) {
		if(maxRepeatPct < 0 || maxRepeatPct > 1) {
			throw new IllegalArgumentException("Max repeat percentage must be between 0 and 1");
		}
		if(!includeLowerCase && !includeN) {
			throw new IllegalArgumentException("Must include at least one: lowercase or N");
		}
		maxPct = maxRepeatPct;
		includesLower = includeLowerCase;
		includesN = includeN;
	}
	
	
	public boolean rejectSequence(String probeSeq) {
		int size = probeSeq.length();
		int numRepeats = 0;
		for(int i=0; i<size; i++) {
			char c = probeSeq.charAt(i);
			if(includesLower && !Character.isUpperCase(c)) {
				numRepeats++;
				continue;
			}
			if(includesN && (c == 'n' || c == 'N')) {
				numRepeats++;
				continue;
			}
		}
		double pct = (double) numRepeats / (double) size;
		return pct > maxPct;
	}
	

	@Override
	public boolean rejectProbe(Probe probe) {
		String probeSeq = probe.getProbeSequence().getSequenceBases();
		int size = probeSeq.length();
		int numRepeats = 0;
		for(int i=0; i<size; i++) {
			char c = probeSeq.charAt(i);
			if(includesLower && !Character.isUpperCase(c)) {
				numRepeats++;
				continue;
			}
			if(includesN && (c == 'n' || c == 'N')) {
				numRepeats++;
				continue;
			}
		}
		double pct = (double) numRepeats / (double) size;
		return pct > maxPct;
	}
	
	private static String INCLUDE_N_FLAG = "N";
	private static String INCLUDE_LOWERCASE_FLAG = "lower";

	
}

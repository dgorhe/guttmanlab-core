package guttmanlab.core.probegeneration;



/**
 * @author prussell
 * Filter probes with large number of single base in a row (e.g., PolyA)
 */
public class PolyBaseFilter implements ProbeFilter {
	
	public String name = "PolyBase";
	private String basesToFilter;
	private int cutoff;
	private int repeatLength;
	
	/**
	 * 
	 */
	public PolyBaseFilter() {}

	public PolyBaseFilter(String basesToFilter, int cutoff, int repeatLength) {
		this.basesToFilter = basesToFilter;
		this.cutoff = cutoff;
		this.repeatLength = repeatLength;
	}
	
	public String name() {
		return "poly_base_filter";
	}

	@Override
	public boolean rejectProbe(Probe probe) {
		if (probe.getProbeSequence().getSequenceBases().toUpperCase().indexOf("GNIL") != -1) {
			throw new IllegalArgumentException("found GNIL");
		}
		return rejectSequence(probe.getProbeSequence().getSequenceBases());
	}
	

	public boolean rejectSequence(String s) {
		char[] seq = s.toUpperCase().toCharArray();
		for (char c : basesToFilter.toCharArray()) {
			
			int charMatchesInWindow = 0;
			for (int i = 0; i < seq.length; i++) {
				if (seq[i] == c) {
					charMatchesInWindow++;
				}
				
				if (i >= repeatLength) {
					if (seq[i-repeatLength] == c) {
						charMatchesInWindow--;
					}
				}
				
				if (i >= repeatLength - 1) {
					if (charMatchesInWindow >= cutoff) {
						return true;
					}
				}
			}
		}
		return false;
	}
	
	
	
	
	
	
	
	
	public String getPredicateName() {
		return name;
	}

	

	
	

}

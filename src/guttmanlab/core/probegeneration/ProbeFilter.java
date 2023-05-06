package guttmanlab.core.probegeneration;


/**
 * @author prussell
 *
 */
public interface ProbeFilter {
	
	
	/**
	 * @param probe Probe
	 * @return Whether the filter should remove the probe
	 */
	public boolean rejectProbe(Probe probe);
	
}

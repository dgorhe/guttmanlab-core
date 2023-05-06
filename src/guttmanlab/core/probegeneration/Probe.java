package guttmanlab.core.probegeneration;

import guttmanlab.core.sequence.Sequence;

/**
 * @author prussell
 * The probe itself without primers
 */
public class Probe implements Comparable<Probe> {
	
	private Sequence parent = null;
	private String parentSequenceId;
	private String id;
	private int start;
	private int size;
	private boolean antisense;
	
	
	/**
	 * @param parentTranscript Parent transcript
	 * @param layout Probe layout that produced the probe
	 * @param probeID Unique probe identifier
	 * @param startPosOnTranscript First transcript position included in probe
	 * @param probeSize Probe length
	 * @param antisenseToTranscript Whether the probe is antisense to transcript
	 */
	public Probe(Sequence parentTranscript, String probeID, int startPosOnTranscript, int probeSize, boolean antisenseToTranscript) {
		parent = parentTranscript;
		id = probeID;
		start = startPosOnTranscript;
		size = probeSize;
		antisense = antisenseToTranscript;
		
	}
	
	
	/**
	 * @return Probe ID
	 */
	public String getID() {
		return id;
	}
	
	/**
	 * @return Parent transcript
	 */
	public Sequence getParentTranscript() {
		return parent;
	}
	
	public String getParentTranscriptId() {
		return (parent == null) ? parentSequenceId : parent.getName();
	}
	
	/**
	 * @return Probe start position on transcript
	 */
	public int getStartPosOnTranscript() {
		return start;
	}
	
	/**
	 * @return First transcript position after probe
	 */
	public int getEndPosOnTranscript() {
		return start + size;
	}
	
	/**
	 * @return Probe size
	 */
	public int getSize() {
		return size;
	}
	
	
	
	/**
	 * @return Probe tiling path
	 */
	public int getTilingPath() {
		return start % size;
	}
	
	/**
	 * @return Whether the probe sequence is antisense to the transcript
	 */
	public boolean isAntisenseToTranscript() {
		return antisense;
	}
	
	/**
	 * @return The sequence of the probe
	 */
	public Sequence getProbeSequence() {
		String probeSeq = parent.getSubSequence(start, start + size);
		if(antisense) {
			probeSeq=Sequence.reverseComplement(probeSeq);
		}
		Sequence rtrn=new Sequence(this.getID(), probeSeq);
		return rtrn;
	}
	
	@Override
	public boolean equals(Object o) {
		if(!o.getClass().equals(getClass())) {
			return false;
		}
		Probe p = (Probe)o;
		return p.hashCode() == hashCode();
	}
	
	@Override
	public int hashCode() {
		String str = parent.getName() + "_" + id + "_" + start + "_" + size + "_" + Boolean.valueOf(antisense).toString();
		return str.hashCode();
	}
	
	@Override
	public int compareTo(Probe other) {
		if(!parent.equals(other.getParentTranscript())) {
			return parent.getName().compareTo(other.getID());
		}
		if(start == other.getStartPosOnTranscript()) {
			if(size == other.getSize()) {
				if(antisense && !other.isAntisenseToTranscript()) return -1;
				if(!antisense && other.isAntisenseToTranscript()) return 1;
			}
			return size - other.getSize();
		}
		return start - other.getStartPosOnTranscript();
	}
	
}

package guttmanlab.core.sequence;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SingleInterval;

import java.util.Collection;
import java.util.Iterator;

/**
 * A nucleotide sequence
 * @author prussell
 *
 */
public class Sequence {
	
	private String sequence;
	private String name;
	
	public Sequence(String name, String seq){
		this.name=name;
		this.sequence=seq;
	}
	
	public Sequence(String seq){
		this.sequence=seq;
	}
	
	/**
	 * @return A new sequence that is the reverse complement of this sequence
	 */
	public Sequence reverseComplement() {
		Sequence tmpSeq = new Sequence(this.sequence);
		tmpSeq.reverse();
		return tmpSeq;
	}
	
	
	private void reverse() {
		String seqBases = getSequenceBases();
		String reversedSeq = "";
		for(int j = seqBases.length() - 1; j >= 0 ; j--) {
			char c = seqBases.charAt(j);
			if('c' == c) {
				reversedSeq+='g';
			}else if ('C' == c) {
				reversedSeq+='G';
			}else if ('G' == c) {
				reversedSeq+='C';
			}else if ('g' == c) {
				reversedSeq+=('c');
			}else if ('a' == c) {
				reversedSeq+=('t');
			}else if ('A' == c) {
				reversedSeq+=('T');
			}else if('t' == c) {
				reversedSeq+=('a');
			}else if('T' == c) {
				reversedSeq+=('A');
			}else if('N'==c){
				reversedSeq+=('N');
			}else if('n'==c){
				reversedSeq+=('n');
			}else {
				reversedSeq+=(c);
			}
		}
		
		this.sequence = reversedSeq;
	}
	
	/**
	 * @return Return the sequence bases
	 */
	public String getSequenceBases() {
		return this.sequence;
	}
	
	public String getName(){
		return this.name;
	}
	
	/**
	 * @return Sequence length
	 */
	public int getLength() {
		return sequence.length();
	}
	
	/**
	 * Get subsequence
	 * @param name Name of new sequence to return
	 * @param start Start position of subsequence
	 * @param end Position after last position to include
	 * @return The subsequence
	 */
	public Sequence getSubSequence(String name, int start, int end) {
		String subSeq = sequence.substring(Math.max(start, 0), Math.min(end, sequence.length()));
		Sequence seq = new Sequence(name, subSeq);
		return seq;
	}
	
	/**
	 * Get the spliced transcribed sequence of an annotation
	 * Bases are reported in 5' to 3' direction
	 * @param annot The annotation
	 * @return Sequence with same name as annotation containing the transcribed sequence
	 */
	public Sequence getSubsequence(Annotation annot) {
		if(!annot.getOrientation().equals(Strand.POSITIVE) && !annot.getOrientation().equals(Strand.NEGATIVE)) {
			throw new IllegalArgumentException("Strand must be known");
		}
		String seq = "";
		Iterator<SingleInterval> blockIter = annot.getBlocks();
		while(blockIter.hasNext()) {
			SingleInterval block = blockIter.next();
			Sequence blockSequence = getSubSequence("", block.getReferenceStartPosition(), block.getReferenceEndPosition());
			String forwardBases = blockSequence.getSequenceBases();
			if(annot.getOrientation().equals(Strand.POSITIVE)) {
				seq += forwardBases;
			} else {
				String rc = new Sequence(forwardBases).reverseComplement().getSequenceBases();
				seq = rc + seq;
			}
		}
		return new Sequence(annot.getName(), seq);
	}
	
	/**
	 * Soft mask the specified regions
	 * Throw an exception if the annotations do not refer to this sequence
	 * @param regions The regions to mask with respect to this sequence
	 * @return A new sequence with the regions soft masked
	 */
	public Sequence softMask(Collection<Annotation> regions) {
		// TODO
		throw new UnsupportedOperationException();
	}
	
	/**
	 * Hard mask the specified regions
	 * Throw an exception if the annotations do not refer to this sequence
	 * @param regions The regions to mask with respect to this sequence
	 * @return A new sequence with the regions hard masked
	 */
	public Sequence hardMask(Collection<Annotation> regions) {
		// TODO
		throw new UnsupportedOperationException("TODO");
	}
	
	@Override
	public boolean equals(Object o) {
		if(!o.getClass().equals(Sequence.class)) {
			return false;
		}
		Sequence otherSeq = (Sequence)o;
		if(!getName().equals(otherSeq.getName()))	{
			return false;
		}
		if(!getSequenceBases().equals(otherSeq.getSequenceBases()))	{
			return false;
		}
		return true;
	}
	
	@Override
	public int hashCode() {
		String s = getName() + ":" + getSequenceBases();
		return s.hashCode();
	}
	
}

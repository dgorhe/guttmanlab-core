package guttmanlab.core.annotation;

import org.apache.commons.lang.builder.HashCodeBuilder;

import guttmanlab.core.annotation.Annotation.Strand;

/**
 * Helper methods for Annotation methods
 * @author prussell
 *
 */
class AnnotationHelper {
	
	static boolean fullyContained(SingleInterval block1, SingleInterval block2){
		//is only fully contained if block2 is equal to or a subset of block1
		if(block1.getReferenceStartPosition()<=block2.getReferenceStartPosition() && block1.getReferenceEndPosition()>=block2.getReferenceEndPosition()){
			return true;
		}
		return false;	
	}

	/**
	 * Helper method to calculate overlaps from single blocks
	 * @param block1
	 * @param block2
	 * @return whether the blocks overlap
	 */
	static boolean overlaps(SingleInterval block1, SingleInterval block2){
		int newStart=Math.max(block1.getReferenceStartPosition(), block2.getReferenceStartPosition());
		int newEnd=Math.min(block1.getReferenceEndPosition(), block2.getReferenceEndPosition());
		
		Strand consensusStrand=Annotation.Strand.consensusStrand(block1.getOrientation(), block2.getOrientation());
		if(newStart<newEnd && block1.getReferenceName().equalsIgnoreCase(block2.getReferenceName()) && !consensusStrand.equals(Annotation.Strand.INVALID)){
			return true;
		}
		
		return false;
	}
	
	/**
	 * Merge two single intervals
	 * @param block1 Interval 1
	 * @param block2 Interval 2
	 * @return Merged interval or null if they don't overlap
	 */
	static SingleInterval merge(SingleInterval block1, SingleInterval block2) {
		if(!overlaps(block1, block2)){return null;}
		
		int newStart=Math.min(block1.getReferenceStartPosition(), block2.getReferenceStartPosition());
		int newEnd=Math.max(block1.getReferenceEndPosition(), block2.getReferenceEndPosition());
		Strand consensus=Annotation.Strand.consensusStrand(block1.getOrientation(), block2.getOrientation());
		return new SingleInterval(block1.getReferenceName(), newStart, newEnd, consensus);
	}
	
	
	private static HashCodeBuilder hashCodeBuilder(Annotation a) {
		return new HashCodeBuilder(31,37).append(a.getName()).append(a.getReferenceName()).append(a.getReferenceStartPosition())
				.append(a.getReferenceEndPosition()).append(a.getOrientation()).append(a.getNumberOfBlocks());
	}
	
	/**
	 * Hash code implementation
	 * @param a An annotation
	 * @return A hash code
	 */
	static int hashCode(Annotation a)
	{
		return hashCodeBuilder(a).toHashCode();
	}

	/**
	 * Equals implementation
	 * @param a An annotation
	 * @param b Another annotation
	 * @return True iff they are equal
	 */
	static boolean equals(Annotation a, Annotation b) {
		String aName = null;
		try {aName = a.getName();} catch(NullPointerException e) {}
		String bName = null;
		try {bName = b.getName();} catch(NullPointerException e) {}
		if(aName == null && bName != null) return false;
		if(bName == null && aName != null) return true;
		if(aName != null && bName != null) {
			if(!a.getName().equals(b.getName())) return false;
		}
		return a.compareTo(b) == 0;
	}

	/**
	 * To string implementation
	 * @param a An annotation
	 * @return A string representation of the annotation
	 */
	static String toString(Annotation a) {
		return a.toBED(0,0,0);
	}

}

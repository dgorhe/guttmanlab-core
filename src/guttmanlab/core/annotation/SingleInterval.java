package guttmanlab.core.annotation;

import guttmanlab.core.annotationcollection.AnnotationCollection;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

/**
 * The SingleInterval class represents a simple contiguous interval.
 * This is the basis for all features.
 * Interval coordinates are zero-based, left-closed and right-open.
 */
public class SingleInterval implements Annotation{

	private String referenceName;
	private int startPos;
	private int endPos;
	private Strand orientation;
	private String featureName;

	/**
	 * Constructs an interval. Start and end coordinates are zero-based, left-closed and right-open.
	 * @param refName is the name of the reference
	 * @param start is the start coordinate
	 * @param end is the end coordinate
	 * @param orientation is the strandedness of the SingleInterval
	 * @param featureName is the name of this feature
	 */
	public SingleInterval(String refName, int start, int end, Strand orientation, String featureName){
		this.referenceName = refName;
		this.startPos = start;
		this.endPos = end;
		this.orientation = orientation;
		this.featureName = featureName;
	}
	
	/**
	 * Constructs an interval. Start and end coordinates are zero-based, left-closed and right-open. The feature name
	 * is simply the empty string.
	 * @param refName is the name of the reference
	 * @param start is the start coordinate
	 * @param end is the end coordinate
	 * @param orientation is the strandedness of the SingleInterval
	 */
	public SingleInterval(String refName, int start, int end, Strand orientation){
		this(refName, start, end, orientation, "");
	}
	
	/**
	 * Constructs an interval. Start and end coordinates are zero-based, left-closed and right-open. The feature name
	 * is simply the empty string. The orientation is set to Strand.UNKNOWN
	 * @param refName is the name of the reference
	 * @param start is the start coordinate
	 * @param end is the end coordinate
	 */
	public SingleInterval(String refName, int start, int end) {
		this(refName, start, end, Strand.UNKNOWN, "");
	}
	
	/**
	 * Copy constructor
	 * @param s is the single interval to copy
	 */
	public SingleInterval(SingleInterval a) {
		this(a.getReferenceName(), a.getReferenceStartPosition(), a.getReferenceEndPosition(), a.getOrientation(), a.getName());
	}
	
	/**
	 * Copy and change name
	 * @param a Single interval to copy
	 * @param name New name
	 */
	public SingleInterval(SingleInterval a, String name) {
		this(a);
		this.featureName = name;
	}

	@Override
	public String getName() {
		return this.featureName;
	}

	@Override
	public String getReferenceName() {
		return this.referenceName;
	}

	@Override
	public int getReferenceStartPosition() {
		return this.startPos;
	}

	@Override
	public int getReferenceEndPosition() {
		return this.endPos;
	}

	/**
	 * Returns an iterator over the contained blocks. (Note: This is a SingleInterval, so
	 * it only has one block. This method exists to maintain/unify the interface of other methods.)
	 * @return an iterator over the one contained block
	 */
	@Override
	public Iterator<SingleInterval> getBlocks() {
		Collection<SingleInterval> rtrn = new ArrayList<SingleInterval>();
		rtrn.add(this);
		return rtrn.iterator();
	}
	
	/**
	 * The size of this interval, measured by simply subtracting the start position from the end position.
	 * @return The size of this interval
	 */
	@Override
	public int size() {
		return endPos - startPos;
	}

	@Override
	public Strand getOrientation() {
		return this.orientation;
	}

	@Override
	public int getNumberOfBlocks() {
		return 1;
	}

	@Override
	//FIXME This should be merged with BlockedAnnotation
	public int getRelativePositionFrom5PrimeOfFeature(int referenceStart) {
		if(referenceStart>=this.getReferenceEndPosition() || referenceStart<this.getReferenceStartPosition()){return -1;} //This start position is past the feature
		int relative=referenceStart-getReferenceStartPosition();
		if(getOrientation().equals(Strand.NEGATIVE)){
			relative=size()-relative;
		}
		return relative;
	}

	@Override
	public void setOrientation(Strand orientation) {
		this.orientation = orientation;
	}

	/**
	 * Trims this block to the relative start and end position provided. Preserves the reference name
	 * and the orientation. SingleIntervals with an orientation of 'unknown', 'both', or 'invalid' are
	 * trimmed as though they have a positive orientation.
	 * @param relativeStartPosition is the new start position, relative to the old
	 * @param relativeEndPosition is the new end position, relative to the old
	 * @return a new SingleInterval with the ends appropriately trimmed
	 */
	@Override
	public SingleInterval trim(int relativeStart, int relativeEnd) {
		if (getOrientation().equals(Strand.NEGATIVE)) {
			int newEnd = getReferenceEndPosition() - relativeStart;
			int newStart = getReferenceEndPosition() - relativeEnd;
			return new SingleInterval(getReferenceName(), newStart, newEnd, getOrientation());
		} else {
			int newEnd = getReferenceStartPosition() + relativeEnd;
			int newStart = getReferenceStartPosition() + relativeStart;
			return new SingleInterval(getReferenceName(), newStart, newEnd, getOrientation());
		}
	}
	
	/**
	 * Determines if this SingleInterval contains (fully) another Annotation. This method will return false if
	 * <li>the consensus sequence is Strand.INVALID</li>
	 * <li>the reference names do not match</li>
	 * <li>the method is pass a null object</li>
	 * @param other is the other interval
	 * @return whether or not this SingleInterval overlaps with another
	 */
	@Override
	public boolean contains(Annotation other) {
		if (other == null) {
			return false;
		}
		boolean namesMatch = getReferenceName().equalsIgnoreCase(other.getReferenceName());
		boolean validConsensusStrand = Strand.consensusStrand(getOrientation(), other.getOrientation()) != Strand.INVALID;
		boolean fivePrimeOK = getReferenceStartPosition() <= other.getReferenceStartPosition();
		boolean threePrimeOK = getReferenceEndPosition() >= other.getReferenceEndPosition();

		return namesMatch && validConsensusStrand && fivePrimeOK && threePrimeOK;
	}
	
	/**
	 * Determines if this SingleInterval overlaps another SingleInterval. This method will return false if
	 * <li>the consensus sequence is Strand.INVALID</li>
	 * <li>the reference names do not match</li>
	 * <li>the method is passed a null object</li>
	 * @param other is the other interval
	 * @return whether or not this SingleInterval overlaps with another
	 */
	public boolean overlaps(SingleInterval other) {
		if (other == null) {
			return false;
		}
		int newStart = Math.max(getReferenceStartPosition(), other.getReferenceStartPosition());
		int newEnd = Math.min(getReferenceEndPosition(), other.getReferenceEndPosition());
		boolean referencesMatch = getReferenceName().equalsIgnoreCase(other.getReferenceName());
		boolean validConsensusStrand = Annotation.Strand.consensusStrand(getOrientation(), other.getOrientation()) != Strand.INVALID;

		return newStart < newEnd && referencesMatch && validConsensusStrand;
	}
	
	/**
	 * This method does not work. Do not use.
	 */
	@Override // FIXME
	public AnnotationCollection<DerivedAnnotation<? extends Annotation>> getWindows(
			int windowSize, int stepSize) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public String toString(){
		return AnnotationHelper.toString(this);
	}
	
	
	
	@Override
	public boolean equals(Object other)
	{
		if(!(other instanceof Annotation)) {
			return false;
		}
		return AnnotationHelper.equals(this, (Annotation)other);
	}
	
	@Override
	public int hashCode()
	{
		return AnnotationHelper.hashCode(this);
	}


	
}

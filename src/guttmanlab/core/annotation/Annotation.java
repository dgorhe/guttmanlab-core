package guttmanlab.core.annotation;

import guttmanlab.core.annotationcollection.AnnotationCollection;

import java.util.Collection;
import java.util.Iterator;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * 
 * @author prussell
 *
 */
public interface Annotation extends Comparable<Annotation> {
	
	/**
	 * An enumeration of strand possibilities
	 * @author mguttman
	 *
	 */
	public enum Strand {
		POSITIVE('+'), NEGATIVE('-'), UNKNOWN('?'), BOTH('*'), INVALID('^');
		private char value;
		
		private Strand(char value) {
			this.value = value;
		}
		
		private Strand(String value) {
			if (value.length() > 1) throw new IllegalArgumentException("Illegal strand string");
			this.value = value.charAt(0);
		}

		public String toString() {
			return "" + value;
		}
		
		public Strand getReverseStrand() {
			if(this.equals(POSITIVE)){
				return NEGATIVE;
			}
			else if(this.equals(NEGATIVE)){
				return POSITIVE;
			}
			else if(this.equals(BOTH)){
				return BOTH;
			}
			else{
				return UNKNOWN;
			}
		}
		
		public static Strand consensusStrand(Strand strand1, Strand strand2) {
			Strand rtrn=INVALID;
			if(strand1.equals(UNKNOWN) || strand2.equals(UNKNOWN)){
				rtrn=UNKNOWN;
			}
			if(strand1.equals(BOTH)){
				rtrn=strand2;
			}
			if(strand2.equals(BOTH)){
				rtrn=strand1;
			}
			if(strand1.equals(strand2)){
				rtrn=strand1;
			}	
			return rtrn;
		}
		
		public static Strand fromString(String value) {
			if (value.equals("+")) return POSITIVE;
			if (value.equals("-")) return NEGATIVE;
			if(value.equals("*")) return BOTH;
			else return UNKNOWN;
		}

		public static Strand valueOf(boolean isNegativeStrandFlag) {
			if(isNegativeStrandFlag){return NEGATIVE;}
			else{return POSITIVE;}
		}
	}
	
	/**
	 * @return The name of this annotation
	 */
	public String getName();
	
	/**
	 * @return The name of the reference
	 */
	public String getReferenceName();
	
	/**
	 * @return The start position on the reference
	 */
	public int getReferenceStartPosition();
	
	/**
	 * @return The end position on the reference
	 */
	public int getReferenceEndPosition();
	
	/**
	 * @param other Another annotation
	 * @return True if any of the block in this annotation overlaps the blocks of the other annotation
	 */
	public boolean overlaps(Annotation other);

	/**
	 * @param other Another annotation
	 * @return True iff this annotation contains the other annotation
	 */
	public boolean contains(Annotation other);
	
	/**
	 * Merge with another annotation
	 * @param other Another annotation
	 * @return A new annotation representing a merging of the two
	 */
	public Annotation merge(Annotation other);
	
	/**
	 * Subtract another annotation from this annotation
	 * @param other Another annotation
	 * @return A new annotation representing the part of this annotation remaining after removing the other annotation
	 */
	public Annotation minus(Annotation other);
	
	/**
	 * Intersect another annotation with this annotation
	 * @param other Another annotation
	 * @return A new annotation representing the overlapping regions of the 2 annotations
	 */
	public Annotation intersect(Annotation other);
	
	/**
	 * Get blocks in the alignment
	 * @return The blocks of the alignment
	 */
	public Iterator<SingleInterval> getBlocks();
	
	
	
	/**
	 * Get blocks as a collection
	 * @return The set of blocks
	 */
	public Collection<Annotation> getBlockSet();
	
	/**
	 * Get number of blocks in the annotation
	 * @return The number of blocks
	 */
	public int getNumberOfBlocks();
	
	/**
	 * Return the size of the annotation
	 * @return size of the blocks in the annotation
	 */
	public int size();
	
	/**
	 * Return the orientation of the annotation
	 * @return
	 */
	public Strand getOrientation();
	
	/**
	 * Return a BED string representation of the Annotation
	 * @return String representation
	 */
	public String toBED();

	/**
	 * Return a BED string representation of the Annotation
	 * @param score BED score
	 * @return String representation
	 */
	public String toBED(double score);

	/**
	 * Return a BED string representation of the Annotation
	 * @param r Red color
	 * @param g Green color
	 * @param b Blue color
	 * @return String representation
	 */
	public String toBED(int r, int g, int b);
	
	/**
	 * Return a BED string representation of the Annotation
	 * @param r Red color
	 * @param g Green color
	 * @param b Blue color
	 * @param score BED score
	 * @return String representation
	 */
	public String toBED(int r, int g, int b, double score);
	
	/**
	 * Return a UCSC format representation of the interval covered by the annotation
	 * @return UCSC string
	 */
	public String toUCSC();
	
	/**
	 * Convert the featureAnnotation from feature space into reference space
	 * @param featureAnnotation region to convert in feature space
	 * @return Region in reference space
	 */
	public Annotation convertToReferenceSpace(Annotation featureAnnotation);
	
	/**
	 * Convert this annotation to the feature space represented by the feature passed
	 * @param feature The feature representing the new space
	 * @return this annotation in the new space
	 */
	public Annotation convert(Annotation feature);
	
	/**
	 * Get the relative coordinate in feature space relative to the 5' of the feature
	 * @param referenceStart Reference position
	 * @return Feature position (from 5' start) or -1 if doesn't overlap
	 */
	public int getRelativePositionFrom5PrimeOfFeature(int referenceStart);
	
	/**
	 * Return the CIGAR string representation of the annotation
	 * @return CIGAR format
	 */
	public String getCigarString();

	/**
	 * Get sliding windows over the annotation
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @return Collection of derived windows with pointers to this annotation
	 */
	public AnnotationCollection<DerivedAnnotation<? extends Annotation>> getWindows(int windowSize, int stepSize);
	
	/**
	 * Set the orientation of this feature
	 * @param orientation
	 */
	public void setOrientation(Strand orientation);

	/**
	 * Returns this Annotation as a SAMRecord
	 * @return A SAMRecord representation of this annotation
	 */
	public SAMRecord getSamRecord(SAMFileHeader header);
	
	/**
	 * Tests whether the object is fully contained within this one
	 * @param other The object to compare to this
	 * @return Whether the annotation is fully contained within this object
	 */
	public boolean fullyContained(Annotation other);
	
	/**
	 * Returns an annotation trimmed to the specified start and end coordinates
	 * @param start the lower bound of the trimmed annotation
	 * @param end the upper bound of the trimmed annotation
	 * @return
	 */
	public Annotation trim(int start, int end);
	
	/**
	 * 
	 * @param other annotation to be compared with
	 * @return -1,1,or 0 based on natural ordering defined in Abstract Annotation
	 */
	int compareTo(Annotation other);
	
}

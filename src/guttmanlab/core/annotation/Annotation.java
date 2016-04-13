package guttmanlab.core.annotation;

import guttmanlab.core.annotationcollection.AnnotationCollection;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.function.BiFunction;

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
	public default boolean overlaps(Annotation other) {
		if (other == null) {
			return false;
		}
		Iterator<SingleInterval> blocks1 = getBlocks();
		while (blocks1.hasNext()) {
			SingleInterval block1 = blocks1.next();
			Iterator<SingleInterval> blocks2 = other.getBlocks();
			while (blocks2.hasNext()){
				SingleInterval block2 = blocks2.next();
				if (block1.overlaps(block2)) {
					return true;
				}
			}	
		}
		return false;	}

	/**
	 * @param other Another annotation
	 * @return True iff this annotation contains the other annotation
	 */
	public default boolean contains(Annotation other) {
		return equals(union(other));
	}
	
	/**
	 * Merge with another annotation
	 * @param other Another annotation
	 * @return A new annotation representing a merging of the two
	 */
	public default Annotation merge(Annotation other) {
		if (other == null) {
			return this;
		}
		Strand consensusStrand = Strand.consensusStrand(getOrientation(), other.getOrientation());
		if (consensusStrand.equals(Strand.INVALID)) {
			return null;
		}
		BlockedAnnotation rtrn = new BlockedAnnotation();
		Iterator<SingleInterval> thisBlocks = getBlocks();
		while (thisBlocks.hasNext()) {
			rtrn.addBlocks(thisBlocks.next());
		}
		Iterator<SingleInterval> otherBlocks = other.getBlocks();
		while (otherBlocks.hasNext()) {
			rtrn.addBlocks(otherBlocks.next());
		}
		return rtrn;
	}
	
	/**
	 * Subtract another annotation from this annotation
	 * @param other Another annotation
	 * @return A new annotation representing the part of this annotation remaining after removing the other annotation
	 */
	public default Annotation minus(Annotation other) {
		return merge(other, (a, b) -> a && !b);
	}
	
	/**
	 * Intersect another annotation with this annotation
	 * @param other Another annotation
	 * @return A new annotation representing the overlapping regions of the 2 annotations
	 */
	public default Annotation intersect(Annotation other) {
		return merge(other, (a, b) -> a && b);
	}
	
	/**
	 * Creates the union of this annotation with another.
	 * @param other is the other annotation
	 * @return a new annotation which consists of all regions contained by either of the two input annotations
	 */
	public default Annotation union(Annotation other) {
		return merge(other, (a, b) -> a || b);
	}
	
	/**
	 * Creates the symmetric difference between this annotation and another.
	 * @param other is the other annotation
	 * @return a new annotation which consists of all regions contained by either of the two input annotations, but not both
	 */
	public default Annotation xor(Annotation other) {
		return merge(other, (a, b) -> a ^ b);
	}
	
	/**
	 * Get blocks in the alignment
	 * @return The blocks of the alignment
	 */
	public Iterator<SingleInterval> getBlocks();
	
	
	
	/**
	 * Get blocks as a collection
	 * @return The set of blocks
	 */
	public default Collection<Annotation> getBlockSet() {
		Iterator<SingleInterval> iter = getBlocks();
		Collection<Annotation> rtrn = new ArrayList<Annotation>();
		while(iter.hasNext()) {
			rtrn.add(iter.next());
		}
		return rtrn;
	}
	
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
	public default String toBED() {
		return toBED(0,0,0);
	}

	/**
	 * Return a BED string representation of the Annotation
	 * @param score BED score
	 * @return String representation
	 */
	public default String toBED(double score) {
		return toBED(0,0,0,score);
	}

	/**
	 * Return a BED string representation of the Annotation
	 * @param r Red color
	 * @param g Green color
	 * @param b Blue color
	 * @return String representation
	 */
	public default String toBED(int r, int g, int b) {
		return toBED(r, g, b, 0.0);
	}
	
	/**
	 * Return a BED string representation of the Annotation
	 * @param r Red color
	 * @param g Green color
	 * @param b Blue color
	 * @param score BED score
	 * @return String representation
	 */
	public default String toBED(int r, int g, int b, double score) {
		if(r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255) {
			throw new IllegalArgumentException("RGB values must be between 0 and 255");
		}
		String rgb = r + "," + g + "," + b;
		Iterator<SingleInterval> exons = getBlocks();
		String rtrn=getReferenceName()+"\t"+getReferenceStartPosition()+"\t"+getReferenceEndPosition()+"\t"+(getName() == null ? toUCSC() : getName())+"\t" + score + "\t"+getOrientation()+"\t"+getReferenceEndPosition()+"\t"+getReferenceEndPosition()+"\t"+rgb+"\t"+getNumberOfBlocks();
		String sizes="";
		String starts="";
		while(exons.hasNext()){
			SingleInterval exon=exons.next();
			sizes=sizes+(exon.size())+",";
			starts=starts+(exon.getReferenceStartPosition()-getReferenceStartPosition())+",";
		}
		rtrn=rtrn+"\t"+sizes+"\t"+starts;
		return rtrn;
	}
	
	/**
	 * Return a UCSC format representation of the interval covered by the annotation
	 * @return UCSC string
	 */
	public default String toUCSC() {
		return getReferenceName()+":"+getReferenceStartPosition()+"-"+getReferenceEndPosition();
	}
	
	
	public default Annotation convertToFeatureSpace(Annotation region) {
		//Ensure that region overlaps feature
		if(overlaps(region)){
			int featureStart=getRelativePositionFrom5PrimeOfFeature(region.getReferenceStartPosition());
			int featureEnd=getRelativePositionFrom5PrimeOfFeature(region.getReferenceEndPosition());
			Annotation interval;
			if(featureStart>-1 && featureEnd>-1){
				if(getOrientation().equals(Strand.NEGATIVE)){
					interval=new SingleInterval(getName(), featureEnd, featureStart); //TODO Check strand orientation
				}
				else{interval=new SingleInterval(getName(), featureStart, featureEnd);}
				return interval;
			}
		}
		return null;
	}

	
	/**
	 * Convert the featureAnnotation from feature space into reference space
	 * @param featureAnnotation region to convert in feature space
	 * @return Region in reference space
	 */
	public default Annotation convertToReferenceSpace(Annotation featureAnnotation) {
		BlockedAnnotation rtrn=new BlockedAnnotation();
		Iterator<SingleInterval> blocks = getBlocks();
		int sumBlocks=0;
		
		while(blocks.hasNext()){
			SingleInterval block=blocks.next();
			int origBlockSize = block.size();
			SingleInterval featureSpaceBlock=new SingleInterval(getName(), sumBlocks, sumBlocks+block.size());

			if(getOrientation().equals(Strand.NEGATIVE))
			{
				featureSpaceBlock= new SingleInterval(getName(), size()-(sumBlocks+block.size()),size()-sumBlocks); //FIXME TEST ME
			}
						
			if(featureAnnotation.overlaps(featureSpaceBlock)){
				//trim it, add it
				int shiftStart=0;
				int shiftEnd=0;
				if(featureAnnotation.getReferenceStartPosition()> featureSpaceBlock.getReferenceStartPosition()){
					shiftStart=featureAnnotation.getReferenceStartPosition()-featureSpaceBlock.getReferenceStartPosition();
				}
				if(featureAnnotation.getReferenceEndPosition()<featureSpaceBlock.getReferenceEndPosition())	{
					shiftEnd=featureSpaceBlock.getReferenceEndPosition()-featureAnnotation.getReferenceEndPosition();
				}
				block=block.trim(shiftStart, featureSpaceBlock.size()-shiftEnd);
				
				rtrn.addBlocks(block);
			}
			sumBlocks=sumBlocks+origBlockSize;
		}
		return rtrn;
	}
	
	/**
	 * Convert this annotation to the feature space represented by the feature passed
	 * @param feature The feature representing the new space
	 * @return this annotation in the new space
	 */
	public default Annotation convert(Annotation feature) {
		//Ensure that region overlaps feature
		if(overlaps(feature)){
			int featureStart=feature.getRelativePositionFrom5PrimeOfFeature(getReferenceStartPosition());
			int featureEnd=feature.getRelativePositionFrom5PrimeOfFeature(getReferenceEndPosition());
			Annotation interval;
			if(featureStart>-1 && featureEnd>-1){
				if(getOrientation().equals(Strand.NEGATIVE)){
					interval=new SingleInterval(feature.getName(), featureEnd, featureStart);
				}
				else{interval=new SingleInterval(feature.getName(), featureStart, featureEnd);}
				return interval;
			}
		}
		return null;
	}
	
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
	public default String getCigarString() {
		Iterator<SingleInterval> blocks=getBlocks();
		String cigar="";
		int distance = 0;
		
		int lastEnd=-1;
		while(blocks.hasNext()){
			SingleInterval block=blocks.next();
			if(lastEnd>0){
				distance=block.getReferenceStartPosition()-lastEnd;
				if (distance > 0)
					cigar+=distance+"N";
			}
			int m = block.size();
			if (distance < 0)
				m = m - distance;
			if ( m < 0 )
				System.err.println("inner read distance is negative.");
			cigar+=m+"M";
			lastEnd=block.getReferenceEndPosition();
		}
		return cigar;
	}

	/**
	 * Get sliding windows over the annotation
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @return Collection of derived windows with pointers to this annotation
	 */
	public AnnotationCollection<DerivedAnnotation<? extends Annotation>> getWindows(int windowSize, int stepSize);
	
	/**
	 * Sets the orientation of this feature and its constituent blocks.
	 * @param orientation is the desired orientation
	 */
	public void setOrientation(Strand orientation);

	/**
	 * Returns this Annotation as a SAMRecord
	 * @return A SAMRecord representation of this annotation
	 */
	public default SAMRecord getSamRecord(SAMFileHeader header) {
		SAMRecord record=new SAMRecord(header);
		record.setAlignmentStart(getReferenceStartPosition()+1);
		record.setCigarString(getCigarString());
		record.setReferenceName(getReferenceName());
		record.setReadName(getName());
		record.setReadNegativeStrandFlag(getOrientation().equals(Strand.NEGATIVE));
		return record;
	}
	
	/**
	 * Tests whether the object is fully contained within this one
	 * @param other The object to compare to this
	 * @return Whether the annotation is fully contained within this object
	 */
	public default boolean fullyContained(Annotation other) {
		//All blocks in other must be in blocks on this
		//Go through all blocks2 and check that they are in this
		Iterator<SingleInterval> blocks2=other.getBlocks();
		while(blocks2.hasNext()){
			SingleInterval block2=blocks2.next();
			boolean isInThis=false;
			//check that each block2 has some overlap with a block1
			Iterator<SingleInterval> blocks1=getBlocks();
			while(blocks1.hasNext()){
				SingleInterval block1=blocks1.next();
				if(block1.overlaps(block2)){
					isInThis=true;
					if(!AnnotationHelper.fullyContained(block1, block2)){
						return false;
					}
				}
			}
			if(!isInThis){return false;} //There are no blocks in this that the block overlapped, cant be fully contained
		}
		
		return true;
	}
	
	/**
	 * Returns an annotation trimmed to the specified start and end coordinates
	 * @param start the lower bound of the trimmed annotation
	 * @param end the upper bound of the trimmed annotation
	 * @return
	 */
	public default Annotation trim(int start, int end) {
		SingleInterval bound = new SingleInterval(this.getReferenceName(),start,end,this.getOrientation());
		Annotation a = this.intersect(bound);
		return a;
	}
	
	/**
	 * 
	 * @param other annotation to be compared with
	 * @return 
	 */
	public default int compareTo(Annotation other) {
		int comp = getReferenceName().compareTo(other.getReferenceName());
		if(comp!=0){return comp;}
		
		//second sort by start coordinate
		comp=getReferenceStartPosition() - other.getReferenceStartPosition();
		if(comp!=0){return comp;}
		
		//third sort by end coordinate
		comp=getReferenceEndPosition()-other.getReferenceEndPosition();
		if(comp!=0){return comp;}
		
		//fourth sort by strand
		comp=getOrientation().compareTo(other.getOrientation());
		if(comp!=0){return comp;}
		
		//Fifth sort by number of blocks
		comp=getNumberOfBlocks()-other.getNumberOfBlocks();
		if(comp!=0){return comp;}
		
		//Sixth sort by position of the blocks (in order scan)
		if(other.getNumberOfBlocks()>1){
			Iterator<SingleInterval> blocks = getBlocks();
			Iterator<SingleInterval> b_blocks = other.getBlocks();
			while(blocks.hasNext()){ //must have same number of blocks
				Annotation ann1=blocks.next();
				Annotation ann2=b_blocks.next();
				comp=ann1.compareTo(ann2);
				if(comp!=0){return comp;}
			}
		}
		return 0;
	}
	
	/**
	 * Converts this annotation into an array of integers corresponding to the
	 * start and end coordinates of the blocks. For example, if the annotation
	 * has two blocks, [3, 10) and [15, 20), the output would be the array
	 * [3, 10, 15, 20]. Other information, such as reference name and orientation,
	 * is lost.
	 * @return a flattened representation of this annotation
	 */
	public default int[] flatten() {
		int[] endpoints = new int[getNumberOfBlocks() * 2];
		int idx = 0;
		Iterator<SingleInterval> thisBlocks = getBlocks();
		while (thisBlocks.hasNext()) {
			SingleInterval block = thisBlocks.next();
			endpoints[idx++] = block.getReferenceStartPosition();
			endpoints[idx++] = block.getReferenceEndPosition();
		}
		return endpoints;
	}
	
	
	
	/**
	 * Merges this annotation with another. The type of merging is determined by the
	 * second parameter, op.
	 * @param other is the other annotation to merge with this one
	 * @param op is the function which defines the type of merge. For example, merging
	 * two annotations to obtain the union will require the argument to be
	 * (a, b) -> a || b
	 * @return the annotation resulting from the merge
	 */
	public default Annotation merge(Annotation other, BiFunction<Boolean, Boolean, Boolean> op) {
		if (other == null) {
			return null;
		}
		Strand consensus = Strand.consensusStrand(this.getOrientation(), other.getOrientation());
		if (consensus.equals(Strand.INVALID)) {
			return null;
		}
		if (!getReferenceName().equals(other.getReferenceName())) {
			return null;
		}
		if (getNumberOfBlocks() == 0 || other.getNumberOfBlocks() == 0) {
			return null;
		}
		
		// Flatten the annotations and add a sentinel value at the end
		int[] thisFlattened = flatten();
		int[] otherFlattened = other.flatten();
		
		int[] thisEndpoints = new int[thisFlattened.length + 1];
		for (int i = 0; i < thisFlattened.length; i++) {
			thisEndpoints[i] = thisFlattened[i];
		}
		int[] otherEndpoints = new int[otherFlattened.length + 1];
		for (int i = 0; i < otherFlattened.length; i++) {
			otherEndpoints[i] = otherFlattened[i];
		}

		int sentinel = Math.max(thisEndpoints[thisEndpoints.length - 2],
								otherEndpoints[otherEndpoints.length - 2]) + 1;
		thisEndpoints[thisEndpoints.length - 1] = sentinel;
		otherEndpoints[otherEndpoints.length - 1] = sentinel;
		
		// Go through the flattened annotations and at each endpoint, determine whether
		// it is in the result
		int thisIdx = 0;
		int otherIdx = 0;
		List<Integer> rtrnEndpoints = new ArrayList<Integer>();
		int scan = Math.min(thisEndpoints[thisIdx], otherEndpoints[otherIdx]);
		while (scan < sentinel) {
			boolean in_this = !((scan < thisEndpoints[thisIdx]) ^ (thisIdx % 2 == 1));
			boolean in_other = !((scan < otherEndpoints[otherIdx]) ^ (otherIdx % 2 == 1));
			boolean in_result = op.apply(in_this, in_other);
			
			if (in_result ^ (rtrnEndpoints.size() % 2 == 1)) {
				rtrnEndpoints.add(scan);
			}
			if (scan == thisEndpoints[thisIdx]) {
				thisIdx++;
			}
			if (scan == otherEndpoints[otherIdx]) {
				otherIdx++;
			}
			scan = Math.min(thisEndpoints[thisIdx], otherEndpoints[otherIdx]);
		}
//		for (int i = 0; i < rtrnEndpoints.size(); i++) {
//			System.out.println(rtrnEndpoints.get(i));
//		}
		// Construct the resulting annotation
		BlockedAnnotation rtrn = new BlockedAnnotation(getReferenceName());
		for (int i = 0; i < rtrnEndpoints.size(); i += 2) {
			rtrn.addBlocks(new SingleInterval(getReferenceName(), rtrnEndpoints.get(i), rtrnEndpoints.get(i + 1)));
		}
		rtrn.setOrientation(consensus);
		
		return rtrn;
	}
	
	/**
	 * An equals method that ignores the feature name
	 * @param other Other annotation
	 * @return True if the annotations are equal except for the name
	 */
	public default boolean equalsIgnoreName(Annotation other) {
		return compareTo(other) == 0;
	}
	
}

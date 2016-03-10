package guttmanlab.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.function.BiFunction;

import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * An abstract class that implements many of the shared features of an annotation
 * @author mguttman
 *
 */
public abstract class AbstractAnnotation implements Annotation {
	
	private static Logger logger = Logger.getLogger(AbstractAnnotation.class.getName());
	
	/**
	 * Gets the blocks of this annotation as a collection
	 * @return The set of blocks
	 */
	public Collection<Annotation> getBlockSet() {
		Iterator<SingleInterval> iter = getBlocks();
		Collection<Annotation> rtrn = new ArrayList<Annotation>();
		while(iter.hasNext()) {
			rtrn.add(iter.next());
		}
		return rtrn;
	}

	////////////////////
	// Set operations //
	////////////////////
	
	/**
	 * Converts this annotation into an array of integers corresponding to the
	 * start and end coordinates of the blocks. For example, if the annotation
	 * has two blocks, [3, 10) and [15, 20), the output would be the array
	 * [3, 10, 15, 20]. Other information, such as reference name and orientation,
	 * is lost.
	 * @returns a flattened representation of this annotation
	 */
	public int[] flatten() {
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

	@Override
	public Annotation merge(Annotation other) {
		Strand consensusStrand = Annotation.Strand.consensusStrand(getOrientation(), other.getOrientation());
		if(consensusStrand.equals(Strand.INVALID)) {
			return null;
		}
		BlockedAnnotation rtrn=new BlockedAnnotation();
		Iterator<SingleInterval> thisBlocks=getBlocks();
		while(thisBlocks.hasNext()) {
			rtrn.addBlocks(thisBlocks.next());
		}
		Iterator<SingleInterval> otherBlocks=other.getBlocks();
		while(otherBlocks.hasNext()) {
			rtrn.addBlocks(otherBlocks.next());
		}
		return rtrn;
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
	public Annotation merge(Annotation other, BiFunction<Boolean, Boolean, Boolean> op) {
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
		for (int i = 0; i < rtrnEndpoints.size(); i++) {
			System.out.println(rtrnEndpoints.get(i));
		}
		// Construct the resulting annotation
		BlockedAnnotation rtrn = new BlockedAnnotation(getReferenceName());
		for (int i = 0; i < rtrnEndpoints.size(); i += 2) {
			rtrn.addBlocks(new SingleInterval(getReferenceName(), rtrnEndpoints.get(i), rtrnEndpoints.get(i + 1)));
		}
		rtrn.setOrientation(consensus);
		
		return rtrn;
	}
	
	@Override
	public Annotation minus(Annotation other) {
		return merge(other, (a, b) -> a && !b);
	}
	
	@Override
	public Annotation union(Annotation other) {
		return merge(other, (a, b) -> a || b);
	}
	
	@Override
	public Annotation intersect(Annotation other) {
		return merge(other, (a, b) -> a && b);
	}
	
	@Override
	public Annotation xor(Annotation other) {
		return merge(other, (a, b) -> a ^ b);
	}
	
	@Override
	public Annotation merge(Annotation other) {
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
	
	@Override
	public boolean overlaps(Annotation other) {
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
		return false;
	}

	@Override
	public boolean contains(Annotation other) {
		return equals(union(other));
	}
	
	protected SingleInterval merge(SingleInterval block1, SingleInterval block2) {
		if(!overlaps(block1, block2)){return null;}
		
		int newStart=Math.min(block1.getReferenceStartPosition(), block2.getReferenceStartPosition());
		int newEnd=Math.max(block1.getReferenceEndPosition(), block2.getReferenceEndPosition());
		Strand consensus=Annotation.Strand.consensusStrand(block1.getOrientation(), block2.getOrientation());
		return new SingleInterval(block1.getReferenceName(), newStart, newEnd, consensus);
	}
	
	/**
	 * Helper method to calculate overlaps from single blocks
	 * @param block1
	 * @param block2
	 * @return whether the blocks overlap
	 */
	private boolean overlaps(SingleInterval block1, SingleInterval block2){
		int newStart=Math.max(block1.getReferenceStartPosition(), block2.getReferenceStartPosition());
		int newEnd=Math.min(block1.getReferenceEndPosition(), block2.getReferenceEndPosition());
		
		Strand consensusStrand=Annotation.Strand.consensusStrand(block1.getOrientation(), block2.getOrientation());
		if(newStart<newEnd && block1.getReferenceName().equalsIgnoreCase(block2.getReferenceName()) && !consensusStrand.equals(Annotation.Strand.INVALID)){
			return true;
		}
		
		return false;
	}
	
	@Override
	public String toString(){
		return toBED(0,0,0);
	}
	
	public String toBED() {
		return toBED(0,0,0);
	}
	
	public String toBED(double score) {
		return toBED(0,0,0,score);
	}
	
	public String toBED(int r, int g, int b){
		return toBED(r, g, b, 0.0);
	}
	
	public String toBED(int r, int g, int b, double score){
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

	public String toUCSC() {
		return getReferenceName()+":"+getReferenceStartPosition()+"-"+getReferenceEndPosition();
	}
	
	public Annotation convertToFeatureSpace(Annotation region){
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
	
	public Annotation convert(Annotation feature){
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
	
	public Annotation convertToReferenceSpace(Annotation featureAnnotation){
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
	
	@Override
	public String getCigarString(){
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
	
	@Override
	public SAMRecord getSamRecord(SAMFileHeader header){
		SAMRecord record=new SAMRecord(header);
		record.setAlignmentStart(getReferenceStartPosition()+1);
		record.setCigarString(getCigarString());
		record.setReferenceName(getReferenceName());
		record.setReadName(getName());
		record.setReadNegativeStrandFlag(getOrientation().equals(Strand.NEGATIVE));
		return record;
	}
	
	public boolean fullyContained(Annotation other){
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
					if(!fullyContained(block1, block2)){
						return false;
					}
				}
			}
			if(!isInThis){return false;} //There are no blocks in this that the block overlapped, cant be fully contained
		}
		
		return true;
	}
	
	private boolean fullyContained(SingleInterval block1, SingleInterval block2){
		//is only fully contained if block2 is equal to or a subset of block1
		if(block1.getReferenceStartPosition()<=block2.getReferenceStartPosition() && block1.getReferenceEndPosition()>=block2.getReferenceEndPosition()){
			return true;
		}
		return false;	
	}
	
	public Annotation trim(int start,int end)
	{
		SingleInterval bound = new SingleInterval(this.getReferenceName(),start,end,this.getOrientation());
		Annotation a = this.intersect(bound);
		return a;
	}
	
	@Override
	public int compareTo (Annotation other) {
		return compareToAnnotation(other);
	}
	
	
	public int compareToAnnotation(Annotation b) {
		return compareToAnnotation(b, true);
	}
	
	public int compareToAnnotation(Annotation b, boolean useOrientation) {
		int comp = getReferenceName().compareTo(b.getReferenceName());
		if(comp!=0){return comp;}
		
		//second sort by start coordinate
		comp=getReferenceStartPosition() - b.getReferenceStartPosition();
		if(comp!=0){return comp;}
		
		//third sort by end coordinate
		comp=getReferenceEndPosition()-b.getReferenceEndPosition();
		if(comp!=0){return comp;}
		
		//fourth sort by strand
		if (useOrientation) {
			comp=getOrientation().compareTo(b.getOrientation());
			if(comp!=0){return comp;}
		}
		
		//Fifth sort by number of blocks
		comp=getNumberOfBlocks()-b.getNumberOfBlocks();
		if(comp!=0){return comp;}
		
		//Sixth sort by position of the blocks (in order scan)
		if(b.getNumberOfBlocks()>1){
			Iterator<SingleInterval> blocks = getBlocks();
			Iterator<SingleInterval> b_blocks = b.getBlocks();
			while(blocks.hasNext()){ //must have same number of blocks
				Annotation ann1=blocks.next();
				Annotation ann2=b_blocks.next();
				comp=ann1.compareTo(ann2);
				if(comp!=0){return comp;}
			}
		}
		return 0;
	}
	
	@Override
	public boolean equals(Object other)
	{
		if (other == null) {
			return false;
		}
		if (other instanceof Annotation) {
			return this.compareTo((Annotation)other) == 0;
		}
		return false;
	}
	
	protected HashCodeBuilder hashCodeBuilder() {
		return new HashCodeBuilder(31,37).append(getReferenceName()).append(getReferenceStartPosition()).append(getReferenceEndPosition()).append(getOrientation()).append(getNumberOfBlocks());
	}
	
	@Override
	public int hashCode()
	{
		return hashCodeBuilder().toHashCode();
	}
}

package guttmanlab.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import org.apache.commons.lang.builder.HashCodeBuilder;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * An abstract class that implements many of the shared features of an annotation
 * @author mguttman
 *
 */
public abstract class AbstractAnnotation implements Annotation {
	
	@Override
	public Annotation intersect(Annotation other) {
		BlockedAnnotation rtrn=new BlockedAnnotation();
		Iterator<SingleInterval> blocks1=getBlocks();
		while(blocks1.hasNext()){
			SingleInterval block1=blocks1.next();
			Iterator<SingleInterval> blocks2=other.getBlocks();
			while(blocks2.hasNext()){
				SingleInterval block2=blocks2.next();
				SingleInterval inter=intersect(block1, block2);
				if(inter!=null){rtrn.addBlocks(inter);}
			}
			
		}
		return rtrn;
	}
	
	/**
	 * Get blocks as a collection
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
	

	/**
	 * Helper method to compute the overlap between single blocks
	 * @param block1 Block1
	 * @param block2 Block2
	 * @return The intersection or null if no intersection exists
	 */
	private SingleInterval intersect(SingleInterval block1, SingleInterval block2) {
		if(!overlaps(block1, block2)){return null;}
		int newStart=Math.max(block1.getReferenceStartPosition(), block2.getReferenceStartPosition());
		int newEnd=Math.min(block1.getReferenceEndPosition(), block2.getReferenceEndPosition());
		Strand consensus=Annotation.Strand.consensusStrand(block1.getOrientation(), block2.getOrientation());
		return new SingleInterval(block1.getReferenceName(), newStart, newEnd, consensus);
	}

	@Override
	public Annotation merge(Annotation other) {
		throw new UnsupportedOperationException("Method broken");
//		BlockedAnnotation rtrn=new BlockedAnnotation();
//		Iterator<SingleInterval> blocks1=getBlocks();
//		while(blocks1.hasNext()){
//			SingleInterval block1=blocks1.next();
//			Iterator<SingleInterval> blocks2=other.getBlocks();
//			while(blocks2.hasNext()){
//				SingleInterval block2=blocks2.next();
//				if(block1.overlaps(block2)){
//					SingleInterval merge=merge(block1, block2);
//					if(merge!=null){rtrn.addBlocks(merge);}
//				}
//			}
//			
//		}
//		return rtrn;
	}

	protected SingleInterval merge(SingleInterval block1, SingleInterval block2) {
		if(!overlaps(block1, block2)){return null;}
		
		int newStart=Math.min(block1.getReferenceStartPosition(), block2.getReferenceStartPosition());
		int newEnd=Math.max(block1.getReferenceEndPosition(), block2.getReferenceEndPosition());
		Strand consensus=Annotation.Strand.consensusStrand(block1.getOrientation(), block2.getOrientation());
		return new SingleInterval(block1.getReferenceName(), newStart, newEnd, consensus);
	}

	@Override
	public Annotation minus(Annotation other) {
		// FIXME Auto-generated method stub
		throw new UnsupportedOperationException("TODO");
	}
	
	@Override
	public boolean overlaps(Annotation other) {
		//TODO This method still needs to be tested to ensure that it does what we expect
		//Check if the blocks overlap
		Iterator<SingleInterval> blocks1=getBlocks();
		while(blocks1.hasNext()){
			SingleInterval block1=blocks1.next();
			Iterator<SingleInterval> blocks2=other.getBlocks();
			while(blocks2.hasNext()){
				SingleInterval block2=blocks2.next();
				if(overlaps(block1, block2)){
					return true;
				}
			}	
		}
		return false;
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
	public boolean contains(Annotation other) {
		// FIXME Auto-generated method stub
		throw new UnsupportedOperationException("TODO");
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
		if(other instanceof Annotation)
			return this.compareToAnnotation((Annotation)other)==0;
		return false;
	}
	
	@Override
	public int hashCode()
	{
		return new HashCodeBuilder(31,37).append(getReferenceName()).append(getReferenceStartPosition()).append(getReferenceEndPosition()).append(getOrientation()).append(getNumberOfBlocks()).toHashCode();
	}
}

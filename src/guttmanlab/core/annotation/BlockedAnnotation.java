package guttmanlab.core.annotation;

import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.datastructures.IntervalTree;

import java.util.Collection;
import java.util.Iterator;



public class BlockedAnnotation extends AbstractAnnotation {

	private IntervalTree<SingleInterval> blocks;
	private String referenceName;
	private int startPosition;
	private int endPosition;
	private int size;
	private boolean started;
	private String name;
	private Strand orientation;
	
	/**
	 * An empty constructor
	 */
	public BlockedAnnotation(){
		this.blocks=new IntervalTree<SingleInterval>();
		this.started=false;
	}
	
	/**
	 * Initialize as empty but with name
	 */
	public BlockedAnnotation(String name){
		this();
		this.name=name;
	}
	
	/**
	 * A blocked annotation is defined by its blocks
	 * @param blocks
	 * @param name 
	 */
	public BlockedAnnotation(Collection<Annotation> blocks, String name){
		this(name);
		for(Annotation block: blocks){
			addBlocks(block);
		}
	}
	
	public BlockedAnnotation(Annotation annot){
		this();
		this.name=annot.getName();
		addBlocks(annot);
	}
	
	public BlockedAnnotation(Annotation annot, String newName){
		this(annot);
		this.name = newName;
	}
	
	/**
	 * Add block to current blocks
	 * If overlaps existing block, merge
	 * Requires blocks have the same reference name and Strand
	 * @param annot Annotation with blocks to add
	 * @return if the block was successfully added
	 */
	public boolean addBlocks(Annotation annot) {
		boolean added=false;
		Iterator<SingleInterval> exons=annot.getBlocks();
		while(exons.hasNext()){
			added &= update(exons.next());
		}
		return added;
	}

	/**
	 * Helper method to add single interval
	 * @param interval The interval
	 * @return whether it was successfully added
	 */
	private boolean update(SingleInterval interval) {
		if(!started){
			this.referenceName=interval.getReferenceName();
			this.startPosition=interval.getReferenceStartPosition();
			this.endPosition=interval.getReferenceEndPosition();
			this.orientation=interval.getOrientation();
			started=true;
		}
		else{
			if(!this.referenceName.equalsIgnoreCase(interval.getReferenceName())){return false;}
			if(!this.orientation.equals(interval.getOrientation())){return false;}
			this.startPosition=Math.min(startPosition, interval.getReferenceStartPosition());
			this.endPosition=Math.max(endPosition, interval.getReferenceEndPosition());
		}
		
		boolean hasOverlappers=blocks.hasOverlappers(interval.getReferenceStartPosition(), interval.getReferenceEndPosition());
		SingleInterval merged=interval;
		if(hasOverlappers){
			//pull, merge, and update
			Iterator<SingleInterval> iter=blocks.overlappingValueIterator(interval.getReferenceStartPosition(), interval.getReferenceEndPosition());
			while(iter.hasNext()){
				SingleInterval e=iter.next();
				blocks.remove(e.getReferenceStartPosition(), e.getReferenceEndPosition()); // Pam added on 12/24/14
				merged=merge(merged, e);
				size-=e.size();
			}
		}
		int end = merged.getReferenceEndPosition();
		int start = merged.getReferenceStartPosition();
		if(start>end) {
			blocks.put(end,start, merged);
		}
		else {
			blocks.put(start, end, merged);
		}
		size+=merged.size();
		
		return true;
	}

	@Override
	public String getName() {
		return this.name;
	}

	public Iterator<SingleInterval> getBlocks() {
		return this.blocks.valueIterator();
	}

	@Override
	public String getReferenceName() {
		return this.referenceName;
	}

	@Override
	public int getReferenceStartPosition() {
		return this.startPosition;
	}

	@Override
	public int getReferenceEndPosition() {
		return this.endPosition;
	}

	@Override
	public int size() {
		return this.size;
	}

	@Override
	public Strand getOrientation() {
		return this.orientation;
	}


	@Override
	public int getNumberOfBlocks() {
		return blocks.size();
	}
	
	//TODO This could actually go in the AbstractAnnotation
	public int getRelativePositionFrom5PrimeOfFeature(int referenceStart){
		if(referenceStart>=this.getReferenceEndPosition() || referenceStart<this.getReferenceStartPosition()){return -1;} //This start position is past the feature
		Iterator<SingleInterval> iter=this.blocks.overlappingValueIterator(this.getReferenceStartPosition(), referenceStart);
		int relativeSize=0;
		while(iter.hasNext()){
			SingleInterval interval=iter.next();
			if(interval.getReferenceEndPosition()<referenceStart){
				relativeSize+=interval.size(); //except when overlapping exactly the referenceStart
			}
			else{
				relativeSize+=(referenceStart-interval.getReferenceStartPosition());
			}
		}
		
		//If strand is neg, then position is from end
		if(getOrientation().equals(Annotation.Strand.NEGATIVE)){
			relativeSize=this.size-relativeSize - 1;
		}
		return relativeSize;
	}
	
	public BlockedAnnotation convertToFeatureSpace(SingleInterval region){
		int featureStart=getRelativePositionFrom5PrimeOfFeature(region.getReferenceStartPosition());
		int featureEnd=getRelativePositionFrom5PrimeOfFeature(region.getReferenceEndPosition());
		BlockedAnnotation interval;
		if(getOrientation().equals(Strand.NEGATIVE)){
			interval=new BlockedAnnotation(new SingleInterval(getName(), featureEnd, featureStart));
		}
		else{interval=new BlockedAnnotation(new SingleInterval(getName(), featureStart, featureEnd));}
		return interval;
	}

	@Override
	public void setOrientation(Strand orientation) {
		this.orientation=orientation;
	}

	@Override
	public AnnotationCollection<DerivedAnnotation<? extends Annotation>> getWindows(int windowSize, int stepSize) {
		
		FeatureCollection<DerivedAnnotation<? extends Annotation>> rtrn = new FeatureCollection<DerivedAnnotation<? extends Annotation>>(null);
		
		boolean plusStrand = getOrientation().equals(Strand.NEGATIVE) ? false : true;
		int featureStart = plusStrand ? 0 : size - windowSize;
		int featureEnd = featureStart + windowSize;
		
		while(featureEnd <= size && featureStart >= 0) {
			Annotation windowFeatureSpace = new SingleInterval(getName(), featureStart, featureEnd);
			Annotation windowReferenceSpace = convertToReferenceSpace(windowFeatureSpace);
			DerivedAnnotation<BlockedAnnotation> windowDerived = new DerivedAnnotation<BlockedAnnotation>(windowReferenceSpace, 
					windowReferenceSpace.toUCSC(), this);
			rtrn.add(windowDerived);
			if(plusStrand) {
				featureStart += stepSize;
				featureEnd += stepSize;
			} else {
				featureStart -= stepSize;
				featureEnd -= stepSize;
			}
		}
		
		return rtrn;
		
	}
}

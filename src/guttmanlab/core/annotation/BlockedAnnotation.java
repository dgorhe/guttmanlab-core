package guttmanlab.core.annotation;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;



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
	public BlockedAnnotation(Collection<? extends Annotation> blocks, String name){
		this(name);
		for(Annotation block: blocks){
			addBlocks(block);
		}
	}
	
	public BlockedAnnotation(Collection<? extends Annotation> blocks){
		this();
		for(Annotation block: blocks){
			addBlocks(block);
		}
	}
	
	public BlockedAnnotation(Annotation annot){
		this();
		this.orientation=annot.getOrientation();
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
	
	
	private IntervalTree<Annotation> getGeneTree() {
		boolean forward=getOrientation().equals(Strand.POSITIVE);
		Iterator <SingleInterval> exons=getSortedExons(forward);
		
		IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
		
		int currentStartPosition=0;
		int currentEndPosition=0;
		while(exons.hasNext()){
			Annotation exon=exons.next();
			currentEndPosition+=exon.size();
			tree.put(currentStartPosition,  currentEndPosition, exon);
			currentStartPosition=currentEndPosition;
		}
		return tree;
	}

	
	private Iterator<SingleInterval> getSortedExons(boolean forward) {
		return this.sortedBlockIterator(forward);
	}

	public Annotation trimByRelativePositions(int relativeStart, int relativeEnd) {
		return trim(relativeStart, relativeEnd, getName(), getGeneTree());
	}
	
	private Annotation trim(int start, int end, String name, IntervalTree<Annotation> tree) {
		
		Collection<Annotation> trimmedBlocks=new TreeSet<Annotation>();
		Iterator<Node<Annotation>> overlappers=tree.overlappers(start, end);
		while(overlappers.hasNext()){
			Node<Annotation> node=overlappers.next();
			Annotation exon=node.getValue();
			
			
			int trimmedStart=exon.getReferenceStartPosition();
			int trimmedEnd=exon.getReferenceEndPosition();
			
			if(node.getStart()< start){
				int beginDiff=start-node.getStart();
				if(getOrientation().equals(Strand.POSITIVE)){
					trimmedStart=exon.getReferenceStartPosition()+beginDiff; //TODO If reverse orientation then we need to take this diff from block end
				}
				else{
					trimmedEnd=exon.getReferenceEndPosition()-beginDiff;
				}
			}
			
			if(node.getEnd()>end){
				int endDiff=node.getEnd()-end;
				if(getOrientation().equals(Strand.POSITIVE)){
					trimmedEnd=exon.getReferenceEndPosition()-endDiff;
				}
				else{
					trimmedStart=exon.getReferenceStartPosition()+endDiff;
				}
			}
			
			Annotation trimmed=new SingleInterval(exon.getReferenceName(), trimmedStart, trimmedEnd, exon.getOrientation());
			trimmedBlocks.add(trimmed);
		}
		
		Annotation rtrn=new BlockedAnnotation(trimmedBlocks, name);
		rtrn.setOrientation(this.getOrientation());
		return rtrn;
	}
	
	public Iterator<SingleInterval> sortedBlockIterator(boolean forward){
		if(forward){return this.blocks.valueIterator();}
		return this.blocks.reverseValueIterator();
	}
	
	public SingleInterval getFirstBlock(){
		return this.blocks.min().getValue();
	}
	
	public SingleInterval getLastBlock(){
		return this.blocks.maxValue();
	}
	
	/**
	 * Return block less than or equal to position
	 * @param position
	 * @return
	 */
	public SingleInterval getMaxBlock(int position){
		Iterator<SingleInterval> iter=this.blocks.overlappingValueIterator(position+1, position+1);
		if(iter.hasNext()){return iter.next();}
		
		
		Node<SingleInterval> interval=this.blocks.max(position+1, position+1);
		if(interval!=null){return interval.getValue();}
		return null;
	}
	
	/**
	 * Return block greater than or equal to position
	 * @param position
	 * @return
	 */
	public SingleInterval getMinBlock(int position){
		Iterator<SingleInterval> iter=this.blocks.overlappingValueIterator(position-1, position);
		if(iter.hasNext()){return iter.next();}
		
		
		Node<SingleInterval> interval=this.blocks.min(position, position);
		if(interval!=null){return interval.getValue();}
		return null;
	}
	
	public boolean isWithinBlock(int position){
		Iterator<SingleInterval> iter=this.blocks.overlappingValueIterator(position-1, position+1);
		return iter.hasNext();
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
	
	public Iterator<SingleInterval> getOverlappingBlocks(Annotation region){
		return this.blocks.overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
	}
	
	public int getNumberOfOverlappingBlocks(Annotation region) {
		return this.blocks.numOverlappers(region.getReferenceStartPosition(), region.getReferenceEndPosition());
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
	public void setName(String name) {
		this.name=name;
	}

	
}

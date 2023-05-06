package guttmanlab.core.annotation;

import guttmanlab.core.annotationcollection.AnnotationCollection;
import htsjdk.variant.variantcontext.VariantContext;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

/**
 * A class that represent a simple contiguous interval
 * This is the basis for all features
 * @author mguttman
 *
 */
public class SingleInterval extends AbstractAnnotation{

	private String referenceName;
	private int startPos;
	private int endPos;
	private Strand orientation;
	private String featureName;
	private int count;
	
	public SingleInterval(String refName, int start, int end, Strand orientation){
		this(refName, start, end, orientation, "");
	}
	
	public SingleInterval(String refName, int start, int end, Strand orientation, String featureName){
		this.referenceName=refName;
		this.startPos=start;
		this.endPos=end;
		this.orientation=orientation;
		this.featureName=featureName;
		this.count=0;
	}
	
	public SingleInterval(String ucscString){
		if(ucscString.endsWith("+") || ucscString.endsWith("(+)")) {
			this.orientation=Strand.POSITIVE;
			if(ucscString.endsWith("(+)")) {ucscString=ucscString.substring(0,ucscString.length()-3);}
			if(ucscString.endsWith("+")) {ucscString=ucscString.substring(0,ucscString.length()-1);}
		}
		if(ucscString.endsWith("-")|| ucscString.endsWith("(-)")) {
			this.orientation=Strand.NEGATIVE;
			if(ucscString.endsWith("(-)")) {ucscString=ucscString.substring(0,ucscString.length()-3);}
			if(ucscString.endsWith("-")) {ucscString=ucscString.substring(0,ucscString.length()-1);}
		}
		this.referenceName=ucscString.split(":")[0];
		this.startPos=new Integer(ucscString.split(":")[1].split("-")[0]);
		this.endPos=new Integer(ucscString.split(":")[1].split("-")[1]);
		this.featureName=ucscString;
		//this.orientation=Strand.BOTH;
		this.count=0;
	}
	
	public int getCount(){return this.count;}
	
	public SingleInterval(String refName, int start, int end) {
		this(refName, start, end, Strand.BOTH, "");
	}
	
	public SingleInterval(VariantContext variant) {
		String chr=variant.getChr();
		int start=variant.getStart();
		int end=variant.getEnd();
		new SingleInterval(chr, start, end, Strand.POSITIVE);
	}

	public boolean isWithinBlock(int position){
		if(position>=getReferenceStartPosition() && position <= getReferenceEndPosition()){return true;}
		return false;
	}

	@Override
	public String getName() {
		return this.featureName;
	}
	
	public void setName(String name) {
		this.featureName=name;
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


	@Override
	public Iterator<SingleInterval> getBlocks() {
		Collection<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		rtrn.add(this);
		return rtrn.iterator();
	}
	
	@Override
	public int size() {
		return endPos-startPos;
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
	
	public SingleInterval getMidPoint() {
		int midPoint=this.getReferenceStartPosition()+(this.getLength()/2);
		return new SingleInterval(this.referenceName, midPoint, midPoint+1);
	}

	@Override
	public void setOrientation(Strand orientation) {
		this.orientation=orientation;
	}
	
	@Override
	public boolean contains(Annotation other) {
		if(other.getReferenceName().equalsIgnoreCase(this.getReferenceName()) && other.getReferenceStartPosition()>=this.getReferenceStartPosition() && other.getReferenceEndPosition()<=this.getReferenceEndPosition()){
			return true;
		}
		return false;
	}

	/**
	 * Trim this block to the relative start and end position provided
	 * @param relativeStartPosition relative start position
	 * @param relativeEndPosition relative end position
	 * @return
	 */
	public SingleInterval trim(int relativeStart, int relativeEnd) {
		if(getOrientation().equals(Strand.NEGATIVE)){
			int newEnd=getReferenceEndPosition()-relativeStart;
			int newStart=getReferenceEndPosition()-relativeEnd;
			return new SingleInterval(getReferenceName(), newStart, newEnd);
		}
		else{
			return new SingleInterval(getReferenceName(), getReferenceStartPosition()+relativeStart, getReferenceStartPosition()+relativeEnd);
		}
	}

	public String getFileName() {
		return getReferenceName()+"_"+getReferenceStartPosition()+"_"+getReferenceEndPosition();
	}

	public int getLength() {
		return getReferenceEndPosition()-getReferenceStartPosition();
	}

	@Override
	public Annotation trimByRelativePositions(int relativeStart, int relativeEnd) {
		SingleInterval rtrn;
		if(getOrientation().equals(Strand.NEGATIVE)){
			int newEnd=this.getReferenceEndPosition()-relativeStart;
			int newStart=this.getReferenceEndPosition()-relativeEnd;
			rtrn=new SingleInterval(this.getReferenceName(), newStart, newEnd);
		}
		else{
			rtrn= new SingleInterval(this.getReferenceName(), this.getReferenceStartPosition()+relativeStart, this.getReferenceStartPosition()+relativeEnd);
		}
		rtrn.setOrientation(this.getOrientation());
		return rtrn;
	}

	public String toNode() {
		return getReferenceName()+"_"+getReferenceStartPosition();
	}

	public void addCount() {
		count++;
		
	}

	public void setCount(int numberOfAnnotationsInWindow) {
		this.count=numberOfAnnotationsInWindow;
	}

	@Override
	public SingleInterval bin(int resolution) {
		int startIndex=getReferenceStartPosition()/resolution;
		int newStart=startIndex*resolution;
		int newEnd=newStart+Math.max(getLength(), resolution);
		SingleInterval newInterval=new SingleInterval(getReferenceName(), newStart, newEnd);
		return newInterval;
	}
	
	
	public Collection<SingleInterval> getBins(int resolution) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		int startIndex=getReferenceStartPosition()/resolution;
		int endIndex=getReferenceEndPosition()/resolution;
		
		if(startIndex==endIndex) {
			int newStart=startIndex*resolution;
			int newEnd=newStart+resolution;
			SingleInterval newInterval=new SingleInterval(getReferenceName(), newStart, newEnd);
			rtrn.add(newInterval);
		}
		
		else {
			for(int i=startIndex; i<=endIndex; i++) {
				int newStart=(i)*resolution;
				int newEnd=newStart+resolution;
				SingleInterval newInterval=new SingleInterval(getReferenceName(), newStart, newEnd);
				rtrn.add(newInterval);
			}
		}
		
		
		return rtrn;
	}
	

	public String toBedgraph(Double score) {
		return getReferenceName()+"\t"+getReferenceStartPosition()+"\t"+getReferenceEndPosition()+"\t"+score;
	}
	
	public String toBedgraph(int score) {
		return getReferenceName()+"\t"+getReferenceStartPosition()+"\t"+getReferenceEndPosition()+"\t"+score;
	}

	public Collection<SingleInterval> getWindowsCollection(int binResolution, int step) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		CloseableIterator<DerivedAnnotation<? extends Annotation>> iter= this.getWindows(binResolution, step).sortedIterator();
		while(iter.hasNext()) {
			rtrn.add(iter.next().getSingleInterval());
		}
		return rtrn;
	}

	public String toShortBED() {
		return this.getReferenceName()+"\t"+this.getReferenceStartPosition()+"\t"+this.getReferenceEndPosition();
	}

	public static SingleInterval bin(SAMRecord record, int binSize) {
		SAMFragment f=new SAMFragment(record);
		SingleInterval binned=f.getSingleInterval().bin(binSize);
		return binned;
	}

	public String toShortBED(String name) {
		return this.getReferenceName()+"\t"+this.getReferenceStartPosition()+"\t"+this.getReferenceEndPosition()+"\t"+name;
	}

	public Collection<SingleInterval> allBins(int resolution) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		int startIndex=getReferenceStartPosition()/resolution;
		
		
		int endIndex=getReferenceEndPosition()/resolution;
		
		for(int i=startIndex; i<=endIndex; i++) {
			int newStart=i*resolution;
			int newEnd=newStart+resolution;
			SingleInterval newInterval=new SingleInterval(getReferenceName(), newStart, newEnd);
			newInterval.setOrientation(this.orientation);
			rtrn.add(newInterval);
			
		}
		
		
		return rtrn;
		
	}
	
	
	public Collection<SingleInterval> allBins(int resolution, Strand strand) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		int startIndex=getReferenceStartPosition()/resolution;
		
		
		int endIndex=getReferenceEndPosition()/resolution;
		
		for(int i=startIndex; i<=endIndex; i++) {
			int newStart=i*resolution;
			int newEnd=newStart+resolution;
			SingleInterval newInterval=new SingleInterval(getReferenceName(), newStart, newEnd);
			newInterval.setOrientation(strand);
			rtrn.add(newInterval);
			
		}
		
		
		return rtrn;
		
	}
	
	public Collection<SingleInterval> allBins(int[] binSizes) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		for(int i=0; i<binSizes.length; i++) {
			rtrn.addAll(allBins(binSizes[i]));
		}
		return rtrn;
	}
	
	public Collection<SingleInterval> allBins(int resolution, int stagger) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		int startIndex=getReferenceStartPosition()/stagger;
		
		
		int endIndex=getReferenceEndPosition()/stagger;
		
		for(int i=startIndex; i<=endIndex; i++) {
			int newStart=i*resolution;
			int newEnd=newStart+resolution;
			SingleInterval newInterval=new SingleInterval(getReferenceName(), newStart, newEnd);
			rtrn.add(newInterval);
			
		}
		
		
		/*for(int i=getReferenceStartPosition()-resolution; i<getReferenceEndPosition(); i++) {
			int start=i;
			int end=i+resolution;
			SingleInterval newInterval=new SingleInterval(getReferenceName(), start, end);
			rtrn.add(newInterval);
		}*/
		
		return rtrn;
		
	}
	
	
	public Collection<SingleInterval> allBins(int resolution, int stagger, Strand strand) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		int startIndex=getReferenceStartPosition()/stagger;
		
		
		int endIndex=getReferenceEndPosition()/stagger;
		
		for(int i=startIndex; i<=endIndex; i++) {
			int newStart=i*resolution;
			int newEnd=newStart+resolution;
			SingleInterval newInterval=new SingleInterval(getReferenceName(), newStart, newEnd);
			newInterval.setOrientation(strand);
			rtrn.add(newInterval);
			
		}
		
		
		/*for(int i=getReferenceStartPosition()-resolution; i<getReferenceEndPosition(); i++) {
			int start=i;
			int end=i+resolution;
			SingleInterval newInterval=new SingleInterval(getReferenceName(), start, end);
			rtrn.add(newInterval);
		}*/
		
		return rtrn;
		
	}
	
	
	public static void main (String[] args) throws IOException {
		SingleInterval i=new SingleInterval(args[0]);
		Collection<SingleInterval> list=i.allBins(1000000);
		for(SingleInterval l: list) {System.err.println(l.toUCSC());}
	}

	public static SingleInterval antisense(SingleInterval intron) {
		SingleInterval rtrn=new SingleInterval(intron.getReferenceName(), intron.getReferenceStartPosition(), intron.getReferenceEndPosition());
		Strand orientation=intron.getOrientation();
		if(intron.getOrientation().equals(Strand.POSITIVE)) {orientation=Strand.NEGATIVE;}
		if(intron.getOrientation().equals(Strand.NEGATIVE)) {orientation=Strand.POSITIVE;}
		rtrn.setOrientation(orientation);
		return rtrn;
	}

	
	
	
}

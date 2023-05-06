package guttmanlab.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.sequence.Sequence;

/**
 * A simple extension of a BlockedAnnotation that keeps track of the start and end of a coding region
 * @author mguttman
 *
 */
public class Gene extends BlockedAnnotation{

	private static final int NO_CDS=-1;
	private int cdsStartPos=NO_CDS;
	private int cdsEndPos=NO_CDS;
	
	private IntervalTree<SingleInterval> intronTree;
	
	public Gene(Collection<Annotation> blocks, int cdsStartPos, int cdsEndPos, String name) {
		super(blocks, name);
		this.cdsEndPos=cdsEndPos;
		this.cdsStartPos=cdsStartPos;
		makeIntronTree();
	}
	
	public Gene(Collection<? extends Annotation> blocks, String name) {
		super(blocks, name);
		makeIntronTree();
	}
	
	public Gene(Annotation annot) {
		super(annot);
		makeIntronTree();
	}
	
	public Gene(Annotation annot, int cdsStartPos, int cdsEndPos) {
		super(annot);
		this.cdsEndPos=cdsEndPos;
		this.cdsStartPos=cdsStartPos;
		makeIntronTree();
	}

	private void makeIntronTree() {
		this.intronTree=new IntervalTree<SingleInterval>();
		for(SingleInterval intron: this.getIntronSet()){
			intronTree.put(intron.getReferenceStartPosition(), intron.getReferenceEndPosition(), intron);
		}
		
	}
	
	public Sequence getSequence(Sequence genomeSeq){
		return genomeSeq.getSubsequence(this);
	}

	/**
	 * @return A BlockedAnnotation representing the coding region of the gene
	 */
	public Annotation getCodingRegion(){
		if(cdsStartPos==NO_CDS || cdsEndPos==NO_CDS){return null;}
		SingleInterval cds=new SingleInterval(getReferenceName(), cdsStartPos, cdsEndPos);
		Annotation rtrn= intersect(cds);
		rtrn.setName(this.getName());
		return rtrn;
	}
	
	/**
	 * @return The 5'-UTR
	 */
	public Annotation get5UTR() {
		if(cdsStartPos==NO_CDS || cdsEndPos==NO_CDS || cdsStartPos == cdsEndPos) {
			return this;
		}
		Strand orientation = getOrientation();
		if(orientation.equals(Strand.POSITIVE)) {
			if(cdsStartPos == getReferenceStartPosition()) {
				return null;
			}
			SingleInterval utr = new SingleInterval(getReferenceName(), getReferenceStartPosition(), cdsStartPos);
			Annotation rtrn = intersect(utr);
			return new Gene(rtrn.getBlockSet(), -1, -1, this.getName() + "_5UTR");
		}
		if(orientation.equals(Strand.NEGATIVE)) {
			if(cdsEndPos == getReferenceEndPosition()) {
				return null;
			}
			SingleInterval utr = new SingleInterval(getReferenceName(), cdsEndPos, getReferenceEndPosition());
			Annotation rtrn =  intersect(utr);
			return new Gene(rtrn.getBlockSet(), -1, -1, this.getName() + "_5UTR");
		}
		throw new IllegalArgumentException("Can't get 5'-UTR because gene strand is unknown.");
	}
	
	
	public int get5UTRSize() {
		Annotation a=get5UTR();
		if(a==null) {return 0;}
		return a.size();
	}
	
	/**
	 * @return The 3'-UTR
	 */
	public Annotation get3UTR() {
		if(cdsStartPos==NO_CDS || cdsEndPos==NO_CDS || cdsStartPos == cdsEndPos) {
			return this;
		}
		Strand orientation = getOrientation();
		if(orientation.equals(Strand.POSITIVE)) {
			if(cdsEndPos == getReferenceEndPosition()) {
				return null;
			}
			SingleInterval utr = new SingleInterval(getReferenceName(), cdsEndPos, getReferenceEndPosition());
			//return intersect(utr);
			Annotation rtrn = intersect(utr);
			return new Gene(rtrn.getBlockSet(), -1, -1, this.getName() + "_3UTR");
		}
		if(orientation.equals(Strand.NEGATIVE)) {
			if(cdsStartPos == getReferenceStartPosition()) {
				return null;
			}
			SingleInterval utr = new SingleInterval(getReferenceName(), getReferenceStartPosition(), cdsStartPos);
			//return intersect(utr);
			Annotation rtrn = intersect(utr);
			return new Gene(rtrn.getBlockSet(), -1, -1, this.getName() + "_3UTR");
		}
		throw new IllegalArgumentException("Can't get 3'-UTR because gene strand is unknown.");
	}
	
	@Override
	public String toString(){
		return toBED(0,0,0);
	}
	
	public String toBED(double score) {
		return toBED(0,0,0,score);
	}
	
	@Override
	public String toBED(int r, int g, int b){
		return toBED(r, g, b, 0.0);
	}
		
	@Override
	public String toBED(int r, int g, int b, double score){
		if(r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255) {
			throw new IllegalArgumentException("RGB values must be between 0 and 255");
		}
		String rgb = r + "," + g + "," + b;
		Iterator<SingleInterval> exons = getBlocks();
		String rtrn=getReferenceName()+"\t"+getReferenceStartPosition()+"\t"+getReferenceEndPosition()+"\t"+getName() +"\t" + score + "\t"+getOrientation()+"\t"+this.cdsStartPos+"\t"+this.cdsEndPos+"\t"+rgb+"\t"+getNumberOfBlocks();
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
	
	//TODO wite tests
	public Collection<SingleInterval> getIntronSet()
	{
		Collection<SingleInterval> introns = new ArrayList<SingleInterval>();
		Iterator<SingleInterval> blocks = getBlocks();
		SingleInterval first;
		if(blocks.hasNext()){ first = blocks.next(); }
		else{ return introns; }
		
		while(blocks.hasNext())
		{
			SingleInterval next = blocks.next();
			SingleInterval intron = new SingleInterval(getReferenceName(),first.getReferenceEndPosition(),next.getReferenceStartPosition());
			intron.setName(getName()+"_intron");
			intron.setOrientation(getOrientation());
			introns.add(intron);
			first = next;
		}
		
		return introns;
	}

	public Collection<Annotation> getExonSet() {
		return getBlockSet();
		
	}

	public int distanceFrom5PrimeSpliceSite(int positionOfPolymerase) {
		
		
		//else
		Annotation exon=get3PrimeExon(positionOfPolymerase);
		
		//System.err.println("3' Exon: "+exon.toUCSC());
		
		//If positionOfPolymerase is in first exon or first intron --> d=-1
		if(exon.equals(getFirstExon())){return -1;}
		
		if(this.getOrientation().equals(Strand.POSITIVE)){
			return positionOfPolymerase-exon.getReferenceStartPosition();
		}
		
		if(this.getOrientation().equals(Strand.NEGATIVE)){
			return exon.getReferenceEndPosition()-positionOfPolymerase;
		}
		
		return -2;
		
		
	}
	
	public int distanceFrom5PrimeEnd(int position) {
		SingleInterval region=new SingleInterval(this.getReferenceName(), position, position+1, this.getOrientation());
		if(this.getOrientation().equals(Strand.POSITIVE)){
			//go from beginning until we find position
			boolean foundPosition=false;
			int length=0;
			Iterator<Annotation> exons=this.getExonSet().iterator();
			
			while(!foundPosition && exons.hasNext()){
				Annotation exon=exons.next();
				if(exon.contains(region)){
					foundPosition=true;
					length+=position-exon.getReferenceStartPosition();
				}
				else{length+=exon.size();}
			}
			return length;
		}
		
		
		if(this.getOrientation().equals(Strand.NEGATIVE)){
			//go from beginning until we find position
			boolean started=false;
			int length=0;
			Iterator<Annotation> exons=this.getExonSet().iterator();
			
			while(exons.hasNext()){
				Annotation exon=exons.next();
				if(exon.contains(region)){
					started=true;
					length+=exon.getReferenceEndPosition()-position;
				}
				else if(started){
					length+=exon.size();
				}
			}
			return length;
		}
		
		return -1;
		
	}
	
	public int distanceFrom3PrimeEnd(int position) {
		SingleInterval region=new SingleInterval(this.getReferenceName(), position, position+1, this.getOrientation());
		if(this.getOrientation().equals(Strand.NEGATIVE)){
			//go from beginning until we find position
			boolean foundPosition=false;
			int length=0;
			Iterator<Annotation> exons=this.getExonSet().iterator();
			
			while(!foundPosition && exons.hasNext()){
				Annotation exon=exons.next();
				if(exon.contains(region)){
					foundPosition=true;
					length+=position-exon.getReferenceStartPosition();
				}
				else{length+=exon.size();}
			}
			return length;
		}
		
		
		if(this.getOrientation().equals(Strand.POSITIVE)){
			//go from beginning until we find position
			boolean started=false;
			int length=0;
			Iterator<Annotation> exons=this.getExonSet().iterator();
			
			while(exons.hasNext()){
				Annotation exon=exons.next();
				if(exon.contains(region)){
					started=true;
					length+=exon.getReferenceEndPosition()-position;
				}
				else if(started){
					length+=exon.size();
				}
			}
			return length;
		}
		
		return -1;
		
	}
	
	
	public int distanceFrom5PrimeSpliceSite(Annotation fragment) {
		int positionOfPolymerase=fragment.get3PrimePosition();
		Annotation exon=get3PrimeExon(positionOfPolymerase);
		
		//If positionOfPolymerase is in first exon or first intron --> d=-1
		if(exon==null || exon.equals(getFirstExon())){return -1;}
		
		if(this.getOrientation().equals(Strand.POSITIVE)){
			int spliceSite=exon.getReferenceStartPosition();
			//If fragment overlaps exon.getStart
			if(fragment.getReferenceStartPosition()<=spliceSite && fragment.getReferenceEndPosition()>=spliceSite){
				return positionOfPolymerase-exon.getReferenceStartPosition();
			}
			
			return -3;
			
		}
		
		if(this.getOrientation().equals(Strand.NEGATIVE)){
			int spliceSite=exon.getReferenceEndPosition();
			if(fragment.getReferenceStartPosition()<=spliceSite && fragment.getReferenceEndPosition()>=spliceSite){
				return exon.getReferenceEndPosition()-positionOfPolymerase;
			}
			return -3;
		}
		
		return -2;
		
		
	}
	
	
	
	

	public Annotation get3PrimeExon(int positionOfPolymerase) {
		if(this.getOrientation().equals(Strand.POSITIVE)){
			return this.getMaxBlock(positionOfPolymerase);
		}
		else if(this.getOrientation().equals(Strand.NEGATIVE)){
			return this.getMinBlock(positionOfPolymerase);
		}
		else{
			throw new IllegalArgumentException("Strand must be Positive or Negative");
		}
	}

	private Annotation getFirstExon() {
		if(this.getOrientation().equals(Strand.POSITIVE)){
			return this.getFirstBlock();
		}
		else if(this.getOrientation().equals(Strand.NEGATIVE)){
			return this.getLastBlock();
		}
		else{
			throw new IllegalArgumentException("Strand must be Positive or Negative");
		}
	}

	public boolean inIntron(int positionOfPolymerase) {
		return !this.isWithinBlock(positionOfPolymerase);
	}

	public boolean overlapsExon(SingleInterval region) {
		return this.getOverlappingBlocks(region).hasNext();
	}

	public Iterator<SingleInterval> getOverlappingIntrons(Annotation fragment) {
		return this.intronTree.overlappingValueIterator(fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition());
	}

	public int getNumberOfOverlappingIntrons(Annotation fragment) {
		return this.intronTree.numOverlappers(fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition());
	}

	public int hasOneOverlappingIntron(Annotation fragment) {
		Iterator<SingleInterval> iter=getOverlappingIntrons(fragment);
		if(iter.hasNext()){
			iter.next();
			if(iter.hasNext()){return 2;}
			return 1;
		}
		return 0;
	}
	
	public int hasOneOverlappingExon(Annotation fragment) {
		Iterator<SingleInterval> iter=getOverlappingBlocks(fragment);
		if(iter.hasNext()){
			iter.next();
			if(iter.hasNext()){return 2;}
			return 1;
		}
		return 0;
	}

	public boolean hasCodingRegion() {
		if(cdsStartPos==NO_CDS || cdsEndPos==NO_CDS){return false;}
		if(cdsStartPos==cdsEndPos){return false;}
		return true;
	}

	public Gene extend5Prime(int extension) {
		SingleInterval extend=new SingleInterval(this.getReferenceName(), (this.get5PrimePosition()-extension), this.get5PrimePosition(), this.getOrientation());
		
		if(this.getOrientation().equals(Strand.NEGATIVE)){
			extend=new SingleInterval(this.getReferenceName(), this.get5PrimePosition(), (this.get5PrimePosition()+extension), this.getOrientation());
		}
		return new Gene(extend);
	}

	public Gene extend3Prime(int extend) {
		SingleInterval extension=new SingleInterval(this.getReferenceName(), this.get3PrimePosition(), (this.get3PrimePosition()+extend), this.getOrientation());
		if(this.getOrientation().equals(Strand.NEGATIVE)){
			extension=new SingleInterval(this.getReferenceName(), (this.get3PrimePosition()-extend), this.get3PrimePosition(), this.getOrientation());
		}
		
		return new Gene(extension);
	}

	public Annotation get5Prime(int length) {
		return get5Prime(0, length);
	}
	
	public Annotation get5Prime(int startOffSet, int length) {
		int start=startOffSet;
		int end=length+start;
		SingleInterval feature=new SingleInterval(this.getName(), start, end, Strand.POSITIVE);
		Annotation reference=convertToReferenceSpace(feature);
		reference.setName(getName());
		reference.setOrientation(getOrientation());
		return reference;
	}

	public Annotation get3Prime(int endOffSet, int length) {
		
		int end=size()-endOffSet;
		int start=end-length;
		SingleInterval feature=new SingleInterval(this.getName(), start, end, Strand.POSITIVE);
		Annotation reference=convertToReferenceSpace(feature);
		reference.setName(getName());
		reference.setOrientation(getOrientation());
		return reference;
	}
	
	public Annotation get3Prime(int length) {
		return get3Prime(0, length);
	}

	public SingleInterval getGenomicRegion() {
		return new SingleInterval(getReferenceName(), getReferenceStartPosition(), getReferenceEndPosition(), getOrientation(), getName());
	}

	public SingleInterval get5PrimeExon() {
		SingleInterval rtrn=this.getLastBlock();
		if(this.getOrientation().equals(Strand.POSITIVE)){rtrn=this.getFirstBlock();}
		rtrn.setName(this.getName());
		rtrn.setOrientation(this.getOrientation());
		return rtrn;
	}

	public int get3UTRSize() {
		Annotation a=get3UTR();
		if(a==null) {return 0;}
		return a.size();
	}

	public String getChrNum() {
		return getReferenceName().replaceFirst("chr", "");
	}

	public Collection<Gene> getJunctions() {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		List<SingleInterval> blocks=getExonList();
		
		for(int i=0; i<blocks.size()-1; i++) {
			Collection<Annotation> exons=new TreeSet<Annotation>();
			exons.add(blocks.get(i));
			exons.add(blocks.get(i+1));
			String name=getName()+"_junction"+i;
			Gene junction=new Gene(exons, name);
			rtrn.add(junction);	
		}
		
		return rtrn;
	}

	private List<SingleInterval> getExonList() {
		List<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		
		Iterator<SingleInterval> iter=this.getBlocks();
		while(iter.hasNext()) {rtrn.add(iter.next());}
		
		return rtrn;
	}

	public String toBED(String name) {
		
		String rgb = 0 + "," + 0 + "," + 0;
		Iterator<SingleInterval> exons = getBlocks();
		String rtrn=getReferenceName()+"\t"+getReferenceStartPosition()+"\t"+getReferenceEndPosition()+"\t"+name +"\t" + 0 + "\t"+getOrientation()+"\t"+this.cdsStartPos+"\t"+this.cdsEndPos+"\t"+rgb+"\t"+getNumberOfBlocks();
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

	public Collection<SingleInterval> getAllWindows(int binSize) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		//rtrn.addAll(getSplicedWindows(binSize));
		
		Iterator<SingleInterval> iter=this.getBlocks();
		while(iter.hasNext()) {
			SingleInterval block=iter.next();
			
			for(int i=block.getReferenceStartPosition(); i<=block.getReferenceEndPosition()-binSize; i++) {
				SingleInterval w=new SingleInterval(block.getReferenceName(), i, i+binSize);
				w.setOrientation(getOrientation());
				rtrn.add(w);
			}
			
		}
		return rtrn;	
	}

	

	/*public Collection<Annotation> getSplicedWindows(int binSize) {
		
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//get relative position of junctions
		Collection<Integer> junctionPositions=getRelativePositionsOfJunctions();
		//for each junction make a bin -bin/2 to +/bin/2
		for(Integer pos: junctionPositions) {
			Annotation bin=trimByRelativePositions(pos-binSize/2, pos+binSize/2);
			rtrn.add(bin);
		}
		return rtrn;
	}*/

	

	

	

	

	
	
}

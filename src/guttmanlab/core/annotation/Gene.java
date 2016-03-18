package guttmanlab.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

/**
 * A simple extension of a BlockedAnnotation that keeps track of the start and end of a coding region
 * @author mguttman
 *
 */
public class Gene extends BlockedAnnotation{

	private static final int NO_CDS = -1;
	protected int cdsStartPos = NO_CDS;
	protected int cdsEndPos = NO_CDS;
	
	public Gene(Collection<Annotation> blocks, int cdsStartPos, int cdsEndPos, String name) {
		super(blocks, name);
		this.cdsEndPos=cdsEndPos;
		this.cdsStartPos=cdsStartPos;
	}
	
	public Gene(Annotation annot) {
		super(annot);
	}

	/**
	 * @return A BlockedAnnotation representing the coding region of the gene
	 */
	public Annotation getCodingRegion(){
		if(cdsStartPos==NO_CDS || cdsEndPos==NO_CDS){return null;}
		SingleInterval cds=new SingleInterval(getReferenceName(), cdsStartPos, cdsEndPos);
		return intersect(cds);
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
			return intersect(utr);
		}
		if(orientation.equals(Strand.NEGATIVE)) {
			if(cdsStartPos == getReferenceStartPosition()) {
				return null;
			}
			SingleInterval utr = new SingleInterval(getReferenceName(), getReferenceStartPosition(), cdsStartPos);
			return intersect(utr);
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
			introns.add(intron);
			first = next;
		}
		
		return introns;
	}
	
}

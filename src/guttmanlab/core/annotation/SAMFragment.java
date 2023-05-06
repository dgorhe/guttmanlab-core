package guttmanlab.core.annotation;

import guttmanlab.core.annotation.predicate.ReadFlag;
import guttmanlab.core.annotationcollection.AnnotationCollection;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.TextCigarCodec;

public class SAMFragment implements MappedFragment{

	private SAMRecord record;
	private boolean strandIsFirstOfPair; 
	private Annotation annotation;
	private Collection<? extends ReadFlag> readFlags;
	public static String SAM_NUM_HITS_TAG = "NH";

	/**
	 * @param record SAM record
	 */
	public SAMFragment(SAMRecord record){
		this(record, false);
	}
	
	/**
	 * 
	 * @param record The SAM Record
	 * @param strandIsFirstOfPair Whether to treat the first of pair read as the fragment strand
	 */
	public SAMFragment(SAMRecord record, boolean strandIsFirstOfPair){
		super();
		this.record=record;
		this.strandIsFirstOfPair=strandIsFirstOfPair;
	}
	
	@Override
	public String getName() {
		return record.getReadName();
	}
	
	@Override
	public Iterator<SingleInterval> getBlocks() {
		return getAnnotation().getBlocks();
	}
	
	public Annotation getAnnotation(){
		if(this.annotation!=null){return this.annotation;}
		else{
			return parseCigar(record.getCigarString(), record.getReferenceName(), record.getAlignmentStart()-1, getOrientation(), getName()); 
		}
	}

	
	public static SingleInterval getSingleInterval(SAMRecord r) {
		return new SingleInterval(r.getReferenceName(), r.getAlignmentStart()-1, r.getAlignmentEnd());
	}
	
	
	@Override
	public String getReferenceName() {
		String chr=record.getReferenceName();
		if(!record.getReferenceName().contains("chr")){chr="chr"+chr;}
		return chr;
	}

	/**
	 * Returns the start position of this annotation in our coordinate space
	 * SAM coordinates are 1-based and inclusive whereas all of our objects are 0-based exclusive
	 */
	@Override
	public int getReferenceStartPosition() {
		//return record.getAlignmentStart()-1;
		return getAnnotation().getReferenceStartPosition();
	}

	@Override
	public int getReferenceEndPosition() {
		return getAnnotation().getReferenceEndPosition();
		//return record.getAlignmentEnd();  //this method uses an incorrect cigar parser
		//return record.getAlignmentStart() + this.size()-1;
	}
	
	/**
	 * Return the SAM Record object
	 * @return Original SAMRecord object
	 */
	public SAMRecord getSamRecord(SAMFileHeader header) {
		return record;
	}
	
	/**
	 * @return The SAM record
	 */
	public SAMRecord getSamRecord(){
		return record;
	}
	
	public int numberOfSplices() {
		Cigar cigar = TextCigarCodec.getSingleton().decode(record.getCigarString());
    	List<CigarElement> elements=cigar.getCigarElements();
		
    	int rtrn=0;
		for(CigarElement element: elements){
			CigarOperator op=element.getOperator();
			
			if(op.equals(CigarOperator.N)){
				rtrn++;
			}
		}
			return rtrn;
	}
	
	 /**
     * Populate an annotation from a Cigar string
     * @param cigarString Cigar string
     * @param chr Fragment reference sequence
     * @param start Fragment start
	 * @param strand Fragment strand
	 * @param name Name of annotation to return
     * @return A blocked annotation
     */
	public static Annotation parseCigar(String cigarString, String chr, int start, Strand strand, String name) {
		//String oldCigar=cigarString;
		//cigarString=compressCIGAR(cigarString);
		
		//System.out.println(oldCigar+"\t"+cigarString);
		
    	Cigar cigar = TextCigarCodec.getSingleton().decode(cigarString);
    	List<CigarElement> elements=cigar.getCigarElements();
		
    	BlockedAnnotation rtrn=new BlockedAnnotation(name);
    	
		int currentOffset = start;
		
		for(CigarElement element: elements){
			CigarOperator op=element.getOperator();
			int length=element.getLength();
			
			//then lets create a block from this
			if(op.equals(CigarOperator.MATCH_OR_MISMATCH)){
				int blockStart=currentOffset;
				int blockEnd=blockStart+length;
				rtrn.addBlocks(new SingleInterval(chr, blockStart, blockEnd, strand, name));
				currentOffset=blockEnd;
			}
			else if(op.equals(CigarOperator.S)){
				
			}
			else if(op.equals(CigarOperator.N)){ //This is spliced
				int blockStart=currentOffset;
				int blockEnd=blockStart+length;
				currentOffset=blockEnd;
			}
			else if(op.equals(CigarOperator.INSERTION) ||  op.equals(CigarOperator.H) || op.equals(CigarOperator.DELETION)|| op.equals(CigarOperator.SKIPPED_REGION)){
				currentOffset+=length;
			}
		}
		
		return rtrn;
	}
	
	
	
	public static SAMRecord compressCIGAR(SAMRecord record) {
		String cigar=record.getCigarString();
		String newCIGAR=compressCIGAR(cigar);
		
		record.setCigarString(newCIGAR);
		record.setReadString("*");
		record.setBaseQualityString("*");
		return record;
	}
	
	private static String compressCIGAR(String cigarString) {
		Cigar cigar = TextCigarCodec.getSingleton().decode(cigarString);
    	List<CigarElement> elements=cigar.getCigarElements();
    	
    	List<CigarElement> oldelements=new ArrayList<CigarElement>();
    	oldelements.addAll(elements);
    	elements=remove(elements, CigarOperator.I);
    	elements=remove(elements, CigarOperator.H);
    	
    	CigarElement first=elements.get(0);
    	CigarElement last=elements.get(elements.size()-1);
    	String lastString="";
    	
    	String rtrn="";
    	if(first.getOperator().equals(CigarOperator.S)) {
    		rtrn+=first.getLength()+"S";
    		elements.remove(0);
    	}
    	
    	if(last.getOperator().equals(CigarOperator.S)) {
    		elements.remove(elements.size()-1);
    		lastString=last.getLength()+"S";
    	}
    	
    	
    	//if N split list into 2
    	//Go through each and when hit an N make a new list
    	List<Integer> indexOfN=new ArrayList<Integer>();
    	for(int i=0; i<elements.size(); i++) {
    		CigarElement e=elements.get(i);
    		if(e.getOperator().equals(CigarOperator.N)) {indexOfN.add(i);}
    	}
    	
    	int start=0;
    	for(Integer index: indexOfN) {
    		int length=merge(start, index, elements); //TODO still need to deal with S
    		rtrn+=length+"M"+elements.get(index).getLength()+"N";
    		start=index+1;
    	}
    	
    	int length=merge(start, elements.size(), elements);
    	rtrn+=length+"M";
    	
    	
    	rtrn+=lastString;
    	
    	//TODO pop first and last if S write
    	
    	return rtrn;
    	
    	/*
    	
    	System.out.println(cigarString+" "+oldelements.size()+" "+elements.size());
		
    	String rtrn="";
		
		int cumLength=0;
		for(int i=0; i<elements.size(); i++) {
			CigarElement element=elements.get(i);
			CigarOperator op=element.getOperator();
			int length=element.getLength();
			
			if(op.equals(CigarOperator.S)){
				if(cumLength>0) {rtrn+=cumLength+"M";}
				rtrn+=length+""+op;
			}
			
			if(op.equals(CigarOperator.MATCH_OR_MISMATCH)){
				cumLength+=length;
				//keep going until hit N or S
				i=i+1;
				if(i<elements.size()) {
					element=elements.get(i);
					op=element.getOperator();
					if(op.equals(CigarOperator.N)) {rtrn+=cumLength+"M"+element.getLength()+"N"; cumLength=0;}
					else if(op.equals(CigarOperator.S)) {rtrn+=cumLength+"M"+element.getLength()+"S"; cumLength=0;}
					else {
						cumLength+=element.getLength();
					}
				}
				else {rtrn+=cumLength+"M";}
				
			}
			if(cumLength>0) {rtrn+=cumLength+"M";}
			
		}
		return rtrn;*/
		
	}

	private static int merge(int start, int end, List<CigarElement> elements) {
		int sum=0;
		for(int i=start; i<end; i++) {
			sum+=elements.get(i).getLength();
		}
		return sum;
	}

	private static List<CigarElement> remove(List<CigarElement> elements, CigarOperator excludeOp) {
		List<CigarElement> rtrn=new ArrayList<CigarElement>();
		
		for(CigarElement e: elements) {
			if(!e.getOperator().equals(excludeOp)) {rtrn.add(e);}
		}
		
		return rtrn;
	}

	@Override
	/**
	 * Use strand info from instantiation
	 */
	public Strand getOrientation() {
		Strand rtrn=Annotation.Strand.POSITIVE;
		if(this.record.getReadNegativeStrandFlag()){rtrn=Annotation.Strand.NEGATIVE;}
		if((this.isPaired() && this.strandIsFirstOfPair && !this.record.getFirstOfPairFlag()) || (this.isPaired() && !this.strandIsFirstOfPair && this.record.getFirstOfPairFlag())){rtrn=rtrn.getReverseStrand();}
		return rtrn;
	}
	
	@Override
	public int getNumberOfBlocks() {
		return getAnnotation().getNumberOfBlocks();
	}

	@Override
	public int size() {
		return getAnnotation().size();
	}

	@Override
	public int getRelativePositionFrom5PrimeOfFeature(int referenceStart) {
		return getAnnotation().getRelativePositionFrom5PrimeOfFeature(referenceStart);
	}

	@Override
	public Collection<? extends ReadFlag> getReadFlags() {
		return readFlags();
	}
	
	private Collection<? extends ReadFlag> readFlags(){
		//IF already parsed, return the collection
		if(readFlags!=null){return readFlags;}
		
		//Else, parse it, save, and return
		else{
			readFlags=parseFlags();
			return readFlags;
		}
	}

	private Collection<? extends ReadFlag> parseFlags() {
		throw new UnsupportedOperationException("TODO");
	}

	/**
	 * @return whether this read is part of a pair
	 */
	public boolean isPaired() {
		return record.getReadPairedFlag();
	}

	/**
	 * @return The location of the alignment start of the mate
	 */
	public int getMateReferenceStart() {
		return record.getMateAlignmentStart()-1;
	}

	/**
	 * The reference of the mate
	 * @return
	 */
	public String getMateReferenceName() {
		return record.getMateReferenceName();
	}

	/**
	 * @return The orientation of the read
	 */
	public Strand getReadOrientation() {
		if(this.getSamRecord().getReadNegativeStrandFlag()){return Strand.NEGATIVE;}
		return Strand.POSITIVE;
	}

	@Override
	public void setOrientation(Strand orientation) {
		throw new UnsupportedOperationException("Defined by parts");
	}
	
	/**
	 * Get the value of an int SAM tag
	 * @param tag
	 * @return From picard documentation: returns the value of a tag or throws RuntimeException if the value is not an Integer type or will not fit in an integer
	 */
	public int getIntTag(String tag) {
		return record.getIntegerAttribute(tag).intValue();
	}
	
	/**
	 * Get the value of a string SAM tag
	 * @param tag
	 * @return The value of the tag
	 */
	public String getStringTag(String tag) {
		return record.getStringAttribute(tag);
	}

	@Override
	public AnnotationCollection<DerivedAnnotation<? extends Annotation>> getWindows(
			int windowSize, int stepSize) {
		throw new UnsupportedOperationException();
	}

	@Override
	public int getNumHits() {
		return getIntTag(SAM_NUM_HITS_TAG);
	}

	@Override
	public int getMappingQuality() {
		return record.getMappingQuality();
	}

	@Override
	public void setName(String name) {
		// TODO Auto-generated method stub
		throw new IllegalStateException("setName() not implemented for SAMFragment");
	}

	@Override
	public Annotation trimByRelativePositions(int relativeStart, int relativeEnd) {
		// TODO Auto-generated method stub
		throw new IllegalStateException("trimByRelativePositions() not implemented for SAMFragment");
	}

	public char getSequenceAtPosition(String chr, int pos) {
		if(!this.getReferenceName().equals(chr)){return ' ';}
		char[] seq=getSamRecord().getReadString().toCharArray();
		for(int i=0; i<seq.length; i++){
			if(getSamRecord().getReferencePositionAtReadPosition(i) == pos){return seq[i];}
		}
		return ' ';
	}
	
	public static char getSequenceAtPosition(htsjdk.samtools.SAMRecord record2, String chr, int pos) {
		if(!record2.getReferenceName().equals(chr)){return ' ';}
		char[] seq=record2.getReadString().toCharArray();
		for(int i=0; i<seq.length; i++){
			if(record2.getReferencePositionAtReadPosition(i) == pos){return seq[i];}
		}
		return ' ';
	}

	public SingleInterval getGenomicInterval() {
		SingleInterval rtrn=new SingleInterval(this.getReferenceName(), this.getReferenceStartPosition(), this.getReferenceEndPosition(), this.getOrientation());
		return rtrn;
	}

	/*public static Collection<Annotation> allBins(SAMRecord read, int resolution, int stepSize) {
		SAMFragment frag=new SAMFragment(read);
		return frag.getWindowCollection(resolution, resolution);
	}*/
	
	public static Collection<SingleInterval> allBins(SAMRecord read, int resolution) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		SAMFragment f=new SAMFragment(read);
		Iterator<SingleInterval> blocks=f.getBlocks();
		while(blocks.hasNext()) {
			rtrn.addAll(blocks.next().allBins(resolution, f.getOrientation()));
		}
		
		return rtrn;
		//return allBins(read, resolution, resolution);
	}
	
	public static Collection<SingleInterval> getAllWindows(SAMRecord read, int binSize) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		SAMFragment f=new SAMFragment(read);
		Iterator<SingleInterval> iter=f.getBlocks();
		while(iter.hasNext()) {
			SingleInterval block=iter.next();
			//System.err.println(block.toUCSC());
			for(int i=block.getReferenceStartPosition()-(binSize-1); i<block.getReferenceEndPosition()+binSize; i++) {
				SingleInterval w=new SingleInterval(block.getReferenceName(), i, i+binSize);
				w.setOrientation(f.getOrientation());
				rtrn.add(w);
			}
			
		}
		return rtrn;
	}
	

	public static Collection<SingleInterval> allBins(SAMRecord read, int resolution, int stepSize) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		SAMFragment f=new SAMFragment(read);
		Iterator<SingleInterval> blocks=f.getBlocks();
		while(blocks.hasNext()) {
			rtrn.addAll(blocks.next().allBins(resolution, stepSize, f.getOrientation()));
		}
		
		return rtrn;
		//return allBins(read, resolution, resolution);
	}

	public static boolean isSpliced(SAMRecord read) {
		return read.getCigarString().contains("N");
	}
	
	@Override
	public String getCigarString(){
		return this.record.getCigarString();
	}
		
	public boolean isSpliced() {
		return this.getCigarString().contains("N");
		//return frag.getNumberOfBlocks()>1;
	}

	public static Strand getOrientation(SAMRecord record) {
		Strand rtrn=Annotation.Strand.POSITIVE;
		if(record.getReadNegativeStrandFlag()){rtrn=Annotation.Strand.NEGATIVE;}
		if((record.getReadPairedFlag() && record.getFirstOfPairFlag())){rtrn=rtrn.getReverseStrand();}
		return rtrn;
	}

	public boolean contains(int position) {
		if(this.getReferenceStartPosition()<position && this.getReferenceEndPosition()>position){return true;}
		return false;
	}
	
	
	

	
	
	//TODO For intersect, merge, and convert --> override and add all ReadFlags to the new object

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

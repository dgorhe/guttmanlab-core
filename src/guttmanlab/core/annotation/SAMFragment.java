package guttmanlab.core.annotation;

import guttmanlab.core.annotation.predicate.ReadFlag;
import guttmanlab.core.annotationcollection.AnnotationCollection;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.TextCigarCodec;

public class SAMFragment extends AbstractAnnotation implements MappedFragment{

	private SAMRecord record;
	private boolean strandIsFirstOfPair; 
	private Annotation annotation;
	private Collection<? extends ReadFlag> readFlags;
	public static String SAM_NUM_HITS_TAG = "NH";

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
	
	private Annotation getAnnotation(){
		if(this.annotation!=null){return this.annotation;}
		else{
			return parseCigar(record.getCigarString(), record.getReferenceName(), record.getAlignmentStart()-1, getOrientation(), getName()); 
		}
	}

	@Override
	public String getReferenceName() {
		return record.getReferenceName();
	}

	/**
	 * Returns the start position of this annotation in our coordinate space
	 * SAM coordinates are 1-based and inclusive whereas all of our objects are 0-based exclusive
	 */
	@Override
	public int getReferenceStartPosition() {
		return record.getAlignmentStart()-1;
	}

	@Override
	public int getReferenceEndPosition() {
		//return record.getAlignmentEnd();  //this method uses an incorrect cigar parser
		return record.getAlignmentStart() + this.size()-1;
	}
	
	/**
	 * Return the SAM Record object
	 * @return Original SAMRecord object
	 */
	@Override
	public SAMRecord getSamRecord(SAMFileHeader header) {
		return record;
	}
	
	public SAMRecord getSamRecord(){
		return record;
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
			else if(op.equals(CigarOperator.N)){
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

	public Strand getReadOrientation() {
		if(this.getSamRecord().getReadNegativeStrandFlag()){return Strand.NEGATIVE;}
		return Strand.POSITIVE;
	}

	@Override
	public void setOrientation(Strand orientation) {
		//Empty implementation since defined by parts
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
	
	//TODO For intersect, merge, and convert --> override and add all ReadFlags to the new object

	
}

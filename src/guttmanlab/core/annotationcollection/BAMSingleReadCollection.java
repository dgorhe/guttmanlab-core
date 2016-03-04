package guttmanlab.core.annotationcollection;


import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.predicate.ContainedByFilter;
import guttmanlab.core.annotation.predicate.OverlapsFilter;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import org.apache.commons.collections15.Predicate;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;

/**
 * This class represents a single-end read collection
 */
public class BAMSingleReadCollection extends AbstractAnnotationCollection<SAMFragment>{

	private SAMFileReader reader;
	private CoordinateSpace referenceSpace;
	private final File bamFile;
	
	/**
	 * Constructs a collection of single-read aligned fragments from a BAM file.
	 * @param bamFile is the BAM file containing the single-read alignments
	 */
	public BAMSingleReadCollection(File bamFile){
		super();
		this.reader = new SAMFileReader(bamFile);
		this.referenceSpace = new CoordinateSpace(reader.getFileHeader());
		this.bamFile = bamFile;
	}

	/**
	 * Constructs a collection of single-read aligned fragments from a BAM file.
	 * @param bamFilePath is the path of the BAM file containing the single-read alignments
	 */
	public BAMSingleReadCollection(String bamFilePath) {
		this(new File(bamFilePath));
	}
	
	@Override
	public CloseableIterator<SAMFragment> sortedIterator() {
		return new FilteredIterator<SAMFragment>(new WrappedIterator(reader.iterator()), getFilters());
	}
	
	@Override
	public CloseableIterator<SAMFragment> sortedIterator(Annotation region, boolean fullyContained) {

		// Create the interval hull of `region`. This step is redundant if `region` just has one block.
		Annotation hull = new SingleInterval(region.getReferenceName(),
				   region.getReferenceStartPosition(),
				   region.getReferenceEndPosition(),
				   region.getOrientation(),
				   region.getName());
		
		// Get the reads that overlap the interval hull. Not all of these reads will necessarily overlap the
		// original blocked interval.
		CloseableIterator<SAMFragment> iter = new WrappedIterator(reader.queryOverlapping(hull.getReferenceName(),
				   												  hull.getReferenceStartPosition() + 1,
				   												  hull.getReferenceEndPosition()));
		
		// Add existing filters. Also add an additional filter depending on 'fullyContained'.
		// Copy the filters ArrayList(), so we don't add filters to the original.
		Collection<Predicate<SAMFragment>> filters = new ArrayList<Predicate<SAMFragment>>(getFilters());
		if (fullyContained) {
			filters.add(new ContainedByFilter<SAMFragment>(region));
		} else {
			filters.add(new OverlapsFilter<SAMFragment>(region));
		}
		
		// TODO Check if we need the StrandFilter provided by the third argument. The overlap/contains methods
		// might deal with strandedness already.
		return new FilteredIterator<SAMFragment>(iter, filters, region.getOrientation());
	}
		
	
	public void writeToFile(String fileName) {
		CloseableIterator<SAMFragment> iter= sortedIterator();
		writeToFile(fileName, iter);
	}
	
	private void writeToFile(String fileName, CloseableIterator<SAMFragment> iter){
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(this.reader.getFileHeader(), false, new File(fileName));
		
		while(iter.hasNext()){
			SAMFragment ann=iter.next();
			writer.addAlignment(ann.getSamRecord());
		}
		writer.close();
		iter.close();
	}
	
	public void writeToFile(String fileName, Annotation region) {
		CloseableIterator<SAMFragment> iter= sortedIterator(region, false);
		writeToFile(fileName, iter);
	}

	@Override
	public CoordinateSpace getReferenceCoordinateSpace() {
		return this.referenceSpace;
	}
	
	/**
	 * Gets the file header of this collection's BAM file
	 * @return the file header of this collection's BAM file
	 */
	public SAMFileHeader getFileHeader() {
		return reader.getFileHeader();
	}

	public PairedMappedFragment<SAMFragment> findReads(SAMFragment fragment) {
		//TODO A few ideas about how to implement this, simplest, just look up alignment start and alignment end and match names
		SAMRecordIterator alignment=this.reader.queryAlignmentStart(fragment.getSamRecord().getReferenceName(), fragment.getSamRecord().getAlignmentStart());
		SAMFragment read1=findRead(alignment, fragment.getName());;
		
		SAMRecordIterator mate=this.reader.queryAlignmentStart(fragment.getSamRecord().getReferenceName(), fragment.getSamRecord().getMateAlignmentStart());
		SAMFragment read2=findRead(mate, fragment.getName());
		
		PairedMappedFragment<SAMFragment> rtrn=new PairedMappedFragment<SAMFragment>(read1, read2);
		return rtrn;
	}

	private SAMFragment findRead(SAMRecordIterator alignment, String name) {
		SAMFragment rtrn=null;
		while(alignment.hasNext()){
			SAMRecord record=alignment.next();
			if(record.getReadName().equalsIgnoreCase(name)){
				rtrn=new SAMFragment(record);
				break;
			}
		}
		alignment.close();
		return rtrn;
	}
	
	public File getBamFile() {
		return bamFile;
	}
	
	/**
	 * Gets a String representation of this collection of reads. Currently this is simply the
	 * basename of the BAM file, e.g., a BAM file "/home/user/test.bam" is represented as
	 * "test".
	 */
	@Override
	public String toString() {
		return bamFile.getName().split("\\.(?=[^\\.]+$)")[0];
	}
	
	/**
	 * A wrapper class for Picard's SAMRecordIterator.
	 */
	public class WrappedIterator implements CloseableIterator<SAMFragment>{

		SAMRecordIterator iter;
		
		/**
		 * Constructor which wraps the input SAMRecordIterator.
		 * @param iter the SAMRecordIterator to wrap
		 */
		public WrappedIterator(SAMRecordIterator iter){
			this.iter=iter;
		}

		@Override
		public boolean hasNext() {
			return iter.hasNext();
		}

		@Override
		public SAMFragment next() {
			return new SAMFragment(iter.next());
		}

		@Override
		public void remove() {
			iter.remove();
			
		}

		@Override
		public void close() {
			iter.close();
		}
	}
	
	/*public BAMSingleReadCollection convert(AnnotationCollection<? extends Annotation> features, boolean fullyContained){
		//Setup BAM File Writer
		CoordinateSpace space=features.getFeatureCoordinateSpace();
		SAMFileHeader header=space.getBAMFileHeader();
		File tmpFile=new File(System.currentTimeMillis()+"_tmpConvereted.bam");
		//tmpFile.deleteOnExit(); TODO Put back?
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(header, false, tmpFile);
				
		
		//Get iterator of pairs if exists
		CloseableIterator<SAMFragment> iter=sortedIterator();
		
		int counter=0;//TODO
		while(iter.hasNext()){
			SAMFragment original=iter.next();
			Collection<SAMRecord> converted=convertCoordinates(features, original, header); //TODO This is wildly inefficient
			for(SAMRecord record: converted){
				writer.addAlignment(record);
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		
		writer.close();
		iter.close();
		
		//Read in the new file and parse an AnnotationCollection
		return new BAMSingleReadCollection(tmpFile);
	}

	private Collection<SAMRecord> convertCoordinates(AnnotationCollection<? extends Annotation> features, SAMFragment original, SAMFileHeader featureHeader) {
		Collection<SAMRecord> rtrn=new ArrayList<SAMRecord>();
		
		//IF mate is null just use original, else get features that overlap both
		//Find features that overlap
		//Find features overlapping the annotation
		CloseableIterator<? extends Annotation> iter=features.sortedIterator(original, false);

		//Adjust the coordinates of the feature as needed in featureSpace (ie as distance from start and end)
		while(iter.hasNext()){
			Annotation feature=iter.next();
			Annotation convertedOriginal=feature.convertToFeatureSpace(original);
			Strand orientation=getOrientation(original, feature);
			
			if(convertedOriginal!=null){
				convertedOriginal.setOrientation(orientation);
				int mateStart=-1;
				if(original.getReferenceName().equalsIgnoreCase(original.getMateReferenceName())){
					mateStart=feature.getRelativePositionFrom5PrimeOfFeature(original.getMateReferenceStart());
				}
				SAMRecord record=buildSAMRecord(original, convertedOriginal, mateStart, featureHeader);
				rtrn.add(record);
			}
		}
		return rtrn;			
	}

	private Strand getOrientation(SAMFragment original, Annotation feature) {
		Strand orientation;
		if(original.getReadOrientation().equals(feature.getOrientation())){
			orientation=Strand.POSITIVE; //If orientations match, then new is Positive
		}
		else{orientation=Strand.NEGATIVE;} //If orientation don't match, then new is Negative
		return orientation;
	}

	private SAMRecord buildSAMRecord(SAMFragment original, Annotation convertedOriginal, int mateStart, SAMFileHeader featureHeader) {
		//Make new SAMRecord
		SAMRecord record=new SAMRecord(featureHeader);
		
		//Set reference name, start, cigar, strand, readname
		record.setReferenceName(convertedOriginal.getReferenceName());
		record.setAlignmentStart(convertedOriginal.getReferenceStartPosition()+1);
		record.setCigarString(convertedOriginal.getCigarString());
		record.setReadNegativeStrandFlag(convertedOriginal.getOrientation().equals(Strand.NEGATIVE));
		record.setReadName(original.getName());
		
		if(convertedOriginal.size()<0){
			System.err.println(convertedOriginal.getCigarString()+" "+original.getName()+" "+original.getReferenceName()+":"+original.getReferenceStartPosition()+"-"+original.getReferenceEndPosition()+" "+convertedOriginal.getReferenceName());
		}
		
		//Set sequence of the read
		//TODO This is really just to debug, it's not really needed
		//IF the feature was neg strand --> reverse compliment
		Sequence seq=new Sequence(original.getSamRecord().getReadString());
		if(original.getOrientation().equals(Strand.NEGATIVE)){seq=seq.reverseComplement();}
		record.setReadString(seq.getSequenceBases());
		
		//Set additional read flags
		record.setDuplicateReadFlag(original.getSamRecord().getDuplicateReadFlag());
		record.setMappingQuality(original.getSamRecord().getMappingQuality());
		record.setNotPrimaryAlignmentFlag(original.getSamRecord().getNotPrimaryAlignmentFlag());
		record.setBaseQualityString(original.getSamRecord().getBaseQualityString());
		
		record.setReadPairedFlag(original.getSamRecord().getReadPairedFlag());
		record.setProperPairFlag(original.getSamRecord().getProperPairFlag());
		record.setFirstOfPairFlag(original.getSamRecord().getFirstOfPairFlag());
		record.setSecondOfPairFlag(original.getSamRecord().getSecondOfPairFlag());
		
		//IF has mate, set mate reference name, start
		if(mateStart>-1){
			record.setMateReferenceName(convertedOriginal.getReferenceName());
			record.setMateAlignmentStart(mateStart+1);
		}
		else{record.setMateUnmappedFlag(true);}
		
		return record;
	}*/
}

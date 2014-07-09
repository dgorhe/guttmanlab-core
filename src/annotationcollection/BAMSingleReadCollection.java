package annotationcollection;


import java.io.File;
import java.util.ArrayList;
import java.util.Collection;

import sequence.Sequence;
import coordinatespace.CoordinateSpace;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;
import annotation.Annotation;
import annotation.Annotation.Strand;
import annotation.DerivedAnnotation;
import annotation.PairedMappedFragment;
import annotation.SAMFragment;
import annotation.SingleInterval;

/**
 * This class represents a single-end read collection
 * @author mguttman
 *
 */
public class BAMSingleReadCollection extends AbstractAnnotationCollection<SAMFragment>{

	private SAMFileReader reader;
	private CoordinateSpace referenceSpace;
	
	public BAMSingleReadCollection(File bamFile){
		super();
		this.reader=new SAMFileReader(bamFile);
		this.referenceSpace=new CoordinateSpace(reader.getFileHeader());
	}

	@Override
	public CloseableIterator<SAMFragment> sortedIterator() {
		return new FilteredIterator<SAMFragment>(new WrappedIterator(reader), getFilters());
	}

	@Override
	public CloseableIterator<SAMFragment> sortedIterator(Annotation region, boolean fullyContained) {
		return new FilteredIterator<SAMFragment>(new WrappedIterator(reader, region), getFilters());
	}

	//TODO Consider whether to delete
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
	
	public class WrappedIterator implements CloseableIterator<SAMFragment>{

		SAMRecordIterator iter;
		
		public WrappedIterator(SAMFileReader reader){
			this.iter=reader.iterator();
		}
		
		public WrappedIterator(SAMFileReader reader, Annotation region) {
			this.iter=reader.queryOverlapping(region.getReferenceName(), region.getReferenceStartPosition()+1, region.getReferenceEndPosition());
			
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
			// TODO Auto-generated method stub
			
		}

		@Override
		public void close() {
			iter.close();
		}
	}

	@Override
	public CoordinateSpace getReferenceCoordinateSpace() {
		return this.referenceSpace;
	}
	
	public SAMFileHeader getFileHeader(){return this.reader.getFileHeader();}

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
	
	public BAMSingleReadCollection convert(AnnotationCollection<? extends Annotation> features, boolean fullyContained){
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
	}
}
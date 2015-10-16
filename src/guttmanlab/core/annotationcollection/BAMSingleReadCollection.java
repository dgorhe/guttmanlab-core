package guttmanlab.core.annotationcollection;


import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.predicate.StrandFilter;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;

/**
 * This class represents a single-end read collection
 * @author mguttman
 *
 */
public class BAMSingleReadCollection extends AbstractAnnotationCollection<SAMFragment>{

	private SAMFileReader reader;
	private CoordinateSpace referenceSpace;
	private String bam;
	
	public BAMSingleReadCollection(File bamFile){
		super();
		this.reader=new SAMFileReader(bamFile);
		this.referenceSpace=new CoordinateSpace(reader.getFileHeader());
		this.bam = bamFile.getName();
	}

	public String getBamFile() {
		return bam;
	}
	
	public String toString() {
		return getBamFile().replaceAll(".bam","");
	}
	
	
	@Override
	public CloseableIterator<SAMFragment> sortedIterator() {
		return new FilteredIterator<SAMFragment>(new WrappedIterator(reader.iterator()), getFilters());
	}

	@Override
	public CloseableIterator<SAMFragment> sortedIterator(Annotation region, boolean fullyContained) {
		CloseableIteratorChain iter_chain = new CloseableIteratorChain(region,fullyContained);
		return new FilteredIterator<SAMFragment>(iter_chain, getFilters(),region.getOrientation());
	}
	
	public class CloseableIteratorChain implements Iterator<SAMFragment>{
		
		private CloseableIterator<SAMFragment> currentIterator;
		private Iterator<SingleInterval> blocks;
		private Annotation region;
		private SAMFragment next;
		private ArrayList<String> splicedReadNames;
		private boolean fullyContained;
		
		public CloseableIteratorChain(Annotation region)
		{
			this.region = region;
			this.blocks = region.getBlocks();
			this.currentIterator = null;
			this.splicedReadNames = new ArrayList<String>();
		}
		
		public CloseableIteratorChain(Annotation region, boolean fullyContained)
		{
			this.region = region;
			this.blocks = region.getBlocks();
			this.currentIterator = null;
			this.splicedReadNames = new ArrayList<String>();
			this.fullyContained = fullyContained;
		}
		
		public boolean hasNext(){
			if(currentIterator == null)
			{
				if(blocks.hasNext())
				{
					Annotation block = blocks.next();
					if(this.fullyContained)
						currentIterator = new WrappedIterator(reader.queryContained(region.getReferenceName(), block.getReferenceStartPosition()+1,block.getReferenceEndPosition()));
					else
						currentIterator = new WrappedIterator(reader.queryOverlapping(region.getReferenceName(), block.getReferenceStartPosition()+1,block.getReferenceEndPosition()));
					return hasNext();
				}
				else //there were no more blocks
					return false;
			}
			else
			{
				findNext();
				if(next!=null)
				{
					return true;
				}
				else //we've reached the end of the current iterator
				{
					currentIterator.close();
					currentIterator = null;
					return hasNext();
				}
			}
		}
		
		private void findNext()
		{
			if(next==null && currentIterator.hasNext())
			{
				SAMFragment n = currentIterator.next();
				//TODO find out what this was for, maybe add back
//				if(n.getNumberOfBlocks()>1)
//				{
//					boolean firstOfPair = n.getSamRecord().getReadPairedFlag() ? n.getSamRecord().getFirstOfPairFlag() : true;
//					String id = ""+n.getName()+firstOfPair;
//					if(!splicedReadNames.contains(id))
//					{
//						this.next = n;
//						splicedReadNames.add(id);
//					}
//					else
//						findNext();
//				}
//				else
					this.next = n;
			}
		}
		
		@Override
		public SAMFragment next() {
			SAMFragment n = next;
			this.next = null;
			return n;
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
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

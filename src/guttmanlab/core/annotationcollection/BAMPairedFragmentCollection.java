package guttmanlab.core.annotationcollection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.predicate.ContainedByFilter;
import guttmanlab.core.annotation.predicate.OverlapsFilter;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.Pair;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeMap;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloseableIterator;

/**
 * This class represents a collection of mapped paired-end fragments.
 */
public class BAMPairedFragmentCollection extends AbstractAnnotationCollection<PairedMappedFragment<SAMFragment>>{

	//TODO Override coordinate conversion by converting the single reads and then parsing this BAM file into the paired fragment
	
	private static final String EXTENSION = ".pe.bam";
	private static final String ALIGNMENT_CIGAR = "aC";    // Note: Keep the lowercase letter in these tags.
	private static final String MATE_CIGAR = "mC";         // Lowercase tags are reserved for custom use, and
	private static final String MATE_MAPPING_QUALITY="mQ"; // will never be used in future implementations.
	private File bamFile;
	private File fragmentFile;
	private BAMSingleReadCollection reads;
	private SpecialBAMPECollection fragmentReader;
	private static Logger logger = Logger.getLogger(BAMPairedFragmentCollection.class.getName());
	
	/**
	 * Constructs a collection of paired-end aligned fragments from a BAM file. Constructing this
	 * object will create a temporary file which is essentially a modified BAM file that will support
	 * iterating by matching read pairs. This temporary "fragment" file may be quite large.
	 * @param bamFile the BAM file containing the paired-end alignments
	 * @throws IOException if the temporary fragment file cannot be written
	 */
	public BAMPairedFragmentCollection(File bamFile) throws IOException {
		reads = new BAMSingleReadCollection(bamFile);
		this.fragmentFile = makeFragmentFile(bamFile);
		this.bamFile = bamFile;
	}
	
	/**
	 * Constructs a collection of paired-end aligned fragments from a BAM file. Constructing this
	 * object will create a temporary file which is essentially a modified BAM file that will support
	 * iterating by matching read pairs. This temporary "fragment" file may be quite large.
	 * @param bamFilePath is the BAM file path
	 * @throws IOException if the temporary fragment file cannot be written
	 */
	public BAMPairedFragmentCollection(String bamFilePath) throws IOException {
		this(new File(bamFilePath));
	}

	/**
	 * Constructs a collection of paired-end aligned fragments from a collection of
	 * single-read aligned fragments. Constructing this object will create a temporary file which is
	 * essentially a modified BAM file that will support iterating by matching read pairs. This\
	 * temporary "fragment" file may be quite large.
	 * @param reads the single-read fragment collection
	 * @throws IOException if the temporary fragment file cannot be written
	 */
	public BAMPairedFragmentCollection(BAMSingleReadCollection reads) throws IOException {
		this.reads = reads;
		this.fragmentFile = makeFragmentFile();
		bamFile = reads.getBamFile();
	}
	
	public File getBamFile() {
		return bamFile;
	}
	
	private File makeFragmentFile() throws IOException {
		File file = File.createTempFile("temp", EXTENSION);
		file.deleteOnExit();
		return file;
	}

	private File makeFragmentFile(File bamFile) throws IOException {
		String baseName = bamFile.getName().split("\\.(?=[^\\.]+$)")[0];
		File file = File.createTempFile(baseName, EXTENSION, bamFile.getAbsoluteFile().getParentFile());
		file.deleteOnExit(); //TODO We might want to cache this, but we need to figure out how to make sure the files are in sync
		return file;
	}
	
	/**
	 * Iterating over a paired-end BAM file is not simple because the read pairs are not linked
	 * in a way such that one can easily retrieve one given the other. As a work-around, a new
	 * BAM file is written in which each read pair has been combined into a single record spanning the
	 * entire fragment. Information specific to each read (for example, the CIGAR strings) are stored
	 * as custom tags so that the original reads can be reconstructed.
	 * 
	 * @return a reader of the paired-end fragment BAM file associated with this collection. If the file
	 * doesn't yet exist, it is created.
	 */
	private SpecialBAMPECollection getPairedEndFragmentFile() {

		if (this.fragmentReader != null) {
			return fragmentReader;
		/*	this code interferes with multiple calls to numOverlappers
		} else if (this.fragmentFile.exists()) {
			logger.info("Fragment file exists: " + fragmentFile.getName());
			fragmentReader = new SpecialBAMPECollection(fragmentFile);
        	return fragmentReader; */
		} else {
			CloseableIterator<PairedMappedFragment<SAMFragment>> iter = null;
			SAMFileWriter writer = null;
			
			try {
				iter = sortedIterator();
				writer = new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reads.getFileHeader(), false, fragmentFile);

				logger.info("Writing fragment file " + fragmentFile.getName());
				while (iter.hasNext()) {
					PairedMappedFragment<SAMFragment> pair = iter.next();
					writer.addAlignment(convertToCustomSAMFormat(pair));
				}
			} finally {
				if (iter != null)   { iter.close(); }
				if (writer != null) { writer.close(); }
			}

			fragmentReader = new SpecialBAMPECollection(fragmentFile);
			return fragmentReader;
		}
	}
	
	/**
	 * Creates a custom SAMRecord from a pair of reads. The custom SAMRecord will be a
	 * single interval which spans the entire length of the fragment.
	 * @param pair is the paired-end fragment to be converted
	 * @return a SAMRecord in the custom format to be written by a SAMFileWriter
	 */
	private SAMRecord convertToCustomSAMFormat(PairedMappedFragment<SAMFragment> pair) {
		Pair<SAMFragment> minMax = getMinMax(pair);
		SAMRecord alignment = minMax.getValue1().getSamRecord();
		SAMRecord mate = minMax.getValue2().getSamRecord();
		
		SAMRecord record = alignment;
		record.setMateAlignmentStart(alignment.getMateAlignmentStart());
		record.setAlignmentStart(alignment.getAlignmentStart());
		
		// The original CIGAR string can't be retained, since some other programs require the
		// length implied by the CIGAR string to match the fragment length. Store it as a custom
		// tag.
		record.setAttribute(ALIGNMENT_CIGAR, record.getCigarString());
		record.setAttribute(MATE_CIGAR, mate.getCigarString());
		record.setAttribute(MATE_MAPPING_QUALITY, mate.getMappingQuality());
		
		int fragmentLength = Math.max(pair.getRead1().getReferenceEndPosition(), pair.getRead2().getReferenceEndPosition())
				- Math.min(pair.getRead1().getReferenceStartPosition(), pair.getRead2().getReferenceStartPosition());
		
		// Add a dummy CIGAR which is the length of the fragment. Some other programs require the
		// length implied by the CIGAR string to match the fragment length.
	    record.setCigarString(fragmentLength + "M");
		record.setInferredInsertSize(fragmentLength);
		return record;
	}
	
	/**
	 * Get the SAMFragments from a PairedMappedFragment as a Pair. The SAMFragment with
	 * the smaller reference start position as element 1.
	 * @param pair is the PairedMappedFragment
	 * @return the Pair of SAMFragments
	 */
	private Pair<SAMFragment> getMinMax(PairedMappedFragment<SAMFragment> pair) {
		SAMFragment read1 = pair.getRead1();
		SAMFragment read2 = pair.getRead2();
		return pair.getRead1().getReferenceStartPosition() < pair.getRead2().getReferenceStartPosition()
			   ? Pair.of(read1, read2)
		       : Pair.of(read2, read1);
	}
	
	@Override
	public CloseableIterator<PairedMappedFragment<SAMFragment>> sortedIterator() {
		return new FilteredIterator<PairedMappedFragment<SAMFragment>>(new PairedIterator(reads.sortedIterator()), getFilters());
	}
	
	@Override
	public CloseableIterator<PairedMappedFragment<SAMFragment>> sortedIterator(Annotation region, boolean fullyContained) {

		// Get this collection's fragment-file. If it doesn't yet exist, make it.
		SpecialBAMPECollection fragments = this.getPairedEndFragmentFile();

		// Get existing filters and add them.
		Collection<Predicate<PairedMappedFragment<SAMFragment>>> c = getFilters();
		Iterator<Predicate<PairedMappedFragment<SAMFragment>>> iter = c.iterator();
		while (iter.hasNext()) {
			fragments.addFilter(iter.next());
		}
		
		// Add an additional filter to exclude non-overlapping or non-contained fragments.
		if (fullyContained) {
			fragments.addFilter(new ContainedByFilter<PairedMappedFragment<SAMFragment>>(region));
		} else {
			fragments.addFilter(new OverlapsFilter<PairedMappedFragment<SAMFragment>>(region));
		}
		
		return fragments.sortedIterator();
	}
	
	@Override
	public CoordinateSpace getReferenceCoordinateSpace() {
		return reads.getReferenceCoordinateSpace();
	}

	public SAMFileHeader getFileHeader() {
		return reads.getFileHeader();
	}
	
	/**
	 * Writes the reads in this collection to file.
	 * @param file is the file to write the reads to
	 */
	public void writeToFile(File file) {
		writeToFile(file, sortedIterator());
	}
	
	/**
	 * Writes the reads in this collection which overlap the given region to file.
	 * @param file is the file to write the reads to.
	 * @param region is an Annotation which may overlap the reads in this file. Only reads which overlap are output.
	 * @param fullyContained determines whether the reads must be fully contained within the region, or merely overlap it
	 */
	public void writeToFile(File file, Annotation region, boolean fullyContained) {
		if (fullyContained) {
			writeToFile(file, sortedIterator(region, true));
		} else {
			writeToFile(file, sortedIterator(region, false));
		}
	}
	
	private void writeToFile(File file, CloseableIterator<PairedMappedFragment<SAMFragment>> iter) {
		SAMFileWriter writer = null;
		int counter = 0;
		try {
			writer = new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reads.getFileHeader(), false, file);
			while (iter.hasNext()) {
				PairedMappedFragment<SAMFragment> ann = iter.next();
				writer.addAlignment(ann.getRead1().getSamRecord());
				writer.addAlignment(ann.getRead2().getSamRecord());
				if (++counter % 10000 == 0) {
					logger.info("Written " + counter + "fragments to " + file.getName());
				}
			}
		} finally {
			if (writer != null) {writer.close();}
			if (iter != null) {iter.close();}
		}
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
	 * Iterator which goes through a paired-end BAM file and returns complete fragments (that is, fragments
	 * with both reads). It does this by going through the BAM file with a CloseableIterator<SAMFragment>.
	 * When a read is examined, but its mate has yet to be seen, it is stored in a Map and the next read is
	 * examined. If its mate has already been seen, the mate is retrieved from the Map, the read and its mate are
	 * used to construct a PairedMappedFragment, and that PairedMappedFragment is returned.
	 */
	private class PairedIterator implements CloseableIterator<PairedMappedFragment<SAMFragment>> {

		CloseableIterator<SAMFragment> iter;
		Pair<SAMFragment> nextPair;
		Map<String, Pair<SAMFragment>> partials;
		String currentReference;

		public PairedIterator(CloseableIterator<SAMFragment> iter) {
			this.iter = iter;
			this.partials = new TreeMap<String, Pair<SAMFragment>>();
			findNext();
		}

		@Override
		public boolean hasNext() {
			return nextPair != null;
		}

		private void findNext() {
			nextPair = null;
			while (iter.hasNext() && nextPair == null) {			
				SAMFragment read = iter.next();
				SAMRecord rec = read.getSamRecord();
				
				// When switching from chromosomes we should clear cache
				if (!read.getReferenceName().equalsIgnoreCase(currentReference)){
					currentReference = read.getReferenceName();
					this.partials = new TreeMap<String, Pair<SAMFragment>>();
				}
	
				boolean isPaired = read.getSamRecord().getReadPairedFlag();
				boolean mateMapped = !read.getSamRecord().getMateUnmappedFlag();
				boolean onSameReference = rec.getReferenceName().equalsIgnoreCase(rec.getMateReferenceName());
	
				if (isPaired && mateMapped && onSameReference) {
					Pair<SAMFragment> pair = partials.containsKey(read.getName())
											 ? partials.remove(read.getName())       // Mate found in a Pair. Get it.
											 : new Pair<SAMFragment>();              // Mate not seen. Make a new pair.

					if (read.getSamRecord().getFirstOfPairFlag()) {
						pair.setValue1(read);
					} else {
						pair.setValue2(read);
					}
	
					if (pair.hasValue1() && pair.hasValue2()) {
						nextPair = pair;                          // Pair is complete! Return it.
					} else {
						partials.put(read.getName(), pair);       // Pair is incomplete. Store the partial Pair.
					}
				}
			}
		}

		@Override
		public PairedMappedFragment<SAMFragment> next() {
			if (!hasNext()) {
				throw new NoSuchElementException("PairedIterator.next() called with no element.");
			} else {
				PairedMappedFragment<SAMFragment> result = new PairedMappedFragment<SAMFragment>(nextPair);
				findNext();
				return result;
			}
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
	
	/**
	 * Iterating over a paired-end BAM file is not simple because the read pairs are not linked
	 * in a way such that one can easily retrieve one given the other. As a work-around, a new
	 * BAM file is written in which each read pair has been combined into a single record spanning the
	 * entire fragment. Information specific to each read (for example, the CIGAR strings) are stored
	 * as custom tags so that the original reads can be reconstructed. This class represents the
	 * collection stored in the custom BAM file.
	 */
	private class SpecialBAMPECollection extends AbstractAnnotationCollection<PairedMappedFragment<SAMFragment>>{

		private SAMFileReader reader;
		private CoordinateSpace referenceSpace;
		
		public SpecialBAMPECollection(File bamFile){
			super();
			this.reader = new SAMFileReader(bamFile);
			this.referenceSpace = new CoordinateSpace(reader.getFileHeader());
		}

		@Override
		public CloseableIterator<PairedMappedFragment<SAMFragment>> sortedIterator() {
			return new FilteredIterator<PairedMappedFragment<SAMFragment>>(new WrappedIterator(reader), getFilters());
		}

		@Override
		public CloseableIterator<PairedMappedFragment<SAMFragment>> sortedIterator(Annotation region, boolean fullyContained) {
			return new FilteredIterator<PairedMappedFragment<SAMFragment>>(new WrappedIterator(reader, region), getFilters());
		}

		@SuppressWarnings("unused")
		public void writeToFile(String fileName) {
			writeToFile(fileName, sortedIterator());
		}

		@SuppressWarnings("unused")
		public void writeToFile(String fileName, Annotation region) {
			writeToFile(fileName, sortedIterator(region, false));
		}
		
		private void writeToFile(String fileName, CloseableIterator<PairedMappedFragment<SAMFragment>> iter){
			SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(fileName));
		
			while(iter.hasNext()){
				PairedMappedFragment<SAMFragment> ann=iter.next();
				writer.addAlignment(ann.getRead1().getSamRecord());
				writer.addAlignment(ann.getRead2().getSamRecord());
			}
			writer.close();
			iter.close();
		}
		
		@Override
		public CoordinateSpace getReferenceCoordinateSpace() {
			return this.referenceSpace;
		}

		@Override
		public CoordinateSpace getFeatureCoordinateSpace() {
			// FIXME Auto-generated method stub
			throw new UnsupportedOperationException("TODO");
		}
		
		/**
		 * A wrapper class for HTSJDK's CloseableIterator. Iterates over the custom fragment file
		 * and returns the combined record as a PairedMappedFragment<SAMFragment> by splitting the record
		 * into its constituent reads.
		 */
		private class WrappedIterator implements CloseableIterator<PairedMappedFragment<SAMFragment>> {

			SAMRecordIterator iter;
			
			public WrappedIterator(SAMFileReader reader){
				this.iter = reader.iterator();
			}
			
			public WrappedIterator(SAMFileReader reader, Annotation region) {
				this.iter = reader.queryOverlapping(region.getReferenceName(),
													region.getReferenceStartPosition() + 1,
													region.getReferenceEndPosition());
			}

			@Override
			public boolean hasNext() {
				return iter.hasNext();
			}

			@Override
			public PairedMappedFragment<SAMFragment> next() {
				SAMRecord record = iter.next();
				PairedMappedFragment<SAMFragment> reads = getSAMRecords(reader.getFileHeader(), record);
				return reads;
			}

			/**
			 * Split the custom SAM record into its two constituent reads.
			 * @param fileHeader is the file header from the parsed BAM file
			 * @param record is the record to split
			 * @return the paired end fragment with the two reads
			 */
			private PairedMappedFragment<SAMFragment> getSAMRecords(SAMFileHeader fileHeader, SAMRecord record) {

				// Make read 1
				SAMRecord read1 = new SAMRecord(fileHeader);
				read1.setAlignmentStart(record.getAlignmentStart());
				read1.setCigarString(record.getAttribute(ALIGNMENT_CIGAR).toString());
				read1.setFirstOfPairFlag(record.getFirstOfPairFlag());
				read1.setReadName(record.getReadName());
				read1.setReadNegativeStrandFlag(record.getReadNegativeStrandFlag());
				read1.setReferenceName(record.getReferenceName());
				read1.setReadPairedFlag(true);
				read1.setMappingQuality(record.getMappingQuality());
				
				// Make read 2
				SAMRecord read2 = new SAMRecord(fileHeader);
				read2.setAlignmentStart(record.getMateAlignmentStart());
				read2.setCigarString(record.getAttribute(MATE_CIGAR).toString());
				read2.setFirstOfPairFlag(!record.getFirstOfPairFlag());
				read2.setReadName(record.getReadName());
				read2.setReadNegativeStrandFlag(record.getMateNegativeStrandFlag());
				read2.setReferenceName(record.getMateReferenceName());
				read2.setReadPairedFlag(true);
				read2.setMappingQuality(new Integer(record.getAttribute(MATE_MAPPING_QUALITY).toString()));
				
				// Add mate info to read 1
				read1.setMateAlignmentStart(record.getMateAlignmentStart());
				read1.setMateNegativeStrandFlag(record.getMateNegativeStrandFlag());
				read1.setMateReferenceName(record.getMateReferenceName());
				read1.setMateUnmappedFlag(record.getMateUnmappedFlag());
				read1.setProperPairFlag(record.getProperPairFlag());
				
				// Add mate info to read 2
				read2.setMateAlignmentStart(record.getAlignmentStart());
				read2.setMateNegativeStrandFlag(record.getReadNegativeStrandFlag());
				read2.setMateReferenceName(record.getReferenceName());
				read2.setMateUnmappedFlag(record.getReadUnmappedFlag());
				read2.setProperPairFlag(record.getProperPairFlag());
				
				SAMFragment frag1 = new SAMFragment(read1);
				SAMFragment frag2 = new SAMFragment(read2);
				
				if (!record.getFirstOfPairFlag()) {
					return new PairedMappedFragment<SAMFragment>(frag1, frag2);
				}
				
				return new PairedMappedFragment<SAMFragment>(frag1, frag2);	
			}

			@Override
			public void remove() { iter.remove(); }

			@Override
			public void close() { iter.close(); }
		}
	}
	
	/*public BAMPairedFragmentCollection convert(AnnotationCollection<? extends Annotation> features, boolean fullyContained){
		//TODO This needs to be rewritten directly use the paired end iterator to write to disk
		return new BAMPairedFragmentCollection(this.reads.convert(features, fullyContained));
	} */
}

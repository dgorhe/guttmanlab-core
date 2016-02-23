package guttmanlab.core.annotationcollection;

import guttmanlab.core.annotation.MappedFragment;

import java.io.File;
import java.io.IOException;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

/**
 * Get bam collection as single or paired end by checking the file
 * @author prussell
 *
 */
public class BAMFragmentCollectionFactory {
	
	private static Logger logger = Logger.getLogger(BAMFragmentCollectionFactory.class.getName());
	
	/**
	 * @param bamFile Bam file
	 * @return True iff the file contains paired reads
	 */
	public static boolean isPairedEnd(String bamFile) {
		return isPairedEnd(new File(bamFile));
	}
	
	/**
	 * @param bamFile Bam file
	 * @return True iff the file contains paired reads
	 */
	public static boolean isPairedEnd(File bamFile) {
		SAMFileReader reader = new SAMFileReader(bamFile);
		SAMRecordIterator iter = reader.iterator();
		if(!iter.hasNext()) {
			reader.close();
			throw new IllegalArgumentException("Bam file is empty");
		}
		SAMRecord firstRecord = iter.next();
		reader.close();
		return firstRecord.getReadPairedFlag();
	}
	
	/**
	 * Get an appropriate implementation of annotation collection depending on whether the reads are paired
	 * @param bamFile Bam file
	 * @return Single or paired end read collection
	 * @throws IOException 
	 */
	public static AbstractAnnotationCollection<? extends MappedFragment> createFromBam(String bamFile) throws IOException {
		return createFromBam(new File(bamFile));
	}
	
	/**
	 * Get an appropriate implementation of annotation collection depending on whether the reads are paired
	 * @param bamFile Bam file
	 * @return Single or paired end read collection
	 * @throws IOException 
	 */
	public static AbstractAnnotationCollection<? extends MappedFragment> createFromBam(File bamFile) throws IOException {
		return createFromBam(bamFile, false);
	}
	
	/**
	 * Get an appropriate implementation of annotation collection depending on whether the reads are paired
	 * Optionally force single end implementation
	 * @param bamFile Bam file
	 * @param forceSingleEnd Treat reads as single end even if they are paired
	 * @return Single or paired end read collection
	 * @throws IOException 
	 */
	public static AbstractAnnotationCollection<? extends MappedFragment> createFromBam(String bamFile, boolean forceSingleEnd) throws IOException {
		return createFromBam(new File(bamFile), forceSingleEnd);
	}
	
	/**
	 * Get an appropriate implementation of annotation collection depending on whether the reads are paired
	 * Optionally force single end implementation
	 * @param bamFile Bam file
	 * @param forceSingleEnd Treat reads as single end even if they are paired
	 * @return Single or paired end read collection
	 * @throws IOException 
	 */
	public static AbstractAnnotationCollection<? extends MappedFragment> createFromBam(File bamFile, boolean forceSingleEnd) throws IOException {
		
		if(forceSingleEnd) {
			logger.info("Reading bam file " + bamFile.getName() + " and forcing single end implementation.");
			return new BAMSingleReadCollection(bamFile);
		}
		
		if(isPairedEnd(bamFile)) {
			logger.info("Reading bam file " + bamFile.getName() + " as paired ends.");
			return new BAMPairedFragmentCollection(bamFile);
		}
		
		logger.info("Reading bam file " + bamFile.getName() + " as single ends.");
		return new BAMSingleReadCollection(bamFile);
		
	}
	
}

package guttmanlab.core.pipeline.util;

import guttmanlab.core.pipeline.LSFJob;
import guttmanlab.core.pipeline.OGSJob;
import guttmanlab.core.pipeline.Scheduler;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.Session;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BamUtils {

	private static Logger logger = Logger.getLogger(BamUtils.class.getName());

	/**
	 * Merge bam files
	 * @param inputBams Names of bam files to merge
	 * @param outputMergedBam Output merged bam file to write
	 * @param scheduler Scheduler
	 * @param drmaaSession DRMAA session if using OGS or null otherwise
	 * @param picardJarDir Directory containing Picard jar files
	 * @param assumeSorted Assume input files are sorted
	 * @throws DrmaaException
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void mergeBams(Collection<String> inputBams, String outputMergedBam, Scheduler scheduler, Session drmaaSession, String picardJarDir, boolean assumeSorted) throws DrmaaException, IOException, InterruptedException {
		logger.info("Creating merged bam file " + outputMergedBam + ".");
		String inputs = "";
		for(String bam : inputBams) {
			inputs += "INPUT=" + bam + " ";
		}
		String output = "OUTPUT=" + outputMergedBam;
		String cmmd = "java -jar " + picardJarDir + "/MergeSamFiles.jar " + inputs + " " + output + " ASSUME_SORTED=" + assumeSorted + " MERGE_SEQUENCE_DICTIONARIES=true";
		logger.info("Running picard command: " + cmmd);
		switch(scheduler) {
		case LSF:
			String jobID = Long.valueOf(System.currentTimeMillis()).toString();
			logger.info("LSF job ID is " + jobID + ".");
			// Submit job
			LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, "merge_bam_files_" + jobID + ".bsub", "week", 8);	
			job.submit();
			job.waitFor();
			break;
        case OGS:
            if(drmaaSession == null) {
                    throw new IllegalArgumentException("DRMAA session is null. Must provide an active DRMAA session to use OGS. There can only be one active session at a time. Session should have been created in the main method of the class calling this method.");
            }
            OGSJob ogsJob = new OGSJob(drmaaSession, cmmd, "merge_bam_files", null);
            ogsJob.submit();
            logger.info("OGS job ID is " + ogsJob.getID() + ".");
            ogsJob.waitFor();
            break;
		default:
			throw new IllegalArgumentException("Scheduler " + scheduler.toString() + " is not supported.");
		}

	}

	
	/**
	 * Split a bam file into several smaller bam files
	 * @param inputBam Bam file to split
	 * @param numFilesToWrite Number of smaller files to write
	 * @return The names of the smaller files
	 */
	public static Collection<String> splitBam(String inputBam, int numFilesToWrite) {
		return splitBam(inputBam, numFilesToWrite, false);
	}
	
	/**
	 * Split a bam file into several smaller bam files
	 * @param inputBam Bam file to split
	 * @param numFilesToWrite Number of smaller files to write
	 * @param getNamesOnly Only get split file names, do not actually write files
	 * @return The names of the smaller files
	 */
	public static Collection<String> splitBam(String inputBam, int numFilesToWrite, boolean getNamesOnly) {
		

		if(numFilesToWrite < 1) {
			throw new IllegalArgumentException("Must write at least 1 file");
		}
		
		logger.info("");
		logger.info("Splitting " + inputBam + " into " + numFilesToWrite + " files.");
		
		SAMFileReader reader = new net.sf.samtools.SAMFileReader(new File(inputBam));
		SAMFileHeader header = reader.getFileHeader();
		
		Collection<String> rtrn = new ArrayList<String>();
		
		ArrayList<SAMFileWriter> writers = new ArrayList<SAMFileWriter>();
		for(int i = 0; i < numFilesToWrite; i++) {
			String name = inputBam + "." + i;
			rtrn.add(name);
			if(!getNamesOnly) {
				SAMFileWriterFactory factory = new SAMFileWriterFactory();
				SAMFileWriter writer = factory.makeBAMWriter(header, false, new File(name));
				writers.add(writer);
			}
		}
		
		if(getNamesOnly) {
			reader.close();
			return rtrn;
		}
		
		int numDone = 0;
		
		SAMRecordIterator iter = reader.iterator();
		while(iter.hasNext()) {
			SAMFileWriter writer = writers.get(numDone % writers.size());
			SAMRecord record = iter.next();
			writer.addAlignment(record);
			numDone++;
			if(numDone % 1000000 == 0) {
				logger.info("Finished " + numDone + " records.");
			}
		}
		
		reader.close();
		
		for(SAMFileWriter writer : writers) {
			writer.close();
		}
		
		logger.info("Done splitting bam file.");
		return rtrn;
		
	}
	
	/**
	 * Convert SAM file to BAM file
	 * @param inputSam SAM file to convert
	 * @param outputBam Output BAM file
	 */
	public static void samToBam(File inputSam, File outputBam) {
		
		logger.info("Converting SAM to BAM");
		SamReader samReader = SamReaderFactory.makeDefault().open(inputSam);
		htsjdk.samtools.SAMFileHeader header = samReader.getFileHeader();
		htsjdk.samtools.SAMFileWriter writer = new htsjdk.samtools.SAMFileWriterFactory().makeBAMWriter(header, false, outputBam);
		logger.info("Reading from " + inputSam.getAbsolutePath() + " and writing to " + outputBam.getAbsolutePath());
		htsjdk.samtools.SAMRecordIterator iter = samReader.iterator();
		while(iter.hasNext()) {
			htsjdk.samtools.SAMRecord record = iter.next();
			writer.addAlignment(record);
		}
		writer.close();
		iter.close();
		logger.info("Done writing bam.");
		
	}
	
	/**
	 * Sort a BAM file
	 * @param inputBam Input BAM
	 * @param outputSorted Output sorted BAM
	 */
	public static void sortBam(File inputBam, File outputSorted) {
		
		logger.info("Sorting BAM file " + inputBam.getAbsolutePath());
		String[] argv = {"I=" + inputBam.getAbsolutePath(), "O=" + outputSorted.getAbsolutePath(), "SORT_ORDER=coordinate"};
		new net.sf.picard.sam.SortSam().instanceMain(argv);
		logger.info("Done sorting BAM file");
		
	}
	
	
	/**
	 * Index BAM file
	 * @param inputBam Input BAM
	 */
	public static void indexBam(File inputBam) {
		
		logger.info("Indexing BAM file " + inputBam.getAbsolutePath());
		String[] argv = {"I=" + inputBam.getAbsolutePath()};
		new net.sf.picard.sam.BuildBamIndex().instanceMain(argv);
		logger.info("Done indexing BAM file");
		
	}
	
	
	

}

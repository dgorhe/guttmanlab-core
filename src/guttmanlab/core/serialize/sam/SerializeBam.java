package guttmanlab.core.serialize.sam;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;

import org.apache.avro.AvroRuntimeException;
import org.apache.avro.Schema;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumWriter;
import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecord.SAMTagAndValue;
import net.sf.samtools.SAMRecordIterator;

/**
 * Convert a bam file to an avro file
 * @author prussell
 *
 */
public class SerializeBam {
	
	private static Logger logger = Logger.getLogger(SerializeBam.class.getName());
	
	/**
	 * Get whether a record passes all of a collection of filters
	 * @param record The record
	 * @param filters The filters
	 * @return True if the record passes all filters or the collection of filters is empty
	 */
	private static <T extends Object> boolean passesAll(T record, Collection<Predicate<T>> filters) {
		return filters.stream().allMatch(predicate -> predicate.test(record));
	}
	
	
	
	/**
	 * Write the avro file
	 * @param schemaFile Avro schema file with .avsc extension
	 * @param inputBam Bam file to serialize
	 * @param outputAvro Avro file to write
	 * @throws IOException
	 */
	public static void serialize(String schemaFile, String inputBam, String outputAvro, Collection<Predicate<SAMRecord>> filters) throws IOException {
		
		logger.info("Serializing " + inputBam + "...");
		
		// Reader for bam file
		SAMFileReader samReader = new SAMFileReader(new File(inputBam));
		SAMRecordIterator samIter = samReader.iterator();
		// Create the schema
		Schema schema = new Schema.Parser().parse(new File(schemaFile));
		// This file will have Avro output data
		File AvroFile = new File(outputAvro);
		// Create a writer to serialize the record
		DatumWriter<GenericRecord> datumWriter = new GenericDatumWriter<GenericRecord>(schema);		         
		DataFileWriter<GenericRecord> dataFileWriter = new DataFileWriter<GenericRecord>(datumWriter);
		dataFileWriter.setCodec(CodecFactory.snappyCodec());
		dataFileWriter.create(schema, AvroFile);
		
		int numDone = 0;		
		// Iterate over bam file and write to Avro output file
		while(samIter.hasNext()) {
			numDone++;
			if(numDone % 100000 == 0) logger.info("Finished " + numDone + " records");
			SAMRecord samRecord = samIter.next();
			if(!passesAll(samRecord, filters)) continue;
			// Create a record to hold sam record
			GenericRecord avroRec = new GenericData.Record(schema);
			avroRec.put("qname", samRecord.getReadName());
			avroRec.put("flag", samRecord.getFlags());
			avroRec.put("rname", samRecord.getReferenceName());
			avroRec.put("pos", samRecord.getAlignmentStart());
			avroRec.put("mapq", samRecord.getMappingQuality());
			avroRec.put("cigar", samRecord.getCigarString());
			avroRec.put("rnext", samRecord.getMateReferenceName());
			avroRec.put("pnext", samRecord.getMateAlignmentStart());
			avroRec.put("tlen", samRecord.getInferredInsertSize());
			avroRec.put("seq", samRecord.getReadString());
			avroRec.put("qual", samRecord.getBaseQualityString());
			List<SAMTagAndValue> tags = samRecord.getAttributes();
			for(SAMTagAndValue tag : tags) {
				String name = "tag" + tag.tag;
				try {
					avroRec.put(name, tag.value);
				} catch(AvroRuntimeException e) {
					samReader.close();
					dataFileWriter.close();
					e.printStackTrace();
					throw new IllegalStateException("Error probably caused by the fact that schema does not contain tag " + tag.tag);
				}
			}
			dataFileWriter.append(avroRec);
		}  // end of for loop

		samReader.close();
		dataFileWriter.close();
		
	}
	
	
}

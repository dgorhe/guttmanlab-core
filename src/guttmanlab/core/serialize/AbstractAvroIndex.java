package guttmanlab.core.serialize;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.TreeMap;

import org.apache.avro.Schema;
import org.apache.avro.file.DataFileReader;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumReader;
import org.apache.log4j.Logger;

/**
 * Implementations of methods that can tolerate the generic type
 * @author prussell
 *
 * @param <T>
 */
public abstract class AbstractAvroIndex<T extends Comparable<T>> implements AvroIndex<T> {
	
	private String avroFileName; // The avro data file based on the schema
	protected String indexFileName; // The index written by this class
	protected TreeMap<T, Long> positionsByKey; // File position of record at beginning of each block
	protected DataFileReader<GenericRecord> reader; // Reader for avro file
	private Schema schema; // Data schema
	protected String indexedField; // The name of the field in the schema that is indexed by this index
	private static Logger logger = Logger.getLogger(AbstractAvroIndex.class.getName());
	private static String VALIDATED_INDEX_SUFFIX = ".VALIDATED";
	
	/**
	 * Get index file name corresponding to an avro file name
	 * @param avroFile Avro file name
	 * @return Index file name
	 */
	public static String getIndexFileName(String avroFile) {
		return avroFile + ".index";
	}
	
	/**
	 * Read index information from index file
	 * @param avroFile Avro file
	 * @param schemaFile Avro schema file
	 * @param indexedFieldName The name of the field in the schema that is indexed in this index
	 */
	public AbstractAvroIndex(String avroFile, String schemaFile, String indexedFieldName) throws IOException {
		logger.info("Loading avro index...");
		avroFileName = avroFile;
		indexFileName = getIndexFileName(avroFileName);		
		File f = new File(indexFileName);
		if(!f.exists()) {
			throw new IllegalStateException("No index exists for avro file " + avroFileName + ".");
		}
		schema = new Schema.Parser().parse(new File(schemaFile));
		indexedField = indexedFieldName;
		//genericRecord = new GenericData.Record(schema);
		DatumReader<GenericRecord> datumReader = new GenericDatumReader<GenericRecord>(schema);
		reader = new DataFileReader<GenericRecord>(new File(avroFileName), datumReader);
		loadIndex(!indexIsValidated());
		logger.info("Done loading and validating index.");
	}
	
	/**
	 * Check if there is a validation file with the index file that is newer than the index
	 * @return True iff the index has been validated
	 * @throws IOException 
	 */
	private boolean indexIsValidated() throws IOException {
		File avroFile = new File(avroFileName).getCanonicalFile();
		File indexFile = new File(indexFileName).getCanonicalFile();
		// Check that index is newer than database
		if(indexFile.lastModified() < avroFile.lastModified()) {
			throw new IllegalStateException("Database " + avroFile.getAbsolutePath() + " was modified after index " + indexFile.getAbsolutePath());
		}
		// Check that validation file exists
		File validatedFile = new File(getValidationFileName());
		if(!validatedFile.exists()) return false;
		// Check that validation file is newer than index
		if(indexFile.lastModified() > validatedFile.lastModified()) {
			validatedFile.delete();
			return false;
		}
		logger.info("Index " + indexFileName + " is validated.");
		return true;
	}

	/**
	 * Write the validation file for the index
	 * @throws IOException
	 */
	public void reportIndexValidated() throws IOException {
		File validatedFile = new File(getValidationFileName());
		@SuppressWarnings("unused")
		boolean created = validatedFile.createNewFile();
		validatedFile.setLastModified(System.currentTimeMillis());
	}
	
	/**
	 * Get validation file name
	 * @return Validation file name
	 * @throws IOException
	 */
	private String getValidationFileName() throws IOException {
		return new File(indexFileName).getCanonicalPath() + VALIDATED_INDEX_SUFFIX;
	}
	
	/**
	 * Get a lower bound on the file position of the first occurrence of a key
	 * @param key
	 * @return The start position of the first block containing the key,
	 * or if the first record containing the key happens to be the first record of a block,
	 * this method will return the beginning of the previous block
	 */
	private long getLowerBoundStartPos(T key) {
		// The file position of the last key before the key of interest
		// Head map is exclusive in case the key extends over multiple blocks
		try {
			long rtrn = positionsByKey.headMap(key, false).lastEntry().getValue().longValue();
			return rtrn;
		} catch(NullPointerException e) {
			if(positionsByKey.firstEntry().getKey().equals(key)) {
				throw new IllegalStateException("Query does not work for the very first barcode in file: " + key + ". TODO: fix.");
			}
			throw e;
		}
	}
	
	/**
	 * Move the data reader to the beginning of the first block containing the key,
	 * or if the first record containing the key happens to be the first record of a block,
	 * this method will move the reader to the beginning of the previous block
	 * @param <D>
	 * @param key The key
	 * @throws IOException 
	 */
	protected void seekToBlockStart(T key) throws IOException {
		reader.seek(getLowerBoundStartPos(key));
	}
	
	/**
	 * @return The current position of the reader
	 * @throws IOException 
	 */
	public long getCurrentFilePosition() {
		try {
			return reader.tell();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
			return -1; // Unreachable
		}
	}
	
	/**
	 * Write a file index with respect to the specified field
	 * @param fieldName Field to index
	 * @param schemaFile Avro schema file
	 * @param inputAvro Avro file to index
	 * @param output Index file to write
	 * @throws IOException
	 */
	public static void writeIndex(String fieldName, String schemaFile, String inputAvro) throws IOException {
		logger.info("");
		logger.info("Writing index for file " + inputAvro + " based on field " + fieldName + "...");
		long start = System.nanoTime();
		// Create a schema
		Schema schema = new Schema.Parser().parse(new File(schemaFile));
		// Create a record using schema
		GenericRecord avroRec = new GenericData.Record(schema);
		DatumReader<GenericRecord> datumReader = new GenericDatumReader<GenericRecord>(schema);
		DataFileReader<GenericRecord> dataFileReader = new DataFileReader<GenericRecord>(new File(inputAvro), datumReader);
		// Write index to a file
		FileWriter w = new FileWriter(getIndexFileName(inputAvro));
		// Go to the first sync point
		dataFileReader.previousSync();
		String prevName = "";
		int numDone = 0;
		while (true) {
			long pos = dataFileReader.tell();
			numDone++;
			if(numDone % 1000 == 0) {
				logger.info("Finished " + numDone + " sync points. Current position: " + pos);
			}
			try {
				avroRec = dataFileReader.next();
			} catch(NoSuchElementException e) {
				break;
			}
			Object val = avroRec.get(fieldName);
			if(val == null) {
				w.close();
				dataFileReader.close();
				throw new IllegalStateException("Indexed value cannot be null");
			}
			String name = val.toString();
			if(name.compareTo(prevName) < 0) {
				w.close();
				dataFileReader.close();
				throw new IllegalStateException("Values out of order: " + prevName + " " + name);
			}
			w.write(name + "\t" + pos + "\n");
			dataFileReader.sync(pos);
			prevName = name;
		}
		w.close();
		long time = (System.nanoTime() - start) / 1000000000;
		logger.info("Wrote index in " + time + " seconds.");
		dataFileReader.close();
	}
	
	
}

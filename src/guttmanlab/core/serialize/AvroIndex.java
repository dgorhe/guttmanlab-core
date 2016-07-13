package guttmanlab.core.serialize;

import java.io.IOException;
import java.util.List;

import org.apache.avro.generic.GenericRecord;

/**
 * An index object for an avro file
 * Avro file should be sorted by the field that is indexed
 * @author prussell
 *
 * @param <T> The type of the indexed field
 */
public interface AvroIndex<T extends Comparable<T>> {
	
	/**
	 * Load the index into memory
	 * @param validate While loading, check that the index is valid, i.e., records are sorted by the proper key and file positions match the avro file
	 * @throws IOException
	 */
	public void loadIndex(boolean validate);
		
	/**
	 * @return The current position of the reader
	 * @throws IOException 
	 */
	public long getCurrentFilePosition();
	
	/**
	 * Seek to the first record with this key in the avro file
	 * @param key The key
	 * @return The first record with this key, which should be the current record for the data reader
	 * @throws IOException
	 */
	public <P extends GenericRecord> P seek(T key);
	
	/**
	 * Get an ordered list of all records with this key
	 * @param key The key
	 * @return List of records
	 * @throws IOException 
	 */
	public <P extends GenericRecord> List<P> get(T key);

}

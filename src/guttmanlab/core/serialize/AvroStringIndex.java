package guttmanlab.core.serialize;

import guttmanlab.core.util.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.avro.AvroRuntimeException;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericData.Record;
import org.apache.avro.generic.GenericRecord;
import org.apache.log4j.Logger;


/**
 * An index for a string field
 * @author prussell
 *
 */
public class AvroStringIndex extends AbstractAvroIndex<String> {

	private static Logger logger = Logger.getLogger(AvroStringIndex.class.getName());
	public static long MAX_RECORDS_TO_GET = Long.MAX_VALUE; //If there are more than this many records with the key, throw an IllegalArgumentException
	
	/**
	 * Read index information from index file
	 * @param avroFileName Avro file
	 * @param schemaFile Avro schema file
	 * @param indexedFieldName The name of the field in the schema that is indexed in this index
	 */
	public AvroStringIndex(String avroFileName, String schemaFile, String indexedFieldName) throws IOException {
		super(avroFileName, schemaFile, indexedFieldName);
	}

	@Override
	public void loadIndex(boolean validate) {
		logger.info("");
		logger.info("Reading index from file " + indexFileName + "...");
		try {
			positionsByKey = new TreeMap<String, Long>();
			FileReader r = new FileReader(indexFileName);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			long previousPos = 0;
			long currentPos = 0;
			String[] previousKey = new String[1];
			Collection<Long> messagesWritten = new TreeSet<Long>();
			while(b.ready()) {
				s.parse(b.readLine());
				String key = s.asString(0);
				currentPos = s.asLong(1);
				// Check that the block numbers are increasing
				if(validate) {
					if(currentPos <= previousPos) {
						b.close();
						throw new IllegalStateException("Index validation error: wrong order for file blocks " + previousPos + ", " + currentPos + ".");
					}
				}
				if(!key.equals(previousKey[0])) {
					// We've found a new key
					previousKey[0] = key;
					if(validate) {
						// Check that the first record in this block has the reported key
						GenericData.Record record = null;
						reader.sync(previousPos);
						try {
							record = new GenericData.Record((Record) reader.next(), true);
						} catch(AvroRuntimeException e) {
							e.printStackTrace();
							logger.warn("Caught exception. Skipping record at position " + currentPos);
							continue;
						}
						String recordKey = record.get(indexedField).toString();
						if(!recordKey.equals(key)) {
							b.close();
							throw new IllegalStateException("At file block " + reader.tell() + ": Index validation error: beginning of block " + currentPos + " in avro file has key " + recordKey + ". Index has " + key + ".");
						}
					}
					Long prevVal = positionsByKey.put(key, Long.valueOf(currentPos));
					if(validate) {
						if(prevVal != null) {
							b.close();
							throw new IllegalStateException("Key " + key + " was already in map with value " + prevVal);
						}
						long nearestMillion = currentPos - (currentPos % 1000000);
						if(nearestMillion % 10000000 == 0 && !messagesWritten.contains(Long.valueOf(nearestMillion))) {
							logger.info("Index position " + nearestMillion);
							messagesWritten.add(Long.valueOf(nearestMillion));
						}				
					}
				}
				previousPos = currentPos;
			}
			b.close();
			if(validate) {
				reportIndexValidated();
			}
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	

	@SuppressWarnings("unchecked")
	@Override
	public GenericRecord seek(String key) {
		try {
			seekToBlockStart(key);
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}
		while(true) {
			try {
				GenericData.Record record = new GenericData.Record((Record) reader.next(), true);
				Object currKey = record.get(indexedField);
				int compare = 0;
				String strCurrKey = currKey.toString();
				compare = strCurrKey.compareTo(key);
				//logger.info(key + "\t" + strCurrKey);
				if(compare == 0) {
					return record;
				}
				if(compare > 0) {
					throw new NoSuchElementException("Key " + key + " not found.");
				}
			} catch(NoSuchElementException e) {
				break;
			}
		}
		throw new NoSuchElementException("Key " + key + " not found.");
	}

	@Override
	public List<GenericRecord> get(String key) {
		try {
			return get(key, null, null);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
			return null; // Unreachable
		}
	}
	
	/**
	 * Get all records with the key
	 * Optionally exclude records with some attribute contained in a set of values to exclude
	 * @param key The key to search for
	 * @param nameOfAttributeForExclusionSet The name of the attribute to check for exclusion, or null if not using
	 * @param attributeValuesToExclude The attribute values to exclude, or null if not using
	 * @return The set of records with the desired key, minus records with the other attribute contained in the exclusion set
	 * @throws IOException
	 */
	public List<GenericRecord> get(String key, String nameOfAttributeForExclusionSet, Collection<String> attributeValuesToExclude) throws IOException {
		List<GenericRecord> rtrn = new ArrayList<GenericRecord>();
		try {
			GenericData.Record first = new GenericData.Record((Record) seek(key), true);
			if(nameOfAttributeForExclusionSet != null && attributeValuesToExclude != null) {
				String attribute = first.get(nameOfAttributeForExclusionSet).toString();
				if(!attributeValuesToExclude.contains(attribute)) {
					rtrn.add(first);
				}
			} else {
				rtrn.add(first);
			}
			while(true) {
				try {
					GenericData.Record record = new GenericData.Record((Record) reader.next(), true);
					Object currKey = record.get(indexedField);
					int compare = 0;
					String strCurrKey = currKey.toString();
					compare = strCurrKey.compareTo(key);
					if(compare == 0) {
						if(nameOfAttributeForExclusionSet != null && attributeValuesToExclude != null) {
							String attribute = record.get(nameOfAttributeForExclusionSet).toString();
							if(!attributeValuesToExclude.contains(attribute)) {
								rtrn.add(record);
							}
						} else {
							rtrn.add(record);
						}
					} else if(compare > 0) {
						return rtrn;
					} else {
						throw new IllegalStateException("Wrong sort order for keys " + key + " " + strCurrKey);
					}
				} catch(NoSuchElementException e) {
					break;
				}
			}
		} catch(IllegalStateException e) {
			logger.warn("Caught exception on query " + key);
			logger.warn("Returned matches will be incomplete for this key.");
			logger.warn(e.getMessage());
		}
		if(rtrn.size() <= MAX_RECORDS_TO_GET) {
			return rtrn;
		}
		throw new IllegalArgumentException("Key " + key + " has " + rtrn.size() + " records, more than the max of " + MAX_RECORDS_TO_GET + ".");
	}
	
	
}

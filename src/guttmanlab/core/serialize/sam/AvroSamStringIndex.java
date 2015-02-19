package guttmanlab.core.serialize.sam;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.serialize.AvroIndex;
import guttmanlab.core.serialize.AvroStringIndex;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.avro.generic.GenericRecord;

public class AvroSamStringIndex implements AvroIndex<String> {

	private AvroStringIndex stringIndex;
	private AnnotationCollection<Gene> excludeRegions;
	
	/**
	 * @param avroFileName Avro database file
	 * @param schemaFile Avro schema file
	 * @param indexedFieldName Name of indexed field
	 * @throws IOException
	 */
	public AvroSamStringIndex(String avroFileName, String schemaFile, String indexedFieldName) throws IOException {
		this(avroFileName, schemaFile, indexedFieldName, null);
	}
	
	/**
	 * @param avroFileName Avro database file
	 * @param schemaFile Avro schema file
	 * @param indexedFieldName Name of indexed field
	 * @param regionsToExclude Exclude matches that overlap any of these annotations, or null if not using
	 * @throws IOException
	 */
	public AvroSamStringIndex(String avroFileName, String schemaFile, String indexedFieldName, AnnotationCollection<Gene> regionsToExclude) throws IOException {
		stringIndex = new AvroStringIndex(avroFileName, schemaFile, indexedFieldName);
		excludeRegions = regionsToExclude;
	}

	@Override
	public void loadAndValidateIndex() throws IOException {
		stringIndex.loadAndValidateIndex();
	}

	@SuppressWarnings("unchecked")
	@Override
	public AvroSamRecord seek(String key) throws IOException {
		//GenericRecord genericRecord = stringIndex.seek(key);
		//return new AvroSamRecord(genericRecord);
		throw new UnsupportedOperationException();
	}

	@SuppressWarnings("unchecked")
	@Override
	public List<AvroSamRecord> get(String key) throws IOException {
		List<GenericRecord> genericRecords = stringIndex.get(key);
		List<AvroSamRecord> rtrn = new ArrayList<AvroSamRecord>();
		for(GenericRecord record : genericRecords) {
			try {
				AvroSamRecord samRecord = new AvroSamRecord(record);
				// Skip the record if it overlaps a region from the exclusion set
				if(excludeRegions != null && excludeRegions.overlaps(samRecord)) {
					continue;
				}
				rtrn.add(samRecord);
			} catch(IllegalStateException e) {
				// Invalid mapping quality
				continue;
			}
		}
		return rtrn;
	}

	@Override
	public long getCurrentFilePosition() throws IOException {
		return stringIndex.getCurrentFilePosition();
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
	public List<AvroSamRecord> get(String key, String nameOfAttributeForExclusionSet, Collection<String> attributeValuesToExclude) throws IOException {
		List<GenericRecord> genericRecords = stringIndex.get(key, nameOfAttributeForExclusionSet, attributeValuesToExclude);
		List<AvroSamRecord> rtrn = new ArrayList<AvroSamRecord>();
		for(GenericRecord record : genericRecords) {
			try {
				AvroSamRecord samRecord = new AvroSamRecord(record);
				// Skip the record if it overlaps a region from the exclusion set
				if(excludeRegions != null && excludeRegions.overlaps(samRecord)) {
					continue;
				}
				rtrn.add(samRecord);
			} catch(IllegalStateException e) {
				// Invalid mapping quality
				continue;
			}
		}
		return rtrn;
	}

	/**
	 * Get all records with the key as an annotation collection
	 * @param key The key to search for
	 * @return The set of records with the desired key
	 * @throws IOException
	 */
	public FeatureCollection<AvroSamRecord> getAsAnnotationCollection(String key) throws IOException {
		return fromList(get(key));
	}

	/**
	 * Get all records with the key as an annotation collection
	 * Optionally exclude records with some attribute contained in a set of values to exclude
	 * @param key The key to search for
	 * @param nameOfAttributeForExclusionSet The name of the attribute to check for exclusion, or null if not using
	 * @param attributeValuesToExclude The attribute values to exclude, or null if not using
	 * @return The set of records with the desired key, minus records with the other attribute contained in the exclusion set
	 * @throws IOException
	 */
	public FeatureCollection<AvroSamRecord> getAsAnnotationCollection(String key, String nameOfAttributeForExclusionSet, Collection<String> attributeValuesToExclude) throws IOException {
		return fromList(get(key, nameOfAttributeForExclusionSet, attributeValuesToExclude));
	}
	
	private static FeatureCollection<AvroSamRecord> fromList(List<AvroSamRecord> list) {
		FeatureCollection<AvroSamRecord> rtrn = new FeatureCollection<AvroSamRecord>(null);
		for(AvroSamRecord record : list) {
			// Skip the record if its mapping quality is too low or not available
			if(!record.mappingQualityIsOk()) {
				continue;
			}
			rtrn.add(record);
		}
		return rtrn;
	}
	
}

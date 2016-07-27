package guttmanlab.core.serialize.sam;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.predicate.ReadFlag;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.util.SAMFlagDecoder;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.sf.samtools.SAMFileHeader;

import org.apache.avro.Schema;
import org.apache.avro.Schema.Field;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang.builder.HashCodeBuilder;

public class AvroSamRecord extends BlockedAnnotation implements GenericRecord, MappedFragment {
	
	private GenericRecord record;
	private Annotation annotation;
	private boolean firstReadTranscriptionStrand;
	
	/**
	 * @param genericRecord Record
	 */
	public AvroSamRecord(GenericRecord genericRecord) {
		this(genericRecord, true);
	}
	
	/**
	 * @param genericRecord Record
	 * @param firstReadIsTranscriptionStrand True if read1 is 5' to 3', false if read2 is
	 */
	public AvroSamRecord(GenericRecord genericRecord, boolean firstReadIsTranscriptionStrand) {
		
		// Get basic attributes
		firstReadTranscriptionStrand = firstReadIsTranscriptionStrand;
		record = genericRecord;
		String cigar = getStringAttributeOrThrow("cigar");
		String chr = getReferenceName();
		int start = getReferenceStartPosition();
		String name = getName();
		
		// Determine strand
		int flag = getIntAttributeOrThrow("flag");
		SAMFlagDecoder decoder = new SAMFlagDecoder(flag);
		boolean isPaired = decoder.readPaired();
		boolean isFirst = decoder.firstInPair();
		boolean plusStrand = !decoder.readReverseStrand();
		Strand strand = null;
		if(isPaired) {
			if((firstReadTranscriptionStrand && isFirst) || (!firstReadTranscriptionStrand && !isFirst)) {
				strand = plusStrand ? Strand.POSITIVE : Strand.NEGATIVE;
			}
			if((firstReadTranscriptionStrand && !isFirst) || (!firstReadTranscriptionStrand && isFirst)) {
				strand = plusStrand ? Strand.NEGATIVE : Strand.POSITIVE;
			}
		} else {
			if(firstReadTranscriptionStrand && plusStrand || !firstReadTranscriptionStrand && !plusStrand) {
				strand = Strand.POSITIVE;
			}
			if(firstReadTranscriptionStrand && !plusStrand || !firstReadTranscriptionStrand && plusStrand) {
				strand = Strand.NEGATIVE;
			}
		}
		
		// Construct annotation
		annotation = SAMFragment.parseCigar(cigar, chr, start, strand, name);
	}
	
	/**
	 * Get the value of a SAM tag
	 * @param attributeName Tag name
	 * @return The value of the SAM tag
	 */
	public Object getAttribute(String attributeName) {
		return record.get(attributeName);
	}
	
	/**
	 * Get the value of a string SAM tag or throw an exception if tag is absent
	 * @param attributeName Tag name
	 * @return The string value
	 */
	public String getStringAttributeOrThrow(String attributeName) {
		return getStringAttribute(attributeName).orElseThrow(() -> new IllegalArgumentException("Record does not contain attribute: " + attributeName));
	}
	
	/**
	 * Get the value of a string SAM tag
	 * @param attributeName Tag name
	 * @return The string value or empty if tag is absent
	 */
	public Optional<String> getStringAttribute(String attributeName) {
		try {
			return Optional.of(record.get(attributeName).toString());
		} catch(NullPointerException e) {
			return Optional.empty();
		}
	}
	
	
	/**
	 * Get the value of a SAM integer tag or throw an exception if tag is absent
	 * @param attributeName Tag name
	 * @return Integer value of the tag
	 */
	public int getIntAttributeOrThrow(String attributeName) {
		return getIntAttribute(attributeName).orElseThrow(() -> new IllegalArgumentException("Record does not contain attribute: " + attributeName)).intValue();
	}
	
	/**
	 * Get the value of a SAM integer tag
	 * @param attributeName Tag name
	 * @return Integer value of the tag or empty if tag is absent
	 */
	public Optional<Integer> getIntAttribute(String attributeName) {
		try {
			return Optional.of(Integer.valueOf((int) record.get(attributeName)));
		} catch(NullPointerException e) {
			return Optional.empty();
		}
	}

	@Override
	public Object get(int attributeName) {
		return record.get(attributeName);
	}

	@Override
	public void put(int attributeName, Object value) {
		throw new UnsupportedOperationException("AvroSamRecord objects are immutable");
	}

	@Override
	public Schema getSchema() {
		return record.getSchema();
	}

	@Override
	public Object get(String attributeName) {
		return record.get(attributeName);
	}

	@Override
	public void put(String attributeName, Object value) {
		throw new UnsupportedOperationException("AvroSamRecord objects are immutable");
	}

	@Override
	public String getName() {
		return getStringAttributeOrThrow("qname");
	}

	@Override
	public String getReferenceName() {
		return getStringAttributeOrThrow("rname");
	}

	@Override
	public int getReferenceStartPosition() {
		return getIntAttributeOrThrow("pos");
	}

	@Override
	public int getReferenceEndPosition() {
		return annotation.getReferenceEndPosition();
	}

	@Override
	public Iterator<SingleInterval> getBlocks() {
		return annotation.getBlocks();
	}

	@Override
	public int getNumberOfBlocks() {
		return annotation.getNumberOfBlocks();
	}

	@Override
	public int size() {
		return annotation.size();
	}

	@Override
	public Strand getOrientation() {
		return annotation.getOrientation();
	}

	@Override
	public int getRelativePositionFrom5PrimeOfFeature(int referenceStart) {
		return annotation.getRelativePositionFrom5PrimeOfFeature(referenceStart);
	}

	@Override
	public AnnotationCollection<DerivedAnnotation<? extends Annotation>> getWindows(int windowSize, int stepSize) {
		return annotation.getWindows(windowSize, stepSize);
	}

	@Override
	public void setOrientation(Strand orientation) {
		throw new UnsupportedOperationException("AvroSamRecord objects are immutable");
	}

	@Override
	public Collection<? extends ReadFlag> getReadFlags() {
		throw new UnsupportedOperationException();
	}

	@Override
	public int getNumHits() {
		return getIntAttributeOrThrow("tagNH");
	}

	@Override
	public int getMappingQuality() {
		return getIntAttributeOrThrow("mapq");
	}
	
	@Override
	public String toString(){
		return toSAM();
	}
	
	private static String toTag(String fieldName) {
		return fieldName.replaceAll("tag", "");
	}
	
	private Optional<String> toStringTagWithValue(String fieldName) {
		Optional<String> val = getStringAttribute(fieldName);
		if(val.isPresent()) return Optional.of(toTag(fieldName) + ":Z:" + val.get());
		else return Optional.empty();
	}
	
	private Optional<String> toIntTagWithValue(String fieldName) {
		Optional<Integer> i = getIntAttribute(fieldName);
		if(i.isPresent()) return Optional.of(toTag(fieldName) + ":i:" + i.get().toString());
		else return Optional.empty();
	}
	
	/**
	 * Get the field value as a SAM tag string for SAM format
	 * @param field The field
	 * @return SAM tag in format TAG:TYPE:VALUE
	 */
	private Optional<String> asTag(Field field) {
		String fieldName = field.name();
		switch(field.schema().getType()) {
		case UNION: // Used for optional tags; union of some type and NULL
			switch(field.schema()
					.getTypes() // List<Schema>
					.stream() // Stream<Schema>
					.filter(schema -> schema.getType() != Schema.Type.NULL) // Stream<Schema>
					.collect(Collectors.toList()) // List<Schema>
					.get(0) // The first non-null field
					.getType()) { // The type of the first non-null field
			
						case STRING:
							return toStringTagWithValue(fieldName);
						case INT:
							return toIntTagWithValue(fieldName);
						default:
							throw new IllegalArgumentException("Field type not supported: " + field.schema().getType().getName());
			}
			
		case STRING:
			return toStringTagWithValue(fieldName);
		case INT:
			return toIntTagWithValue(fieldName);
		default:
			throw new IllegalArgumentException("Field type not supported: " + field.schema().getType().getName());
		}
	}
	
	/**
	 * @return Formatted SAM record
	 */
	public String toSAM() {
		return getStringAttributeOrThrow("qname") + "\t"
				+ getStringAttributeOrThrow("flag") + "\t"
				+ getStringAttributeOrThrow("rname") + "\t"
				+ getStringAttributeOrThrow("pos") + "\t"
				+ getStringAttributeOrThrow("mapq") + "\t"
				+ getStringAttributeOrThrow("cigar") + "\t"
				+ getStringAttributeOrThrow("rnext") + "\t"
				+ getStringAttributeOrThrow("pnext") + "\t"
				+ getStringAttributeOrThrow("tlen") + "\t"
				+ getStringAttributeOrThrow("seq") + "\t"
				+ getStringAttributeOrThrow("qual") + "\t"
				+ record
					.getSchema() // https://avro.apache.org/docs/1.7.6/api/java/org/apache/avro/Schema.html
					.getFields()
					.stream()
					.filter(field -> field.name().startsWith("tag")) // Stream of fields starting with "tag"
					.map(field -> asTag(field))
					.filter(opt -> opt.isPresent())
					.map(opt -> opt.get())
					.collect(Collectors.joining("\t"));
	}
	
	/**
	 * Get the records formatted as a SAM file, one record per line
	 * @param records Stream of records
	 * @return Formatted SAM records, one per line
	 */
	public static String toSAM(Stream<AvroSamRecord> records) {
		return records.map(record -> record.toSAM()).collect(Collectors.joining("\n"));
	}
	
	/**
	 * Write records to a SAM file
	 * @param records Stream of records
	 * @param header Header to use
	 * @param outputSam Output sam file to write to
	 */
	public static void writeToSAM(Stream<AvroSamRecord> records, SAMFileHeader header, File outputSam) {
		try {
			FileWriter writer = new FileWriter(outputSam);
			writer.write(header.getTextHeader());
			writer.write(toSAM(records));
			writer.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	@Override
	public boolean equals(Object other)
	{
		if(!(other instanceof Annotation)) {
			return false;
		}
		Annotation b = (Annotation)other;
		if(!getName().equals(b.getName())) return false;
		return compareTo(b) == 0;
	}
	
	@Override
	public int hashCode()
	{
		return new HashCodeBuilder(31,37).append(getName()).append(getReferenceName()).append(getReferenceStartPosition())
				.append(getReferenceEndPosition()).append(getOrientation()).append(getNumberOfBlocks()).toHashCode();
	}


	
}

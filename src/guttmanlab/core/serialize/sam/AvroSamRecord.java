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

import java.util.Collection;
import java.util.Iterator;

import net.sf.samtools.SAMRecord;

import org.apache.avro.Schema;
import org.apache.avro.generic.GenericRecord;
import org.apache.commons.lang.builder.HashCodeBuilder;

public class AvroSamRecord extends BlockedAnnotation implements GenericRecord, MappedFragment {
	
	private GenericRecord record;
	private Annotation annotation;
	private boolean firstReadTranscriptionStrand;
	public static int MIN_MAPPING_QUALITY = 0;
	public static int MAX_MAPPING_QUALITY = Integer.MAX_VALUE; // mapq=255 means mapping quality not available
	
	public AvroSamRecord(GenericRecord genericRecord) {
		this(genericRecord, true);
	}
	
	public AvroSamRecord(GenericRecord genericRecord, boolean firstReadIsTranscriptionStrand) {
		
		// Get basic attributes
		firstReadTranscriptionStrand = firstReadIsTranscriptionStrand;
		record = genericRecord;
		if(!mappingQualityIsOk()) {
			throw new IllegalStateException("Mapping quality not valid: " + getMappingQuality());
		}
		String cigar = getStringAttribute("cigar");
		String chr = getReferenceName();
		int start = getReferenceStartPosition();
		String name = getName();
		
		// Determine strand
		int flag = getIntAttribute("flag");
		SAMFlagDecoder decoder = new SAMFlagDecoder(flag);
		boolean isPaired = decoder.templateHasMultipleSegmentsInSequencing();
		boolean isFirst = decoder.firstSegmentInTemplate();
		boolean plusStrand = !decoder.seqIsReverseComplemented();
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
	 * @param mapq Mapping quality
	 * @return True iff the mapping quality is at least the minimum and less than 255 (unknown mapq)
	 */
	private static boolean mappingQualityIsOk(int mapq) {
		return mapq >= MIN_MAPPING_QUALITY && mapq <= MAX_MAPPING_QUALITY;
	}
	
	/**
	 * @return True iff the mapping quality is at least the minimum and less than 255 (unknown mapq)
	 */
	public boolean mappingQualityIsOk() {
		return mappingQualityIsOk(getMappingQuality());
	}
	
	/**
	 * @param record Sam record
	 * @return True iff the mapping quality is at least the minimum and less than 255 (unknown mapq)
	 */
	public static boolean mappingQualityIsOk(SAMRecord record) {
		return mappingQualityIsOk(record.getMappingQuality());
	}
	
	public Object getAttribute(String attributeName) {
		return record.get(attributeName);
	}
	
	public String getStringAttribute(String attributeName) {
		return record.get(attributeName).toString();
	}
	
	public int getIntAttribute(String attributeName) {
		return (int) record.get(attributeName);
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
		return getStringAttribute("qname");
	}

	@Override
	public String getReferenceName() {
		return getStringAttribute("rname");
	}

	@Override
	public int getReferenceStartPosition() {
		return getIntAttribute("pos");
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
		return getIntAttribute("tagNH");
	}

	@Override
	public int getMappingQuality() {
		return getIntAttribute("mapq");
	}
	
	@Override
	public String toString(){
		return toBED(0,0,0);
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

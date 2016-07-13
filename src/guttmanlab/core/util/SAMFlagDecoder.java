package guttmanlab.core.util;

/**
 * Utility class for decoding an integer SAM flag and getting the encoded information
 * @author prussell
 *
 */
public class SAMFlagDecoder {
	
	private String binaryString;
	private int stringLength;
	
	private int posMultipleSegments = 0;
	private int posEachSegmentProperlyAligned = 1;
	private int posSegmentUnmapped = 2;
	private int posNextSegmentInTemplateUnmapped = 3;
	private int posSeqReverseComplemented = 4;
	private int posSeqOfNextSegmentInTemplateReversed = 5;
	private int posFirstSegmentInTemplate = 6;
	private int posLastSegmentInTemplate = 7;
	private int posSecondaryAlignment = 8;
	private int posNotPassingQualityControls = 9;
	private int posPcrOrOpticalDuplicate = 10;
	private int posSupplementaryAlignment = 11;
	
	/**
	 * @param flag The int flag
	 */
	public SAMFlagDecoder(int flag) {
		if(flag < 0) {
			throw new IllegalArgumentException("Flag must be >= 0");
		}
		binaryString = Integer.toBinaryString(flag);
		stringLength = binaryString.length();
	}
	
	/**
	 * @return True if there are multiple segments (not single end)
	 */
	public boolean templateHasMultipleSegmentsInSequencing() {
		return bitIsTrue(posMultipleSegments);
	}
	
	/**
	 * @return True if each segment is properly aligned
	 */
	public boolean eachSegmentProperlyAligned() {
		return bitIsTrue(posEachSegmentProperlyAligned);
	}
	
	/**
	 * @return True if this segment is unmapped
	 */
	public boolean segmentUnmapped() {
		return bitIsTrue(posSegmentUnmapped);
	}
	
	/**
	 * @return True if next segment in template is unmapped
	 */
	public boolean nextSegmentInTemplateUnmapped() {
		return bitIsTrue(posNextSegmentInTemplateUnmapped);
	}
	
	/**
	 * @return True if this segment is mapped to minus strand
	 */
	public boolean seqIsReverseComplemented() {
		return bitIsTrue(posSeqReverseComplemented);
	}
	
	/**
	 * @return True if next segment in template is mapped to minus strand
	 */
	public boolean seqOfNextSegmentInTemplateIsReverseComplemented() {
		return bitIsTrue(posSeqOfNextSegmentInTemplateReversed);
	}
	
	/**
	 * @return True if this is the first segment in template
	 */
	public boolean firstSegmentInTemplate() {
		return bitIsTrue(posFirstSegmentInTemplate);
	}
	
	/**
	 * @return True if this is the last segment in template
	 */
	public boolean lastSegmentInTemplate() {
		return bitIsTrue(posLastSegmentInTemplate);
	}
	
	/**
	 * @return True if this is not the primary alignment
	 */
	public boolean secondaryAlignment() {
		return bitIsTrue(posSecondaryAlignment);
	}
	
	/**
	 * @return True if this mapping is marked as a PCR or optical duplicate
	 */
	public boolean pcrOrOpticalDuplicate() {
		return bitIsTrue(posPcrOrOpticalDuplicate);
	}
	
	/**
	 * @return True if this mapping does not pass quality controls
	 */
	public boolean notPassingQualityControls() {
		return bitIsTrue(posNotPassingQualityControls);
	}
	
	/**
	 * @return True if this is a supplementary alignment
	 */
	public boolean supplementaryAlignment() {
		return bitIsTrue(posSupplementaryAlignment);
	}
	
	private boolean bitIsTrue(int pos) {
		if(stringLength <= pos) {
			return false;
		}
		char c = binaryString.charAt(stringLength - 1 - pos);
		return c == '1';
	}
	
}
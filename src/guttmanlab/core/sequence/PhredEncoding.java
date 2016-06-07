package guttmanlab.core.sequence;

/**
 *  A representation of a Phred quality score encoding scheme.
 *  <p><p>
 *  PhredEncodings are mainly used to convert format-specific quality strings to
 *  a more general representation of a quality score. Currently, each encoding
 *  has an offset value as a member variable, and this offset is subtracted
 *  from each character in the quality string (more precisely, from each
 *  character's ASCII value) to get the actual quality score. As an example,
 *  Sanger's encoding maps the characters in the range from '!' to 'I' (or
 *  from ASCII 33 to ASCII 73) to a Phred score range from 0 to 40. The
 *  Sanger offset is therefore 33. Supported encodings are:
 *  <p><ul>
 *  <li>SANGER (Sanger format)
 *  <li>ILLUMINA_13 (Illumina 1.3 format)
 *  </ul><p>
 */
public enum PhredEncoding {

    SANGER(33),
    ILLUMINA_13(64);
    
    private final int offset;

    private PhredEncoding(int offset) {
        this.offset = offset;
    }
    
    /**
     * Returns the offset value of this PhredEncoding
     */
    public final int offset() {
        return this.offset;
    }
}
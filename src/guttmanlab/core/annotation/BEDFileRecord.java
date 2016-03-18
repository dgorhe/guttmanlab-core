package guttmanlab.core.annotation;

import java.awt.Color;
import java.util.Iterator;

import guttmanlab.core.annotationcollection.AnnotationCollection;

/**
 * The {@code BEDFileRecord} class represents a feature given in a line of a BED file. Lines
 * from a BED file have three required fields (chromosome, start position, end position) and
 * nine additional optional fields. A {@code BEDFileRecord} object has values for all twelve
 * fields (sensible defaults for the nine optional fields, if they aren't provided).
 */
public class BEDFileRecord extends Gene implements AnnotationFileRecord {

	private final double score;
	private final int thickStart;
	private final int thickEnd;
	private final Color color;

	private final static double DEFAULT_SCORE = 0;
	private final static Color DEFAULT_COLOR = Color.BLACK;
	private static final int MAX_FIELDS = 12;
	
	/**
	 * Constructs a BED record from a given {@code Annotation}. The score field will be set to zero, and
	 * the color field will be set to black. The thickStart and thickEnd fields will both default
	 * to the start coordinate of the annotation so that there is no thick region.
	 * @param annot is the annotation to base this BED record on
	 */
	public BEDFileRecord(Annotation annot) {
		super(annot);
		thickStart = annot.getReferenceStartPosition();
		thickEnd = thickStart;
		score = DEFAULT_SCORE;
		color = DEFAULT_COLOR;
	}
	
	/**
	 * Constructs a BED record from a given {@code Gene}. The score field will be set to zero, and
	 * the color field will be set to black. The thickStart and thickEnd fields will be set to the
	 * start and end coordinate of the coding region. If there is no coding region, they will
	 * instead be set to the start coordinate of the annotation so that there is no thick region.
	 * @param gene
	 */
	public BEDFileRecord(Gene gene) {
		super(gene);
		score = DEFAULT_SCORE;
		color = DEFAULT_COLOR;
		Annotation cds = gene.getCodingRegion();
		if (cds == null) {
			cdsStartPos = -1;
			cdsEndPos = -1;
			thickStart = gene.getReferenceStartPosition();
			thickEnd = gene.getReferenceStartPosition();
		} else {
			cdsStartPos = cds.getReferenceStartPosition();
			cdsEndPos = cds.getReferenceEndPosition();
			thickStart = cdsStartPos;
			thickEnd = cdsEndPos;
		}
	}
	
	private BEDFileRecord(BEDBuilder builder) {
		super(builder.annot);
		score = builder.score;
		thickStart = builder.thickStart;
		thickEnd = builder.thickEnd;
		cdsStartPos = thickStart == thickEnd ? -1 : thickStart;
		cdsEndPos = thickStart == thickEnd ? -1 : thickEnd;
		color = builder.color;
	}
	
	/**
	 * @return the score of this BED record
	 */
	public double score() {
		return score;
	}
	
	public int thickStart() {
		return thickStart;
	}
	
	public int thickEnd() {
		return thickEnd;
	}
	
	/**
	 * Returns the color of this BED record
	 * @return the color of this BED record
	 */
	public Color color() {
		return color;
	}
	
	/**
	 * Converts this BEDFileRecord to a formatted (i.e., tab-delimited) string suitable for outputting
	 * to a BED file. This string does not have a trailing newline.
	 * @param numFields is the number of fields to limit the output string to. Must be either 3, 4, 5, 6, 8, 9, 12.
	 * @return the String representation of this BEDFileRecord, as would be found in a BED file
	 */
	public String toFormattedString(int numFields) {
		if (numFields < 3 || numFields > MAX_FIELDS || numFields == 7 || numFields == 10 || numFields == 11) {
			throw new IllegalArgumentException("Attempted to convert BED record to string, but requested numFields " + numFields
					+ ". Number of fields must be either 3, 4, 5, 6, 8, 9 or 12.");
		}

		// Required fields
		StringBuilder sb = new StringBuilder();
		sb.append(getReferenceName() + "\t" + getReferenceStartPosition() + "\t" + getReferenceEndPosition());
		int currentFieldNum = 3;
		if (currentFieldNum >= numFields) {
			return sb.toString();
		}
		
		// Name
		String name = getName();
		if (name.isEmpty()) {
			name = ".";
		}
		sb.append("\t" + name);
		currentFieldNum = 4;
		if (currentFieldNum >= numFields) {
			return sb.toString();
		}
		
		// Score
		sb.append("\t" + score);
		currentFieldNum = 5;
		if (currentFieldNum >= numFields) {
			return sb.toString();
		}
		
		Strand strand = getOrientation();
		if (strand != Strand.POSITIVE && strand != Strand.NEGATIVE) {
			sb.append("\t" + ".");
		} else {
			sb.append("\t" + strand.toString());
		}
		currentFieldNum = 6;
		if (currentFieldNum >= numFields) {
			return sb.toString();
		}
		
		sb.append("\t" + thickStart + "\t" + thickEnd);
		currentFieldNum = 8;
		if (currentFieldNum >= numFields) {
			return sb.toString();
		}
		
		sb.append("\t" + color.getRed() + "," + color.getGreen() + "," + color.getBlue());
		currentFieldNum = 9;
		if (currentFieldNum >= numFields) {
			return sb.toString();
		}
		  
		sb.append("\t" + getNumberOfBlocks() + "\t");
		Iterator<SingleInterval> blocks = getBlocks();
		while (blocks.hasNext()) {
			Annotation block = blocks.next();
			sb.append(block.size() + ","); // trailing comma after last block is OK
		}
		sb.append("\t");
		blocks = getBlocks();
		while (blocks.hasNext()) {
			Annotation block = blocks.next();
			sb.append((block.getReferenceStartPosition() - getReferenceStartPosition()) + ","); // trailing comma after last is OK
		}
		return sb.toString();
	}		
	
	/**
	 * Converts this BEDFileRecord to a formatted (i.e., tab-delimited) string suitable for outputting
	 * to a BED file. This string does not have a trailing newline. All twelve fields will be present
	 * in this string. If you wish to limit the number of fields output, see {@link #toFormattedString(int) toFormattedString(int)}
	 * @return the String representation of this BEDFileRecord, as would be found in a BED file
	 */
	@Override
	public String toFormattedString() {
		return toFormattedString(MAX_FIELDS);
	}
	
	/**
	 * Parses a String into a BEDFileRecord. The String should be whitespace-delimited and follow the same
	 * format as a single line from a BED file. Any trailing whitespace, including newline, is removed
	 * during parsing.
	 * @param s is the String to parse to a BEDFileRecord
	 * @return the BEDFileRecord which corresponds to the input String
	 * @throws IllegalArgumentException if the input has an invalid number of fields after splitting on tab
	 * @throws IllegalArgumentException if there is disagreement between the blockCount field and the numbers
	 * of blocks implied by the blockStarts and blockSizes fields
	 * @throws IllegalArgumentException if the orientation is present but is neither positive nor negative
	 * @throws IllegalArgumentException if the color is present but is not a valid RGB color triple
	 * @throws IllegalArgumentException if thickEnd occurs before thickStart, or if either fall outside of the
	 * interval defined by chromStart and chromEnd
	 */
	public static BEDFileRecord fromFormattedString(String s) {
		BEDBuilder bb;
		String[] fields = s.trim().split("\\s+");
		int numFields = fields.length;
		
		if (numFields < 3 || numFields == 7 || numFields == 10 || numFields == 11 || numFields > 12) {
			throw new IllegalArgumentException("fromFormattedString() was passed a String with "
					+ numFields + " fields. A properly formatted BED String must have between three and"
					+ " twelve fields, and cannot have seven, ten or eleven fields.");
		}
		
		// All fields present. This is an Annotation with multiple blocks
		if (numFields == 12) {
			String chrom = fields[0];
			int chromStart = Integer.parseInt(fields[1]);
			String name = fields[3];
			int blockCount = Integer.parseInt(fields[9]);
			int[] blockSizes = parseCommaSeparatedString(fields[10]);
			int[] blockStarts = parseCommaSeparatedString(fields[11]);
			
			if (blockStarts.length != blockCount || blockSizes.length != blockCount) {
				throw new IllegalArgumentException("Malformed BED String. blockCount = " + blockCount
						+ ", blockSizes = " + blockSizes.length + ", and blockStarts = " + blockStarts.length +
						". All should be equal.");
			}
			BlockedAnnotation annot = new BlockedAnnotation(name);
			for (int i = 0; i < blockCount; i++) {
				annot.addBlocks(new SingleInterval(chrom, chromStart + blockStarts[i], chromStart + blockStarts[i] + blockSizes[i]));
			}
			bb = new BEDBuilder(annot);
		
		// This is an Annotation with one block and a name.
		} else if (numFields >= 4) {
			bb = new BEDBuilder(new SingleInterval(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), Strand.UNKNOWN, fields[3]));
		
		// This is the simplest BED Annotation: one block with no name.
		} else {
			bb = new BEDBuilder(new SingleInterval(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2])));
			return bb.build();
		}
		
		// Annotation has been constructed. Add the rest of the fields.
		
		// Nothing else to add.
		if (numFields <= 4) {
			return bb.build();
		}
		
		// Add score
		if (numFields >= 5) {
			bb.score(Double.parseDouble(fields[4]));
			if (numFields == 5) {
				return bb.build();
			}
		}
		
		// Add orientation
		if (numFields >= 6) {
			if (!fields[5].equals("+") && !fields[5].equals("-") && !fields[5].equals(".")) {
				throw new IllegalArgumentException("The strand field of the formatted BED string must either"
						+ " be a '+', a '-', or a '.'");
			}
			Strand str = Strand.fromString(fields[5]);
			bb.strand(str);
			if (numFields == 6) {
				return bb.build();
			}
		}
		
		// Add line thickness.
		if (numFields >= 8) {
			bb.thickStart(Integer.parseInt(fields[6]));
			bb.thickEnd(Integer.parseInt(fields[7]));
			if (numFields == 8) {
				return bb.build();
			}
		}
		
		// Add color.
		if (numFields >= 9) {
			if (fields[8].equals(".")) {
				bb.color(DEFAULT_COLOR);
			} else {
				int[] colorVals = parseCommaSeparatedString(fields[8]);
				if (colorVals.length != 3) {
					throw new IllegalArgumentException("Formatted BED string has invalid color value: " + fields[8]);
				}
				bb.color(new Color(colorVals[0], colorVals[1], colorVals[2]));
			}
		}
		
		return bb.build();
	}
	
	// Helper method for fromFormattedString(). Handles parsing of the blockStarts, blockSizes and color BED fields.
	private static int[] parseCommaSeparatedString(String s) {
		String[] tmp = s.endsWith(",") ? s.substring(0, s.length() - 1).split(",") : s.split(",");
		int[] rtrn = new int[tmp.length];
		for (int i = 0; i < rtrn.length; i++) {
			rtrn[i] = Integer.parseInt(tmp[i]);
		}
		return rtrn;
	}
	
	/**
	 *  Returns a String representation of this BEDFileRecord by calling toString() on its component
	 *  Annotation. If you need a formatted line of text for a BED file, see {@link #toFormattedString() toFormattedString()};
	 */
	@Override
	public String toString() {
		// Method only overridden to add Javadoc. Any better ideas?
		return toString();
	}
	
	/**
	 * A builder class for creating {@code BEDFileRecord} objects. Below is an example of how to use this class.
	 * <p>
	 * {@code BEDFileRecord b = new BEDBuilder("chr1", 100, 300).score(20).strand(Strand.POSITIVE).color(Color.RED).build();}
	 */
	public static class BEDBuilder {
		// Required fields
		private final Annotation annot;
		
		// Optional fields
		private double score;
		private int thickStart;
		private int thickEnd;
		private Color color;
		
		/**
		 * Constructs a BEDBuilder based on an {@code Annotation} (the copied characteristics being the reference
		 * coordinates, the blocks, the orientation, and the name). The thickStart and thickEnd fields are set to
		 * the start coordinate of the annotation. The score field is set to zero, and the color field is set to black.
		 * This constructor makes a copy of the annotation passed to it.
		 * @param annot is the annotation to copy
		 */
		public BEDBuilder(final Annotation annot) {
			this.annot = new BlockedAnnotation(annot);
			thickStart = this.annot.getReferenceStartPosition();
			thickEnd = this.annot.getReferenceStartPosition();
			score = DEFAULT_SCORE;
			color = DEFAULT_COLOR;
		}

		/**
		 * Constructs a BEDBuilder based on a {@code Gene}. The thickStart and thickEnd fields are set to
		 * the corresponding coordinates of the gene's coding region, if it is specified. Otherwise, they
		 * are set to the start coordinate of the gene. The score field is set to zero, and the color field 
		 * is set to black. This constructor makes a copy of the gene passed to it.
		 * @param gene is the annotation to copy
		 */
		public BEDBuilder(final Gene gene) {
			annot = new BlockedAnnotation(gene);
			Annotation cds = gene.getCodingRegion();
			thickStart = cds == null ? gene.getReferenceStartPosition() : cds.getReferenceStartPosition();
			thickEnd = cds == null ? gene.getReferenceStartPosition() : cds.getReferenceEndPosition();
			score = DEFAULT_SCORE;
			color = DEFAULT_COLOR;
		}

		/**
		 * Constructs a BEDBuilder based on a contiguous annotation defined by the parameters. The thickStart
		 * and thickEnd fields are set to the start coordinate of the annotation. The score field is set to
		 * zero, and the color field is set to black. This constructor is identical to
		 * {@code BEDBuilder(new SingleInterval(chrom, chromStart, chromEnd))}.
		 * @param chrom is the reference name of the annotation
		 * @param chromStart is the start coordinate of the annotation
		 * @param chromEnd is the end coordinate of the annotation
		 */		
		public BEDBuilder(final String chrom, final int chromStart, final int chromEnd) {
			annot = new SingleInterval(chrom, chromStart, chromEnd);
			thickStart = annot.getReferenceStartPosition();
			thickEnd = annot.getReferenceStartPosition();
			score = DEFAULT_SCORE;
			color = DEFAULT_COLOR;		}
		
		public BEDBuilder score(double score) {
			this.score = score;
			return this;
		}
		
		public BEDBuilder strand(Strand str) {
			this.annot.setOrientation(str);
			return this;
		}
		
		public BEDBuilder thickStart(int thickStart) {
			this.thickStart = thickStart;
			return this;
		}
		
		public BEDBuilder thickEnd(int thickEnd) {
			this.thickEnd = thickEnd;
			return this;
		}
		
		public BEDBuilder color(Color color) {
			this.color = color;
			return this;
		}
		
		public BEDBuilder color(int r, int g, int b) {
			this.color = new Color(r, g, b);
			return this;
		}
		
		public BEDBuilder color(int rgb) {
			this.color = new Color(rgb);
			return this;
		}
		
		/**
		 * @return a new {@code BEDFileRecord} object with values matching those of this builder
		 */
		public BEDFileRecord build() {
			if (thickStart < annot.getReferenceStartPosition()) {
				throw new IllegalArgumentException("Attempted to create BED record with thickStart " + thickStart + " less than chromStart " + annot.getReferenceStartPosition() + ".");
			}
			if (thickStart > annot.getReferenceEndPosition()) {
				throw new IllegalArgumentException("Attempted to create BED record with thickStart " + thickStart + " greater than chromEnd " + annot.getReferenceEndPosition() + ".");
			}
			if (thickEnd < annot.getReferenceStartPosition()) {
				throw new IllegalArgumentException("Attempted to create BED record with thickEnd " + thickEnd + " less than chromStart " + annot.getReferenceStartPosition() + ".");
			}
			if (thickEnd > annot.getReferenceEndPosition()) {
				throw new IllegalArgumentException("Attempted to create BED record with thickEnd " + thickEnd + " greater than chromEnd " + annot.getReferenceEndPosition() + ".");
			}
			if (thickStart > thickEnd) {
				throw new IllegalArgumentException("Attempted to create BED record with thickStart " + thickStart + " greater than thickEnd " + thickEnd + ".");
			}
			return new BEDFileRecord(this);
		}
	}
	
	@Override
	public int getRelativePositionFrom5PrimeOfFeature(int referenceStart) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public AnnotationCollection<DerivedAnnotation<? extends Annotation>> getWindows(int windowSize, int stepSize) {
		// TODO Auto-generated method stub
		return null;
	}
}
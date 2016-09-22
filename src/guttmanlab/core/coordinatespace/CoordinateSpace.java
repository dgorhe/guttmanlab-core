package guttmanlab.core.coordinatespace;

import guttmanlab.core.annotation.Annotation;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;

/**
 * A coordinate space that is a subset of a background genome
 * @author prussell
 *
 */
public class CoordinateSpace {
	
	public static CoordinateSpace MM9 = new CoordinateSpace(GenomeSize.MM9);
	public static CoordinateSpace MM10 = new CoordinateSpace(GenomeSize.MM10);
	public static CoordinateSpace HG19 = new CoordinateSpace(GenomeSize.HG19);
	
	private enum Precomputed {
				
		MM9("mm9"), MM10("mm10"), HG19("hg19");
		
		private String name;
		
		Precomputed(String name) {this.name = name;}
		
		public static String commaSeparatedString = Arrays.asList(values()).stream().map(c -> c.toString()).collect(Collectors.joining(", "));
		
		public String toString() {return name;}
		
		public static Precomputed fromString(String genomeName) {
			for(Precomputed p : values()) {
				if(p.toString().equals(genomeName)) return p;
			}
			throw new IllegalArgumentException("Supported genome names: " + commaSeparatedString);
		}
		
	}
	
	/**
	 * Get a coordinate space by passing an assembly name
	 * @param assembly Assembly e.g. "mm10"
	 * @return The coordinate space for the assembly
	 */
	public static CoordinateSpace forGenome(String assembly) {
		if(assembly.equals(Precomputed.MM9.toString())) return MM9;
		else if(assembly.equals(Precomputed.MM10.toString())) return MM10;
		else if(assembly.equals(Precomputed.HG19.toString())) return HG19;
		else throw new IllegalArgumentException("Assembly " + assembly + " not supported. Options: " + Precomputed.commaSeparatedString);
	}
	
	/**
	 * Description of reference sequences in this coordinate space
	 * Key is reference name, value is reference length
	 */
	private Map<String, Integer> refSizes;
		
	/**
	 * 
	 * @param referenceSizesFile File containing reference names and lengths
	 */
	public CoordinateSpace(String referenceSizesFile){
		this.refSizes=getRefSeqLengthsFromTable(referenceSizesFile);
	}

	/**
	 * @param sizes Map of reference name to reference size
	 */
	public CoordinateSpace(Map<String, Integer> sizes){
		this.refSizes=sizes;
	}
	
	/**
	 * Create from the reference dictionary in a SAM header
	 * @param fileHeader SAM header
	 */
	public CoordinateSpace(SAMFileHeader fileHeader) {
		this.refSizes=getRefSeqLengthsFromSamHeader(fileHeader);
	}
	
	/**
	 * Tests whether a region is present within this CoordinateSpace
	 * @param annotation The annotation to test
	 * @return boolean true=contained in space, false=not contained
	 */
	public boolean contains(Annotation annotation){
		//Is contained if annotation.getReferenceName() is in references and positions are within size
		if(this.refSizes.containsKey(annotation.getReferenceName())){
			int size=refSizes.get(annotation.getReferenceName());
			if(annotation.getReferenceEndPosition()<size){return true;}
		}
		return false;
	}

	@Override
	public boolean equals(Object o){
		// FIXME Auto-generated method stub
		throw new UnsupportedOperationException("TODO");
	}
	
	/**
	 * @return Map associating each reference name with sequence length
	 */
	public Map<String, Integer> getRefSeqLengths() {
		return refSizes;
	}
	
	/**
	 * @return Total length of all reference sequences
	 */
	public long getTotalReferenceLength() {
		long rtrn = 0;
		for(String chr : refSizes.keySet()) {
			rtrn += refSizes.get(chr).intValue();
		}
		return rtrn;
	}
	
	/**
	 * Get the lengths of the reference sequences from a SAM file header
	 * @param header SAM file header
	 * @return Map associating each reference name with sequence length
	 */
	private Map<String, Integer> getRefSeqLengthsFromSamHeader(SAMFileHeader header) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
		if (records.size() > 0) {
			for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
				String chr = rec.getSequenceName();
				int size=rec.getSequenceLength();
				rtrn.put(chr, size);
		    	}
		    }
		return rtrn;
	}
	
	/**
	 * Parse the reference sizes file
	 * @param referenceSizesFile Tab-delimited file with reference names (ie chromosomes) and lengths
	 * @return Map associating each reference name with sequence length
	 */
	private Map<String, Integer> getRefSeqLengthsFromTable(String referenceSizesFile) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		try{	
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(referenceSizesFile)));
			String nextLine;
			while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
				String[] tokens=nextLine.split("\t| +");
				rtrn.put(tokens[0], new Integer(tokens[1]));
			}
			reader.close();
		}catch(IOException ex){ex.printStackTrace();}
		return rtrn;
	}

	/**
	 * Convert the CoordinateSpace into a SAMFileHeader
	 * @return The SAMFileHeader
	 */
	public SAMFileHeader getBAMFileHeader() {
		SAMFileHeader header=new SAMFileHeader();

		//add sequences
		for(String refSeq: this.refSizes.keySet()){
			int size=this.refSizes.get(refSeq)+1; //TODO I think this is correct for the header
			SAMSequenceRecord seq=new SAMSequenceRecord(refSeq, size);
			header.addSequence(seq);
		}

		//Set sort order
		header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
				
		return header;
	}
	
	
	
}

package guttmanlab.core.simulation;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;

/**
 * This class represents the coordinate-ranges of a genome assembly.
 * <p>
 * Currently, the primary use of the class is to construct headers for SAM and
 * BAM files.
 */
public final class CoordinateSpace {

    public final static CoordinateSpace MM9 = new CoordinateSpace(AssemblySize.MM9);
    public final static CoordinateSpace MM10 = new CoordinateSpace(AssemblySize.MM10);
    public final static CoordinateSpace HG19 = new CoordinateSpace(AssemblySize.HG19);
    public final static CoordinateSpace HG38 = new CoordinateSpace(AssemblySize.HG38);

    private final Map<String, Integer> refSizes;

    /**
     * Class constructor.
     * <p>
     * Constructs an instance of a <code>CoordinateSpace</code> with the given
     * chromosome-to-size mapping.
     * 
     * @param sizes - a mapping from chromosome-name to size-in-bp
     * @throws NullPointerException if the mapping is null
     */
    public CoordinateSpace(Map<String, Integer> sizes) {
        Objects.requireNonNull(sizes, "Attempted to construct a coordinate"
                + " space with a null size mapping.");
        this.refSizes = sizes;
    }
    
   
    
    
    
   
    
    
    public static double getGenomeLength(SAMFileHeader header) {
    	double length=0;
    	List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
        if (records.size() > 0) {
            for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                String chr = rec.getSequenceName();
                int size = rec.getSequenceLength();
                length+=size;
            }
        }
        return length;
    }
    
    
    public static Map<String, Integer> getGenomeLengths(SAMFileHeader header) {
    	double length=0;
    	Map<String, Integer> rtrn=new TreeMap<String, Integer>();
    	List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
        if (records.size() > 0) {
            for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                String chr = rec.getSequenceName();
                int size = rec.getSequenceLength();
                //length+=size;
                rtrn.put(chr,  size);
            }
        }
        return rtrn;
    }
    
    /**
     * Constructs and returns a new htsjdk <code>SAMFileHeader</code> object
     * that corresponds with this coordinate space.
     * 
     * @returns the <code>SAMFileHeader</code> represented by this coordinate
     * space.
     */
    public SAMFileHeader getSAMFileHeader() {
        SAMFileHeader header = new SAMFileHeader();

        for (Map.Entry<String, Integer> entry : refSizes.entrySet()) {
            int size = entry.getValue();
            SAMSequenceRecord seq = new SAMSequenceRecord(entry.getKey(), size);
            header.addSequence(seq);
        }

        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
                
        return header;
    }
    
    /**
     * Constructs and returns a new net.sf.samtools <code>SAMFileHeader</code> object
     * that corresponds with this coordinate space.
     * 
     * @returns the <code>SAMFileHeader</code> represented by this coordinate
     * space.
     */
    public net.sf.samtools.SAMFileHeader getSAMFileHeaderOldVersion() {
        net.sf.samtools.SAMFileHeader header = new net.sf.samtools.SAMFileHeader();

        for (Map.Entry<String, Integer> entry : refSizes.entrySet()) {
            int size = entry.getValue();
            net.sf.samtools.SAMSequenceRecord seq = new net.sf.samtools.SAMSequenceRecord(entry.getKey(), size);
            header.addSequence(seq);
        }

        header.setSortOrder(net.sf.samtools.SAMFileHeader.SortOrder.coordinate);
                
        return header;
    }
    
    /**
     * Returns a string representation of this coordinate space.
     * <p>
     * The exact details of this representation are unspecified and subject to
     * change, but it will typically include a listing of the contained
     * references and the size of each.
     * 
     * @return a string representation of this coordinate space
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Entry<String, Integer> pair : refSizes.entrySet()) {
            sb.append(pair.getKey() + ": " + pair.getValue() +
                    System.lineSeparator());
        }
        return sb.toString();
    }
    
    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }
        
        if (!(o instanceof CoordinateSpace)) {
            return false;
        }
        
        CoordinateSpace other = (CoordinateSpace) o;
        return other.refSizes.equals(this.refSizes);
    }
    
    @Override
    public int hashCode() {
        return refSizes.hashCode();
    }

	public Map<String, Integer> getRefSizes() {
		return this.refSizes;
	}

	
	public List<String> getBins(int binResolution, String chr) {
		List<String> rtrn=new ArrayList<String>();
		
		
			int size=this.getRefSizes().get(chr);
			int length=size/binResolution;
			for(int i=0; i<length; i++){
				int start=i*binResolution;
				int end=start+binResolution;
				SingleInterval region=new SingleInterval(chr, start, end);
				//System.err.println(region.toUCSC());
				rtrn.add(region.toUCSC());
			}
		
		return rtrn;
	}
	
	public List<String> getBins(int binResolution, SingleInterval includeRegion) {
		List<String> rtrn=new ArrayList<String>();
		
			String chr=includeRegion.getReferenceName();
			int size=this.getRefSizes().get(chr);
			int length=size/binResolution;
			for(int i=0; i<length; i++){
				int start=i*binResolution;
				int end=start+binResolution;
				SingleInterval region=new SingleInterval(chr, start, end);
				//System.err.println(region.toUCSC());
				if(region.overlaps(includeRegion)){
					rtrn.add(region.toUCSC());
				}
			}
		
		return rtrn;
	}
	
	public List<String> getBins(int binResolution) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String chr: this.getRefSizes().keySet()){
			int size=this.getRefSizes().get(chr);
			int length=size/binResolution;
			for(int i=0; i<=length; i++){
				int start=i*binResolution;
				int end=start+binResolution;
				SingleInterval region=new SingleInterval(chr, start, end);
				//System.err.println(region.toUCSC());
				rtrn.add(region.toUCSC());
			}
		}
		return rtrn;
	}
	
	public Map<String, IntervalTree<SingleInterval>> getBinTree(int binResolution) {
		Map<String, IntervalTree<SingleInterval>> rtrn=new TreeMap<String, IntervalTree<SingleInterval>>();
		
		for(String chr: this.getRefSizes().keySet()){
			IntervalTree<SingleInterval> tree=new IntervalTree<SingleInterval>();
			int size=this.getRefSizes().get(chr);
			//int length=size/binResolution;
			for(int i=0; i<size; i++){
				int start=i;
				int end=start+binResolution;
				SingleInterval region=new SingleInterval(chr, start, end);
				tree.put(start, end, region);
			}
			rtrn.put(chr, tree);
		}
		return rtrn;
	}

	
	
}
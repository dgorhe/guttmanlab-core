package guttmanlab.core.sharp;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class FilterBam {

	
	/***
	 * Remove blacklisted regions from bam file
	 * @param bamFile
	 * @param geneTree blacklist
	 * @param save
	 */
	public FilterBam(File bamFile, Map<String, IntervalTree<Gene>> geneTree, String save) {
		//Filter reads overlapping genes
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save));
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			
			if(!overlapsGene(record, geneTree)) {
				writer.addAlignment(record);
			}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		writer.close();
		
	}
	
	
	public FilterBam(File bamFile, String save) {
		//Filter reads overlapping genes
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save));
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			
			if(record.getFirstOfPairFlag()) {
				writer.addAlignment(record);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		writer.close();
		
	}
	
	
	private static void excludeReads(File bamFile, Map<String, IntervalTree<Gene>> geneTree, String save) {
		//Filter reads overlapping genes
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save));
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			
			if(!overlapsSingleExon(record, geneTree)) {
				writer.addAlignment(record);
			}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		writer.close();
		
	}
	
	public FilterBam(File bamFile, String chr, String save) {
		//Filter reads overlapping genes
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save));
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(record.getReferenceName().equals(chr)) {
				writer.addAlignment(record);
			}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		writer.close();
		
	}
	
	private static boolean overlapsSingleExon(SAMRecord record, Map<String, IntervalTree<Gene>> geneTree) {
		SAMFragment frag=new SAMFragment(record);
		if(frag.isSpliced()) {return false;}
		IntervalTree<Gene> tree=geneTree.get(record.getReferenceName());
		if(tree!=null){
			Iterator<Gene> genes=tree.overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
			while(genes.hasNext()){
				Gene gene=genes.next();
				Strand o=frag.getOrientation();
				if(o.equals(gene.getOrientation())){
					if(gene.overlaps(frag) || frag.overlaps(gene)) {
					//if(gene.getReferenceStartPosition()>=record.getAlignmentStart() && gene.getReferenceEndPosition()<=record.getAlignmentEnd()) {
						return true;
					}
				}
			}
		}
			return false;	
	}

	private boolean overlapsGene(SAMRecord record, Map<String, IntervalTree<Gene>> geneTree) {
		IntervalTree<Gene> tree=geneTree.get(record.getReferenceName());
		if(tree!=null){
			Iterator<Gene> genes=tree.overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
			while(genes.hasNext()){
				Gene gene=genes.next();
				Strand o=SAMFragment.getOrientation(record);
				SAMFragment fragment=new SAMFragment(record);
				if(o.equals(gene.getOrientation())){
					if(fragment.overlaps(gene)) {
					//if(record.getAlignmentStart()>=gene.getReferenceStartPosition() && record.getAlignmentEnd()<=gene.getReferenceEndPosition()) {
						return true;
					}
				}
			}
		}
			return false;	
	}
	
	private static List<File> getBAMs(String string) {
		List<File> rtrn=new ArrayList<File>();
		
		File[] files=new File(string).listFiles();
		for(int i=0; i<files.length; i++) {
			if(files[i].getName().endsWith(".bam")) {rtrn.add(files[i]);}
	
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			File bam=new File(args[0]);
			Map<String, IntervalTree<Gene>> tree=BEDFileIO.loadTree(args[1]);
			String save=args[2];
			new FilterBam(bam, tree, save);
		}
		else {System.err.println(usage);}
		
		
		
	}
	static String usage=" args[0]=bam file \n args[1]=BED to exclude \n args[2]=save";
}

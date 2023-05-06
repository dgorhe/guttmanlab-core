package guttmanlab.core.clap.old;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.simulation.CoordinateSpace;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.CloseableIterator;

public class ExcludeReadsInList {

	public ExcludeReadsInList(File bamFile, Map<String, Collection<Gene>> excludeRegions, String save){
		BAMSingleReadCollection bam1=new BAMSingleReadCollection(bamFile); //TODO Get fragments
		Map<String, IntervalTree<Gene>> excludeTree=makeGeneTree(excludeRegions);
		SAMFileHeader header=bam1.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter alignmentWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save));
		
		
		CloseableIterator<SAMFragment> reads=bam1.sortedIterator();
		
		
		
		int counter=0;
		while(reads.hasNext()){
			SAMFragment read=reads.next();
			if(excludeTree.containsKey(read.getReferenceName())){
				boolean overlaps=excludeTree.get(read.getReferenceName()).hasOverlappers(read.getReferenceStartPosition(), read.getReferenceEndPosition());
				if(!overlaps){
					alignmentWriter.addAlignment(read.getSamRecord());
				}
			}
			counter++;
			if(counter%1000000 ==0){System.err.println("V2 "+counter);}
		}
		
		reads.close();
		alignmentWriter.close();
	}
	
	private Map<String, IntervalTree<Gene>> makeGeneTree(Map<String, Collection<Gene>> regions) {
		Map<String, IntervalTree<Gene>> rtrn=new TreeMap<String, IntervalTree<Gene>>();
		
		for(String chr: regions.keySet()){
			IntervalTree<Gene> tree=new IntervalTree<Gene>();
			for(Gene gene: regions.get(chr)){
				tree.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
			}
			rtrn.put(chr, tree);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File bam=new File(args[0]);
			Map<String, Collection<Gene>> repeats=BEDFileIO.loadRegionsFromFileByChr(args[1]);
			String save=args[2];
			new ExcludeReadsInList(bam, repeats, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam \n args[1]=excluded regions (BED) \n args[2]=save (.bam)";
}

package guttmanlab.core.orfcloning;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.simulation.CoordinateSpace;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ReadsOverlapping {

	public ReadsOverlapping(File bamFile, Map<String, IntervalTree<Gene>> tree, String save){
		SamReader reader=SamReaderFactory.makeDefault().open((bamFile));
		SAMFileHeader originalHeader=reader.getFileHeader();
		
		SAMFileWriter BAMWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(originalHeader, false, new File(save));
		
		
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			
			if(overlaps(tree, read)){
				BAMWriter.addAlignment(read);
			}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		BAMWriter.close();
		reads.close();
		System.err.println("Done");
		
		
	}

	private boolean overlaps(Map<String, IntervalTree<Gene>> tree, SAMRecord read) {
		String chr=read.getReferenceName();
		if(tree.containsKey(chr)){
			return tree.get(read.getReferenceName()).hasOverlappers(read.getAlignmentStart(), read.getAlignmentEnd());
		}
		return false;
	}
	
	public static void main(String[] args) throws IOException{
		File bam=new File(args[0]);
		Map<String, IntervalTree<Gene>> tree=BEDFileIO.loadTree(args[1]);
		String save=args[2];
		new ReadsOverlapping(bam, tree, save);
	}
	
}

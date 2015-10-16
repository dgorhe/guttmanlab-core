package guttmanlab.core.examples;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.predicate.SecondReadFilter;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;

public class findTSSbetter {

    public static void main(String[] args) throws IOException, Exception {
    	
    	//use read start positions instead of read count
    	boolean useReadStarts = true;
    	
    	// IO setup
    	BAMPairedFragmentCollection fivePrimeReads = new BAMPairedFragmentCollection(new File("/storage/Users/cburghard/Projects/071515_TSS/cage.newhead.cleaned.bam"));
    	fivePrimeReads.addFilter(new SecondReadFilter<PairedMappedFragment<SAMFragment>>());
    	FileWriter writer = new FileWriter("/storage/Users/cburghard/Projects/071515_TSS/gencode_startSites.bed");
        String featureFile = ( "/storage/Users/cburghard/Projects/071515_TSS/gencode_mm9_liftover_head.bed" ); ///storage/Annotations/RefSeq_lincV3/MouseLincRNA.v3.ChromatinDefinedWithNames.bed";
		BEDFileIO io =  new BEDFileIO( "/storage/shared/CoreTestData/refspace.txt" ); 
		AnnotationCollection<Gene> features = io.loadFromFile( featureFile );
		CloseableIterator<Gene> fiter = features.sortedIterator();
		
		
		if(useReadStarts){
			while ( fiter.hasNext() )
			{
				HashMap<Integer,Integer> startSites = new HashMap<Integer,Integer>();
				Gene gene = fiter.next();
				int max = 0;
				//Extend the gene 10kb before 5' end
				BlockedAnnotation region = extendGene(gene);
				
				CloseableIterator<PairedMappedFragment<SAMFragment>> iter = fivePrimeReads.sortedIterator( region,true );
				//for each extended gene
				while ( iter.hasNext() )
				{
					PairedMappedFragment<SAMFragment> s = iter.next();
					
					int count = 1; 
					int start = 0;
					
					if ( region.getOrientation().equals(Strand.NEGATIVE) )
						start = s.getReferenceEndPosition();
					else start = s.getReferenceStartPosition();
					
					if(startSites.containsKey(start))
						count = startSites.get(start)+1;
					startSites.put(start,count);
					if (count > max)
						max = count;
				}
				
				
				//set cutoff, print start sites
				int regionSum = 0;
				for(int val : startSites.values())
				{
					regionSum += val;
				}
				
				for(int key : startSites.keySet())
				{
					//if ( startSites.get(key)/regionSum >= .1 )
						writer.write(region.getReferenceName()+"\t"+region.getName()+"\t"+(key)+"\t"+(key+1)+"\t"+startSites.get(key)+"\n"); //output for paraclu
				}
				iter.close();
			}
		}
		else {

			
		}
		fiter.close();
    	writer.close();
    	System.out.println("done");
    }

    private static BlockedAnnotation extendGene(Gene gene)
    {
		BlockedAnnotation region = new BlockedAnnotation(gene.getName());
    	int numBlocks = gene.getNumberOfBlocks();
		Iterator<SingleInterval> blocks = gene.getBlocks();
		int countBlocks=0;
		Strand orientation = gene.getOrientation();
		
		while(blocks.hasNext()){
			
			countBlocks++;
			SingleInterval block = blocks.next();
			
			if(countBlocks==1){
				if(orientation.equals(Strand.POSITIVE))
					block = new SingleInterval(block.getReferenceName(),block.getReferenceStartPosition()-10000,block.getReferenceEndPosition(),block.getOrientation());
			}
			if(countBlocks==numBlocks){
				if(orientation.equals(Strand.NEGATIVE)) {
					block = new SingleInterval(block.getReferenceName(),block.getReferenceStartPosition(),block.getReferenceEndPosition()+10000,block.getOrientation());
					//System.out.println(block.toBED());
					}
			}
			
			region.addBlocks(block);
		}
		return region;
    }
}

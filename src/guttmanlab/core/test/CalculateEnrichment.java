package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import net.sf.samtools.util.CloseableIterator;

public class CalculateEnrichment {

	public CalculateEnrichment(BAMPairedFragmentCollection bam, Collection<Annotation> regions, String save, Annotation exclude) throws IOException{
		calculateChromosomeCounts(bam, exclude);
		
		Map<Annotation, Double> ratioMap=new TreeMap<Annotation, Double>();
		
		Set<String> excludedSet=new TreeSet<String>();
		CloseableIterator<PairedMappedFragment<SAMFragment>> excludedReads=bam.sortedIterator(exclude, true);
		while(excludedReads.hasNext()){
			PairedMappedFragment<SAMFragment> reads=excludedReads.next();
			excludedSet.add(reads.getName());
		}
		excludedReads.close();
		
		for(Annotation region: regions){
			CloseableIterator<PairedMappedFragment<SAMFragment>> reads=bam.sortedIterator(region, true);
			double numberOfReads=count(reads, excludedSet);
			//double numberOfReads=bam.numOverlappers(region, true);
			int size=region.size();
			double ratio=(double)numberOfReads/(double)size;
			ratioMap.put(region, numberOfReads);
		}
		
		int totalNumberOfReads=bam.getNumAnnotations();
		long genomeLength=bam.getReferenceCoordinateSpace().getTotalReferenceLength();
		double denominator=(double)totalNumberOfReads/genomeLength;
		
		write(save, ratioMap, denominator);
		
	}

	private void calculateChromosomeCounts(BAMPairedFragmentCollection bam, Annotation exclude) {
		int numberOfReads=bam.numOverlappers(exclude, false);
		
		Map<String, Integer> chrLengths=bam.getReferenceCoordinateSpace().getRefSeqLengths();
		for(String chr: chrLengths.keySet()){
			Annotation chrAnnotation=new SingleInterval(chr, 0, chrLengths.get(chr));
			int chrCount=bam.numOverlappers(chrAnnotation, false);
			int chrLength=chrLengths.get(chr);
			double ratio=(double)chrCount/(double)chrLength;
			
			int adjustedCounts=chrCount-numberOfReads;
			int adjustedChrLength=chrLength-exclude.size();
			double adjustedRatio=(double)adjustedCounts/(double)adjustedChrLength;
			
			System.out.println(chr+"\t"+chrCount+"\t"+chrLength+"\t"+ratio+"\t"+adjustedCounts+"\t"+adjustedChrLength+"\t"+adjustedRatio);
		}
	}

	private double count(CloseableIterator<PairedMappedFragment<SAMFragment>> overlappingReads, Set<String> excludedSet) {
		double count=0.0;
		while(overlappingReads.hasNext()){
			PairedMappedFragment<SAMFragment> reads=overlappingReads.next();
			if(!excludedSet.contains(reads.getName())){count++;}
		}
		overlappingReads.close();
		return count;
	}

	private void write(String save, Map<Annotation, Double> ratioMap, double denominator) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Annotation region: ratioMap.keySet()){
			double numberOfReads=ratioMap.get(region);
			double ratio=numberOfReads/region.size();
			double normalizedRatio=ratio/denominator;
			writer.write(region.toUCSC()+"\t"+numberOfReads+"\t"+ratio+"\t"+normalizedRatio+"\n");
		}
		
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
		File bamFile=new File(args[0]);
		AnnotationCollection<Gene> regions=BEDFileIO.loadFromFile(args[1]);
		String save=args[2];
		
		String ucsc=args[3];
		String chr=ucsc.split(":")[0];
		int start=new Integer(ucsc.split(":")[1].split("-")[0]);
		int end=new Integer(ucsc.split(":")[1].split("-")[1]);
		Annotation exclude=new SingleInterval(chr, start, end);
		
		Collection<Annotation> list=new TreeSet<Annotation>();
		Iterator<Gene> iter=regions.sortedIterator();
		while(iter.hasNext()){
			Annotation region=iter.next();
			list.add(region);
		}
		
		BAMPairedFragmentCollection bam=new BAMPairedFragmentCollection(bamFile);
		
		new CalculateEnrichment(bam, list, save, exclude);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=bed file \n args[2]=save \n args[3]=region to exclude";
}

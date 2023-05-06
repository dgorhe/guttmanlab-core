package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.predicate.MinimumLengthFilter;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import net.sf.samtools.util.CloseableIterator;

public class AggregatePlot {
	int windowSize=100;
	int minGeneLength=20000;
	
	
	public AggregatePlot(BAMSingleReadCollection bam, AnnotationCollection<Gene> genes, String save, int minGeneLength, int binSize, boolean spliced) throws IOException{
		this.minGeneLength=minGeneLength;
		
		Map<Integer, Pair<Integer>> distanceCounts3Prime=new TreeMap<Integer, Pair<Integer>>();
		Map<Integer, Pair<Integer>> distanceCounts5Prime=new TreeMap<Integer, Pair<Integer>>();
		
		Map<Integer, Pair<Integer>> distanceCounts5PrimeUpStream=new TreeMap<Integer, Pair<Integer>>();
		Map<Integer, Pair<Integer>> distanceCounts3PrimeDownStream=new TreeMap<Integer, Pair<Integer>>();
		
		AnnotationCollection<Gene> extend5P=extend(genes, true);
		AnnotationCollection<Gene> extend3P=extend(genes, false);
		
		//go through each read and see if it overlaps a gene, if so, what bin is it in --> augment bin counts
		//genes.addFilter(new MinimumLengthFilter<Gene>(minGeneLength));
		CloseableIterator<SAMFragment> reads=bam.sortedIterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMFragment read=reads.next();
			//if overlaps gene, compute distance from 5' end and distance from 3' end
			distance(read, genes, distanceCounts3Prime, true, spliced);
			distance(read, genes, distanceCounts5Prime, false, spliced);
			
			//TODO Get upstream/downstream
			distance(read, extend5P, distanceCounts5PrimeUpStream, false, spliced);
			distance(read, extend3P, distanceCounts3PrimeDownStream, false, spliced);
			
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		
		
		Map<Integer, Pair<List<Integer>>> binned=bin(distanceCounts5Prime, distanceCounts3Prime, distanceCounts5PrimeUpStream, distanceCounts3PrimeDownStream, binSize);
		
		write(save, binned);
		
		/*write(save+".3PrimeUS.distance", distanceCounts3PrimeDownStream, binSize);
		write(save+".5PrimeUS.distance", distanceCounts5PrimeUpStream, binSize);
		
		write(save+".3Prime.distance", distanceCounts3Prime, binSize);
		write(save+".5Prime.distance", distanceCounts5Prime, binSize);*/
	}
	
	

	private void write(String save, Map<Integer, Pair<List<Integer>>> binned) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Integer position: binned.keySet()){
			Pair<List<Integer>> vals=binned.get(position);
			int posVal=sum(vals.getValue1());
			int negVal=sum(vals.getValue2());
			writer.write(position+"\t"+posVal+"\t"+negVal+"\n");
		}
		
		writer.close();
	}



	private int sum(List<Integer> list) {
		int sum=0;
		
		for(Integer v: list){sum+=v;}
		
		return sum;
	}



	private Map<Integer, Pair<List<Integer>>> bin(Map<Integer, Pair<Integer>> dist5Prime, Map<Integer, Pair<Integer>> dist3Prime, Map<Integer, Pair<Integer>> dist5PUS, Map<Integer, Pair<Integer>> dist3PDS, int binSize) {
		Map<Integer, Pair<List<Integer>>> rtrn=new TreeMap<Integer, Pair<List<Integer>>>();
		
		Map<Integer, Pair<List<Integer>>> binned5P=bin(dist5Prime, binSize);
		Map<Integer, Pair<List<Integer>>> binned3P=bin(dist3Prime, binSize);
		Map<Integer, Pair<List<Integer>>> binned3PDS=bin(dist3PDS, binSize);
		Map<Integer, Pair<List<Integer>>> binned5PUS=bin(dist5PUS, binSize);
		
		
		//Merge and put distance in order
		return merge(binned5P, binned3P, binned5PUS, binned3PDS);
		
	}



	private Map<Integer, Pair<List<Integer>>> merge(Map<Integer, Pair<List<Integer>>> binned5p,
			Map<Integer, Pair<List<Integer>>> binned3p, Map<Integer, Pair<List<Integer>>> binned5pus,
			Map<Integer, Pair<List<Integer>>> binned3pds) {
		
		Map<Integer, Pair<List<Integer>>> rtrn=new TreeMap<Integer, Pair<List<Integer>>>();
		
		int counter=0;
		for(Integer position: binned5p.keySet()){
			rtrn.put(position, binned5p.get(position));
			counter=position;
		}
		
		counter=counter+2;
		
		//Need to invert order
		int size=binned3p.size();
		for(Integer position: binned3p.keySet()){
			int relativePosition=(counter+(size-position));
			rtrn.put(relativePosition, binned3p.get(position));
			//counter=Math.max(counter, relativePosition);
		}
		
		counter=counter+size+1;
		
		//counter++;
		
		size=binned3pds.size();
		for(Integer position: binned3pds.keySet()){
			int relativePosition=(counter+position);
			rtrn.put(relativePosition, binned3pds.get(position));
			//counter=Math.max(counter, relativePosition);
		}
		
		counter=counter+size+1;
		
		//Make negative
		for(Integer position: binned5pus.keySet()){
			int relativePosition=-1*(position+1);
			rtrn.put(relativePosition, binned5pus.get(position));
		}
		
		return rtrn;
	}



	private Map<Integer, Pair<List<Integer>>> bin(Map<Integer, Pair<Integer>> distanceMap, int binSize) {
		
		Map<Integer, Pair<List<Integer>>> rtrn=new TreeMap<Integer, Pair<List<Integer>>>();
		
		List<Integer> posCounts=new ArrayList<Integer>();
		List<Integer> negCounts=new ArrayList<Integer>();
		
		int currentPosition=0;
		
		for(int i=0; i<(this.minGeneLength/2); i++){
			int bin=i/binSize;
			if(currentPosition!=bin){
				currentPosition=bin;
				rtrn.put(currentPosition, new Pair<List<Integer>>(posCounts, negCounts));
				posCounts=new ArrayList<Integer>();
				negCounts=new ArrayList<Integer>();
			}
			else{
				if(distanceMap.containsKey(i)){
					posCounts.add(distanceMap.get(i).getValue1());
					negCounts.add(distanceMap.get(i).getValue2());
				}
			}
		}
		
		return rtrn;
	}



	private AnnotationCollection<Gene> extend(AnnotationCollection<Gene> genes, boolean extend5P) {
		FeatureCollection<Gene> rtrn=new FeatureCollection<Gene>();
		
		CloseableIterator<Gene> iter=genes.sortedIterator();
		while(iter.hasNext()){
			Gene gene=iter.next();
			Gene extension=gene.extend5Prime(this.minGeneLength);
			if(!extend5P){
				extension=gene.extend3Prime(this.minGeneLength);
			}
			
			boolean overlapsGene=overlapsGene(extension, genes);
			if(!overlapsGene){
				rtrn.add(extension);
			}
		}
		
		return rtrn;
	}

	private boolean overlapsGene(Gene extension, AnnotationCollection<Gene> genes) {
		Iterator<Gene> overlappers=genes.sortedIterator(extension, false);
		if(overlappers.hasNext()){return true;}
		/*while(overlappers.hasNext()){
			Gene gene=overlappers.next();
			if(gene.getOrientation().equals(extension.getOrientation())){return true;}
		}*/
		return false;
	}



	private void distance(SAMFragment read, AnnotationCollection<Gene> genes,  Map<Integer, Pair<Integer>> distanceCounts, boolean is3Prime, boolean spliced) {
		CloseableIterator<Gene> overlappingGenes=genes.sortedIterator(read, false);
		while(overlappingGenes.hasNext()){
			Gene gene=overlappingGenes.next();
			
			boolean informativeGene=overlaps(read, gene, spliced); //check size, overlap
			boolean strandMatch=gene.getOrientation().equals(read.getOrientation());
			
			if(informativeGene){
				//TODO If spliced, we need to compute distance relative to spliced mRNA
				int distance=distance(read, gene, is3Prime, spliced);
				//System.err.println(read.toUCSC(read.getOrientation())+" "+gene.toUCSC(gene.getOrientation())+" "+distance);
				if(distance<= (minGeneLength/2)){
					Pair<Integer> count=new Pair<Integer>(0,0);
					if(distanceCounts.containsKey(distance)){count=distanceCounts.get(distance);}
					if(strandMatch){
						int posCount=count.getValue1()+1;
						count.setValue1(posCount);
					}
					else if(!strandMatch){
						int negCount=count.getValue2()+1;
						count.setValue2(negCount);
					}
					
					//System.err.println(read.toUCSC()+" "+gene.getName()+" "+distance+" "+is3Prime+" "+count);
					distanceCounts.put(distance, count);
				}
			}
		}
		overlappingGenes.close();
	}



	private boolean overlaps(SAMFragment read, Gene gene, boolean spliced) {
		int geneLength=gene.size();
		if(!spliced){geneLength=gene.getGenomicLength();}
		
		boolean lengthCutoff=geneLength>=this.minGeneLength;
		
		if(!spliced){return lengthCutoff;}
		
		//TODO If spliced, we need to ensure that read overlaps and exon
		boolean overlapsExon=gene.overlaps(read);
		return overlapsExon && lengthCutoff;
	}



	//Given a gene plot 5' and 3' and upstream and downstream by strand
	/*public AggregatePlot(File bamFile, Collection<Gene> regions, String save) throws IOException{
		BAMSingleReadCollection bam1=new BAMSingleReadCollection(bamFile);
		Map<Integer, Integer> updownStreamCounter=new TreeMap<Integer, Integer>();
		
		//Go through each gene and get all reads +/- extend
		int counter=0;
		for(Gene gene: regions){
			//get Upstream and downstream
			Map<Integer, Integer> positionCounts=getUpDownstreamWindows(gene, bam1, updownStreamCounter);
			
			//get first kbs
			//get last kbs
			counter++;
			if(counter%100 ==0){System.err.println(counter+" "+gene.getName()+" "+gene.toUCSC());}
		}
		
		write(save, updownStreamCounter);
	}*/

	private int distance(SAMFragment read, Gene gene, boolean is3Prime, boolean spliced) {
		
		if(!spliced){
			return distanceUnspliced(read, gene, is3Prime);
		}
		
		return distanceSpliced(read, gene, is3Prime);
		
	}



	private int distanceSpliced(SAMFragment read, Gene gene, boolean is3Prime) {
		int distance=gene.distanceFrom3PrimeEnd(read.get3PrimePosition());
		if(!is3Prime){
			distance=gene.distanceFrom5PrimeSpliceSite(read.get5PrimePosition());
		}
		return distance;
	}



	private int distanceUnspliced(SAMFragment read, Gene gene, boolean is3Prime) {
		if(is3Prime){
			int distance=(gene.get3PrimePosition()-read.get3PrimePosition());
			if(gene.getOrientation().equals(Strand.NEGATIVE)){
				distance=read.get3PrimePosition()-gene.get3PrimePosition();
			}
			return distance;
		}
		
		int distance=Math.abs(read.get5PrimePosition()-gene.get5PrimePosition());
		if(gene.getOrientation().equals(Strand.NEGATIVE)){
			distance=gene.get5PrimePosition()-read.get5PrimePosition();
		}
		
		return distance;
	}



	private void write(String save, Map<Integer, Pair<Integer>> distanceMap, int binSize) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int posSum=0;
		int negSum=0;
		int currentPosition=0;
		for(int i=0; i<(this.minGeneLength/2); i++){
			int bin=i/binSize;
			int posCount=0;
			int negCount=0;
			if(distanceMap.containsKey(i)){
				posCount=distanceMap.get(i).getValue1();
				negCount=distanceMap.get(i).getValue2();
			}
			if(currentPosition!=bin){
				writer.write(currentPosition+"\t"+posSum+"\t"+negSum+"\n");
				currentPosition=bin;
				posSum=0;
				negSum=0;
			}
			else{
				posSum+=posCount;
				negSum+=negCount;
			}
		}
		
		writer.close();
	}
	

	public static void main (String[] args) throws IOException{
		if(args.length>5){
			File file=new File(args[0]);
			AnnotationCollection<Gene> regions=BEDFileIO.loadFromFile(args[1]);
			String save=args[2];
			int minGeneLength=new Integer(args[3]);
			int binSize=new Integer(args[4]);
			boolean spliced=new Boolean(args[5]);
			new AggregatePlot(new BAMSingleReadCollection(file), regions, save, minGeneLength, binSize, spliced);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=Annotation BED \n args[2]=save \n args[3]=min gene length \n args[4]=binSize \n args[5]=use spliced cDNA?";
	
}

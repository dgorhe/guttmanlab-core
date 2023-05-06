package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import net.sf.samtools.util.CloseableIterator;

public class Plot5PrimeReadPileUp {

	/*public Plot5PrimeReadPileUp(BAMSingleReadCollection reads, String save, Collection<SingleInterval> regions) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: regions){
			System.err.println("Computing for "+region.toUCSC());
			Map<Integer, Integer> positionCount=countForRegion(region, reads);
			writeBEDGraph(positionCount, region, writer);
		}
		
		writer.close();
		
		
		
	}*/
	
	
	/*public Plot5PrimeReadPileUp(BAMSingleReadCollection reads, String save, Collection<Gene> genes, boolean useSpliced) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		//CloseableIterator<DerivedAnnotation<? extends Annotation>> regions=gene.getGenomicWindows(windowSize, windowSize).sortedIterator();
		
		for(Gene gene: genes){
			if(!useSpliced){
				gene=gene.getGenomicRegion();
			}
			//DerivedAnnotation<? extends Annotation> region=regions.next();
			//int count=reads.numOverlappers(region, true);
			//double expected=(double)count/(double)windowSize;
			//System.err.println("Computing for "+region.toUCSC()+" "+count+" "+expected);
			Map<Integer, Integer> positionCount=countForRegion(gene, reads);
			writeBEDGraph(positionCount, gene, writer);
		}
		
		writer.close();
		
		
		
	}*/

	
	public Plot5PrimeReadPileUp(BAMSingleReadCollection reads, String save) throws IOException{
		
		Map<String, Integer> positionCountPositive=new TreeMap<String, Integer>();
		Map<String, Integer> positionCountNegative=new TreeMap<String, Integer>();
		CloseableIterator<SAMFragment> iter=reads.sortedIterator();
		int counter=0;
		while(iter.hasNext()){
			SAMFragment read=iter.next();
			
			
			if(!read.getReferenceName().equals("*") && read.getSamRecord().getSecondOfPairFlag()){
				Map<String, Integer> positionCount=positionCountNegative;
				int position=read.getReferenceEndPosition();
				if(read.getOrientation().equals(Strand.POSITIVE)){
					position=read.getReferenceStartPosition();
					positionCount=positionCountPositive;
				}
				
				int count=0;
				String region=read.getReferenceName()+"\t"+position+"\t"+(position+1);
				
				
				if(positionCount.containsKey(region)){count=positionCount.get(region);}
				count++;
				positionCount.put(region, count);
			}
			counter++;
			if(counter%100000 ==0 ){System.err.println(counter);}
		}
		iter.close();
		write(save+".pos.bedgraph", positionCountPositive);
		write(save+".neg.bedgraph", positionCountNegative);
	}
	
	private void write(String save, Map<String, Integer> positionCount) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String position: positionCount.keySet()){
			writer.write(position+"\t"+positionCount.get(position)+"\n");
		}
		
		writer.close();
	}


	private Map<Integer, Integer> countForRegion(Annotation region, BAMSingleReadCollection reads) {
		Map<Integer, Integer> positionCount=new TreeMap<Integer, Integer>();
		
		//First base position of read 2
		CloseableIterator<SAMFragment> iter;
		if(region==null){
			iter=reads.sortedIterator();
		}
		else {iter=reads.sortedIterator(region, true);}
		
		while(iter.hasNext()){
			SAMFragment read=iter.next();
			//if(!read.getReferenceName().equals("*") && read.getSamRecord().getFirstOfPairFlag()){
			if(!read.getReferenceName().equals("*") && read.getSamRecord().getSecondOfPairFlag()){
				int position=read.getReferenceEndPosition();
				if(read.getOrientation().equals(Strand.POSITIVE)){
					position=read.getReferenceStartPosition();
				}
				
				//SingleInterval interval=new SingleInterval(read.getReferenceName(), position, position+1);
				
				
				
				int count=0;
				if(positionCount.containsKey(position)){count=positionCount.get(position);}
				count++;
				positionCount.put(position, count);
			}
		}
		
		iter.close();
		return positionCount;
	}

	private void writeBEDGraph(Map<Integer, Integer> positionCount, Annotation region, FileWriter writer) throws IOException {
		double total=avgCounts(positionCount);
		
		System.err.println(total);
		
		for(int i=region.getReferenceStartPosition(); i<region.getReferenceEndPosition(); i++){
			int count=0;
			if(positionCount.containsKey(i)){
				count=positionCount.get(i);
			}
			double norm=(double)count/total;
			writer.write(region.getReferenceName()+"\t"+(i)+"\t"+(i+1)+"\t"+norm+"\n");
		}
		
		
	}
	
	private void writeBEDGraph(Map<Integer, Integer> positionCount, DerivedAnnotation<? extends Annotation> region, double expected, FileWriter writer) throws IOException {
		double total=avgCounts(positionCount);
		
		System.err.println(region.toUCSC()+" "+total+" "+expected);
		
		for(int i=region.getReferenceStartPosition(); i<region.getReferenceEndPosition(); i++){
			int count=0;
			if(positionCount.containsKey(i)){
				count=positionCount.get(i);
			}
			double norm=(double)count/expected;
			writer.write(region.getReferenceName()+"\t"+(i)+"\t"+(i+1)+"\t"+count+"\n");
		}
		
		
	}
	
	private double avgCounts(Map<Integer, Integer> positionCount) {
		double total=0.0;
		double counter=0;
		for(Integer dist: positionCount.keySet()){
			int count=positionCount.get(dist);
			total+=count;
			counter++;
		}
		return total/counter;
	}
	
	private double totalCounts(Map<Integer, Integer> positionCount) {
		double total=0.0;
		for(Integer dist: positionCount.keySet()){
			int count=positionCount.get(dist);
			total+=count;
		}
		return total;
	}

	private static SingleInterval parse(String region) {
		String chr=region.split(":")[0];
		int start=new Integer(region.split(":")[1].split("-")[0]);
		int end=new Integer(region.split(":")[1].split("-")[1]);
		return new SingleInterval(chr, start, end);
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>1){
			BAMSingleReadCollection bam1=new BAMSingleReadCollection(new File(args[0]));
			String save=args[1];
			
			new Plot5PrimeReadPileUp(bam1, save);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=bam file \n args[1]=BED file (genes) \n args[2]=save (.bedgraph) \n args[3]=use spliced (optional, default=true)";
	
}

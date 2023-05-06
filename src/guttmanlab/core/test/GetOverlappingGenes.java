package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import net.sf.samtools.util.CloseableIterator;

public class GetOverlappingGenes {
	int minSize=50;
	String fastaFile="/groups/guttman/genomes/mm10/mm10withchr.fa";

	public GetOverlappingGenes(Collection<SingleInterval> regions, Map<String, IntervalTree<Gene>> trees, String save, Map<String, Sequence> sequences) throws IOException{
		
		Collection<Gene> rtrn=new TreeSet<Gene>();
		Collection<SingleInterval> introns=new TreeSet<SingleInterval>();
		Collection<SingleInterval> exons=new TreeSet<SingleInterval>();
		
		
		for(SingleInterval region: regions){
			Iterator<Node<Gene>> iter= trees.get(region.getReferenceName()).overlappers(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			while(iter.hasNext()){
				Gene gene=iter.next().getValue();
				if(gene.hasCodingRegion()){
					rtrn.add(gene);
					introns.addAll(gene.getIntronSet());
					exons.add(gene.get5PrimeExon());
				}
			}
		}
		
		
		System.err.println(introns.size());
		
		exons=minExon(exons);
		write(save, exons);
		
		write(save+".fa", exons, sequences);
		
	}
	
	

	private void write(String string, Collection<SingleInterval> exons, Map<String, Sequence> sequences) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(SingleInterval exon: exons){
			System.err.println(exon.toUCSC());
			Sequence seq=sequences.get(exon.getReferenceName()).getSubsequence(exon);
			seq.setName(exon.getName()+"_"+exon.toUCSC()+"_"+exon.getOrientation());
			writer.write(seq.toFasta()+"\n");
		}
		
		writer.close();
	}



	private Collection<SingleInterval> minExon(Collection<SingleInterval> exons) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		FeatureCollection<SingleInterval> collection=new FeatureCollection<SingleInterval>();
		for(SingleInterval exon: exons){collection.addAnnotation(exon);}
		
		for(SingleInterval exon: exons){
			CloseableIterator<SingleInterval> iter=collection.sortedIterator(exon, false);
			SingleInterval minExon=min(iter);
			if(minExon!=null){rtrn.add(minExon);}
			//else{rtrn.add(exon);}
		}
		
		return rtrn;
	}



	private SingleInterval min(CloseableIterator<SingleInterval> iter) {
		int start=0;
		int end=Integer.MAX_VALUE;
		String chr="";
		Strand strand=Strand.UNKNOWN;
		String name="";
		
		boolean hasVal=false;
		while(iter.hasNext()){
			SingleInterval next=iter.next();
			start=Math.max(start, next.getReferenceStartPosition());
			end=Math.min(end,  next.getReferenceEndPosition());
			chr=next.getReferenceName();
			strand=next.getOrientation();
			name=next.getName();
			hasVal=true;
		}
		int length=end-start;
		if(!hasVal || length<minSize){return null;}
		
		return new SingleInterval(chr, start, end, strand, name);
	}



	private void write(String save, Collection<? extends Annotation> rtrn) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Annotation gene: rtrn){
			writer.write(gene.toBED()+"\n");
		}
		
		writer.close();
		
	}



	public static void main(String[] args) throws IOException{
		
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[0]);
		Map<String, IntervalTree<Gene>> trees=BEDFileIO.loadTree(args[1]);
		String save=args[2];
		Map<String, Sequence> sequences=FastaFileIOImpl.readFromFileByName("/groups/guttman/genomes/mm10/mm10withchr.fa");
		
		
		new GetOverlappingGenes(regions, trees, save, sequences);
	}
	
}

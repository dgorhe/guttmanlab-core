package guttmanlab.core.probegeneration;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class GenerateProbesForSpecificRegion {
	Collection<SingleInterval> probeRegions;
	
	public GenerateProbesForSpecificRegion(SingleInterval region, int probeLength) throws IOException{
		probeRegions=getProbes(region, probeLength);
	}

	private Collection<SingleInterval> getStrand(Collection<SingleInterval> probeRegions2, Map<String, IntervalTree<Gene>> genes) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		//If overlaps gene then make it antisense to gene
		for(SingleInterval probe: probeRegions2){
			Iterator<Gene> iter=genes.get(probe.getReferenceName()).overlappingValueIterator(probe.getReferenceStartPosition(), probe.getReferenceEndPosition());
			Strand strand=getStrand(iter);
			if(!strand.equals(Strand.UNKNOWN)){probe.setOrientation(strand);}
			rtrn.add(probe);
		}
		
		
		return rtrn;
	}

	private Strand getStrand(Iterator<Gene> iter) {
		Collection<Strand> list=new ArrayList<Strand>();
		while(iter.hasNext()){
			Strand s=iter.next().getOrientation();
			if(!list.contains(s)){list.add(s);}
		}
		if(list.size()==1){return list.iterator().next();}
		return Strand.UNKNOWN;
	}

	private Collection<SingleInterval> filterRepeats(Map<String, IntervalTree<SingleInterval>> repeats, Collection<SingleInterval> probeRegions2) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval probe: probeRegions2){
			boolean hasRepeat=repeats.get(probe.getReferenceName()).hasOverlappers(probe.getReferenceStartPosition(), probe.getReferenceEndPosition());
			if(!hasRepeat){rtrn.add(probe);}
		}
		
		System.err.println(probeRegions2.size()+" "+rtrn.size());
		return rtrn;
	}

	
	

	private Collection<SingleInterval> filterProbes(Map<String, IntervalTree<Gene>> genes, Collection<SingleInterval> probeRegions2) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval probe: probeRegions2){
			boolean hasOverlapper=genes.get(probe.getReferenceName()).hasOverlappers(probe.getReferenceStartPosition(), probe.getReferenceEndPosition());
			if(hasOverlapper){rtrn.add(probe);}
		}
		return rtrn;
	}

	private SingleInterval getFullRegion(Collection<SingleInterval> regions) {
		SingleInterval first=getFirst(regions);
		SingleInterval last=getLast(regions);
		return new SingleInterval(first.getReferenceName(), first.getReferenceStartPosition(), last.getReferenceEndPosition());
	}

	private SingleInterval getLast(Collection<SingleInterval> regions) {
		SingleInterval rtrn=null;
		for(SingleInterval region: regions){
			rtrn=region;
		}
		return rtrn;
	}

	private SingleInterval getFirst(Collection<SingleInterval> regions) {
		return regions.iterator().next();
	}
	

	private static void write(String string, Collection<SingleInterval> regions, Map<String, Sequence> sequenceByChr, String primer1, String primer2, String name) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		int counter=0;
		for(SingleInterval region: regions){
			Sequence probeSeq=sequenceByChr.get(region.getReferenceName()).getSubsequence(region);
			boolean rejectProbe=rejectProbe(probeSeq.getSequenceBases());
			if(!rejectProbe){
				String fullProbe=primer1+probeSeq.getSequenceBases()+Sequence.reverseComplement(primer2);
				writer.write(name+"\t"+region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+region.getOrientation()+"\t"+probeSeq.getSequenceBases()+"\t"+fullProbe+"\n");	
			}
			else{
				System.err.println(counter+" "+probeSeq.getSequenceBases());
			}
				counter++;
		}
		writer.close();
	}
	
	
	
	
	private void writeBED(String string, Collection<Gene> genesToUse) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Gene gene: genesToUse){
			writer.write(gene.toBED()+"\n");
		}
		
		writer.close();
	}

	private Collection<Gene> getGenes(Collection<SingleInterval> regions, Map<String, IntervalTree<Gene>> geneTree) {
		Collection<Gene> genes=new TreeSet<Gene>();
		
		for(SingleInterval region: regions){
			Iterator<Gene> iter=geneTree.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			while(iter.hasNext()){genes.add(iter.next());}
		}
		
		return genes;
	}

	/*private void write(String save, Collection<Sequence> probes) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Sequence probe: probes){
			writer.write(probe.getName()+"\t"+probe.getSequenceBases()+"\n");
		}
		
		writer.close();
	}*/

	private Collection<Sequence> getProbes(Collection<Gene> genesToUse, Map<String, Sequence> sequenceByChr, int probeLength, int spacing) {
		Collection<Sequence> rtrn=new TreeSet<Sequence>();
		int counter=0;
		for(Gene gene: genesToUse){
			System.err.println(gene.getName()+" "+rtrn.size()+" "+counter+" "+genesToUse.size());
			Collection<Sequence> sequences=getProbes(gene, sequenceByChr, probeLength, spacing);
			rtrn.addAll(sequences);
			counter++;
		}
		return rtrn;
	}
	
	
	
	private Collection<SingleInterval> getProbes(Collection<Gene> genesToUse, int probeLength, int spacing) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		int counter=0;
		for(Gene gene: genesToUse){
			System.err.println(gene.getName()+" "+rtrn.size()+" "+counter+" "+genesToUse.size());
			Collection<SingleInterval> sequences=getProbes(gene, probeLength, spacing);
			rtrn.addAll(sequences);
			counter++;
		}
		return rtrn;
	}
	
	private Collection<SingleInterval> getProbes(SingleInterval region, int probeLength) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		net.sf.samtools.util.CloseableIterator<DerivedAnnotation<? extends Annotation>> iter=region.getWindows(probeLength, probeLength).sortedIterator();
		
		Strand orientation=Strand.POSITIVE;
		while(iter.hasNext()){
			DerivedAnnotation<? extends Annotation> window=iter.next();
			SingleInterval interval=new SingleInterval(window.getReferenceName(), window.getReferenceStartPosition(), window.getReferenceEndPosition());
			interval.setOrientation(orientation);
			rtrn.add(interval);
			//orientation=updateStrand(orientation);
		}
		
		return rtrn;
		
	}
	
	private Strand updateStrand(Strand orientation) {
		if(orientation.equals(Strand.POSITIVE)){return Strand.NEGATIVE;}
		return Strand.POSITIVE;
	}

	private Collection<SingleInterval> getProbes(Gene gene, int probeLength, int spacing) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		net.sf.samtools.util.CloseableIterator<DerivedAnnotation<? extends Annotation>> iter=gene.getGenomicRegion().getWindows(probeLength, spacing).sortedIterator();
		
		while(iter.hasNext()){
			DerivedAnnotation<? extends Annotation> window=iter.next();
			SingleInterval interval=new SingleInterval(window.getReferenceName(), window.getReferenceStartPosition(), window.getReferenceEndPosition());
			rtrn.add(interval);
		}
		
		return rtrn;
		
	}

	private Collection<Sequence> getProbes(Gene gene, Map<String, Sequence> sequenceByChr, int probeLength, int spacing) {
		Collection<Sequence> rtrn=new TreeSet<Sequence>();
		net.sf.samtools.util.CloseableIterator<DerivedAnnotation<? extends Annotation>> iter=gene.getGenomicRegion().getWindows(probeLength, spacing).sortedIterator();
		
		while(iter.hasNext()){
			DerivedAnnotation<? extends Annotation> window=iter.next();
			window.setName(window.toUCSC());
			Sequence probeSeq=sequenceByChr.get(window.getReferenceName()).getSubsequence(window);
			if(gene.getOrientation().equals(Strand.POSITIVE)){
				probeSeq=probeSeq.reverseComplement();
			}
			boolean rejectProbe=rejectProbe(probeSeq.getSequenceBases());
			if(!rejectProbe){
				rtrn.add(probeSeq);
			}
		}
		
		return rtrn;
		
	}
	
	private static boolean rejectProbe(String probe) {
		LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 15, 12);
		//RepeatFilter rf=new RepeatFilter(0.07, true, true);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe);
	}
	
	private static Collection<Pair<String>> getPrimers(String file) throws IOException {
		Collection<Pair<String>> rtrn=new ArrayList<Pair<String>>();
		
		Collection<String> lines=BEDFileIO.loadLines(file);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			Pair<String> pair=new Pair<String>(tokens[1], tokens[2]);
			pair.setName(tokens[0]);
			rtrn.add(pair);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>5){
			SingleInterval region=new SingleInterval(args[0]);
			int probeLength=Integer.parseInt(args[1]);
			String save=args[2];
			Map<String, Sequence> sequenceByChr=FastaFileIOImpl.readFromFileByName(args[3]);
			String primer1=args[4];
			String primer2=args[5];
			
			GenerateProbesForSpecificRegion g=new GenerateProbesForSpecificRegion(region,probeLength);
			Collection<SingleInterval> probes=g.probeRegions;
					
			
			write(save, probes, sequenceByChr, primer1, primer2, region.toUCSC());
			//writer.close();
		}
		else{System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=region \n args[1]=probe length \n args[2]=save \n args[3]=fasta file \n args[4]=primer 1 \n args[5]=primer 2";
}

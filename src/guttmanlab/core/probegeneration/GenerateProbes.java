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

public class GenerateProbes {
	Collection<SingleInterval> probeRegions;
	int rnaFragmentLength=25;

	public GenerateProbes(Collection<SingleInterval> regions, Map<String, IntervalTree<Gene>> genes, int probeLength, int spacing, VCFFileReader vcf, String sp1, String sp2, Map<String, IntervalTree<SingleInterval>> repeats, int fragmentLength) throws IOException{
		SingleInterval region=getFullRegion(regions);
		probeRegions=getProbes(region, probeLength, spacing);
		this.rnaFragmentLength=fragmentLength;
		
		//probeRegions=filterProbes(genes, probeRegions);
		
		probeRegions=filterProbes(vcf, sp1, sp2, region, probeRegions);
		
		//TODO Ensure strand is aligned with gene
		//TODO Extend gene 3' in same strand and 5' in antisense strand
		//TODO filter against known repeats
		probeRegions=filterRepeats(repeats, probeRegions);
		
		probeRegions=getStrand(probeRegions, genes);
		
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

	private Collection<SingleInterval> filterProbes(VCFFileReader vcf, String sp1, String sp2, SingleInterval region, Collection<SingleInterval> probeRegions2) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		IntervalTree<SingleInterval> variantTree=getVariants(vcf, sp1, sp2, region);
		
		
		for(SingleInterval probe: probeRegions2){
			boolean overlapsSNP=variantTree.hasOverlappers(probe.getReferenceStartPosition(), probe.getReferenceEndPosition());
			if(overlapsSNP){rtrn.add(probe);}
		}
		
		System.err.println(probeRegions2.size()+" "+rtrn.size());
		return rtrn;
	}

	private IntervalTree<SingleInterval> getVariants(VCFFileReader vcf, String sp1, String sp2, SingleInterval region) {
		CloseableIterator<VariantContext> iter=vcf.query(region.getReferenceName().replace("chr", ""), region.getReferenceStartPosition(), region.getReferenceEndPosition());
		
		IntervalTree<SingleInterval> variantTree=new IntervalTree<SingleInterval>();
		
		while(iter.hasNext()){
			VariantContext variant=iter.next();
			boolean isHet=isHetero(variant, sp1, sp2);
			if(isHet){
				//generate probe
				String chr="chr"+variant.getChr();
				int start=variant.getStart()-rnaFragmentLength;
				int end=variant.getEnd()+rnaFragmentLength;
				SingleInterval variantRegion=new SingleInterval(chr, start, end);
				variantTree.put(start, end, variantRegion);
			}
		}
		
		iter.close();
		return variantTree;
	}

	private boolean isHetero(VariantContext variant, String sp1, String sp2) {
		return !variant.getGenotype(sp1).sameGenotype(variant.getGenotype(sp2));
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
	

	private static void write(String string, Map<String, Collection<SingleInterval>> regions, Map<String, Sequence> sequenceByChr, Collection<Pair<String>> primers) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		Iterator<Pair<String>> primerIter=primers.iterator();
		
		//Pair<String> allPrimer=primerIter.next();
		int counter=0;
		for(String name: regions.keySet()){
			Pair<String> specificPrimer=primerIter.next();
			for(SingleInterval region: regions.get(name)){
				////String seq=sequenceByChr.get(region.getReferenceName()).getSubsequence(region).getSequenceBases();
				Sequence probeSeq=sequenceByChr.get(region.getReferenceName()).getSubsequence(region);
				//boolean rejectProbe=false;
				boolean rejectProbe=rejectProbe(probeSeq.getSequenceBases());
				if(!rejectProbe){
					String fullProbe=specificPrimer.getValue1()+probeSeq.getSequenceBases()+Sequence.reverseComplement(specificPrimer.getValue2());
					writer.write(name+"\t"+region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+region.getOrientation()+"\t"+probeSeq.getSequenceBases()+"\t"+specificPrimer.getName()+"\t"+fullProbe+"\n");
					
				}
				else{
					System.err.println(counter+" "+probeSeq.getSequenceBases());
				}
				counter++;
			}
			
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
	
	private Collection<SingleInterval> getProbes(SingleInterval region, int probeLength, int spacing) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		net.sf.samtools.util.CloseableIterator<DerivedAnnotation<? extends Annotation>> iter=region.getWindows(probeLength, spacing).sortedIterator();
		
		Strand orientation=Strand.POSITIVE;
		while(iter.hasNext()){
			DerivedAnnotation<? extends Annotation> window=iter.next();
			SingleInterval interval=new SingleInterval(window.getReferenceName(), window.getReferenceStartPosition(), window.getReferenceEndPosition());
			interval.setOrientation(orientation);
			rtrn.add(interval);
			orientation=updateStrand(orientation);
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
		if(args.length>9){
			File[] files=new File(args[0]).listFiles();
			Map<String, IntervalTree<Gene>> geneTree=BEDFileIO.loadTree(args[1]);
			String save=args[3];
			int probeLength=new Integer(args[4]);
			int spacing=new Integer(args[5]);
			VCFFileReader vcf=new VCFFileReader(new File(args[6]));
			String sp1=args[7];
			String sp2=args[8];
			Map<String, IntervalTree<SingleInterval>> repeats=BEDFileIO.loadRepeatTree(args[9]);
			Collection<Pair<String>> primers=getPrimers(args[10]);
			int fragmentLength=new Integer(args[11]);
			
			Map<String, Collection<SingleInterval>> allProbes=new TreeMap<String, Collection<SingleInterval>>();
			
			//FileWriter writer=new FileWriter(save);
			int sum=0;
			for(int i=0; i<files.length; i++){
				System.err.println(files[i].getAbsolutePath());
				Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(files[i].getAbsolutePath());
				if(!regions.isEmpty()){
					GenerateProbes g=new GenerateProbes(regions, geneTree,probeLength, spacing, vcf, sp1, sp2, repeats, fragmentLength);
					Collection<SingleInterval> probes=g.probeRegions;
					allProbes.put(files[i].getName(), probes);
					sum+=probes.size();
				}
			}
			System.err.println(allProbes.size()+" "+sum);
			repeats=null;
			System.gc();
			Map<String, Sequence> sequenceByChr=FastaFileIOImpl.readFromFileByName(args[2]);
			write(save, allProbes, sequenceByChr, primers);
			//writer.close();
		}
		else{System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=regions \n args[1]=Genes \n args[2]=Fasta file \n args[3]=save \n args[4]=probe length \n args[5]=probe spacing \n args[6]=vcf file \n args[7]=species 1 \n args[8]=species 2 \n args[9]]=repeats \n args[10]=primers \n args[11]=fragment length";
}

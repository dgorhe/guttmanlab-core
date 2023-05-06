package guttmanlab.core.sars;

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
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.CloseableIterator;

public class SplicingRatios {

	Map<Gene, Pair<Integer>> counts;
	
	//TODO Filter by 3' splice site
	//TODO Remove fragments that are in multiple non-overlapping genes (long inserts)
	//TODO Remove fragments where exon/intron is ambiguous across annotations
	
	public SplicingRatios(File bamFile, Collection<Gene> genes, String save) throws IOException{
		Map<String, IntervalTree<Gene>> geneTree=makeTree(genes);
		
		//Identify reads that span splice junctions
		BAMPairedFragmentCollection reads=readsOverlappingSpliceJunction(bamFile, geneTree);
		
		//Classify each read as spliced, unspliced, ambiguous
		this.counts=assignSpliceStates(reads, geneTree, save);
		
	}
	
	public SplicingRatios(File bamFile, Collection<Gene> genes) throws IOException{
		Map<String, IntervalTree<Gene>> geneTree=makeTree(genes);
		
		//Identify reads that span splice junctions
		BAMPairedFragmentCollection reads=readsOverlappingSpliceJunction(bamFile, geneTree);
		
		//Classify each read as spliced, unspliced, ambiguous
		this.counts=assignSpliceStates(reads, geneTree);
		
	}

	private Map<Gene, Pair<Integer>> assignSpliceStates(BAMPairedFragmentCollection bamData, Map<String, IntervalTree<Gene>> geneTree, String save) throws IOException {
		Map<Gene, Integer> spliced=new TreeMap<Gene, Integer>();
		Map<Gene, Integer> unspliced=new TreeMap<Gene, Integer>();
		Map<Gene, Integer> ambiguos=new TreeMap<Gene, Integer>();
		
		
		SAMFileWriter ambWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(bamData.getFileHeader(), false, new File(save+".amb.bam"));
		SAMFileWriter exonWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(bamData.getFileHeader(), false, new File(save+".exon.bam"));
		SAMFileWriter intronWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(bamData.getFileHeader(), false, new File(save+".intron.bam"));
		
		
		CloseableIterator<PairedMappedFragment<SAMFragment>> reads=bamData.sortedIterator();
		
		Collection<Gene> allGenes=new TreeSet<Gene>();
		
		int counter=0;
		while(reads.hasNext()){
			PairedMappedFragment<SAMFragment> fragment=reads.next();
			Collection<Gene> genes=findGene(fragment, geneTree);
			Collection<String> states=new TreeSet<String>();
			for(Gene gene: genes){
				allGenes.add(gene);
				String state=assignState(fragment, gene);
				states.add(state);
				add(gene, state, spliced, unspliced, ambiguos);
				//write(fragment, state, exonWriter, intronWriter, ambWriter);
			}
			String consensusState=getConsensus(states);
			//add(genes, consensusState, spliced, unspliced, ambiguos);
			write(fragment, consensusState, exonWriter, intronWriter, ambWriter);
			
			//if(states.size()>1){write(fragment, "amb", exonWriter, intronWriter, ambWriter);}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		reads.close();
		
		//write(save, spliced, unspliced);
		
		ambWriter.close();
		exonWriter.close();
		intronWriter.close();
		
		Map<Gene, Pair<Integer>> rtrn=new TreeMap<Gene, Pair<Integer>>();
		for(Gene gene: allGenes){
			Pair<Integer> pair=new Pair<Integer>(get(spliced, gene), get(unspliced,gene));
			rtrn.put(gene, pair);
		}
		return rtrn;
	}
	
	private Map<Gene, Pair<Integer>> assignSpliceStates(BAMPairedFragmentCollection bamData, Map<String, IntervalTree<Gene>> geneTree) throws IOException {
		Map<Gene, Integer> spliced=new TreeMap<Gene, Integer>();
		Map<Gene, Integer> unspliced=new TreeMap<Gene, Integer>();
		Map<Gene, Integer> ambiguos=new TreeMap<Gene, Integer>();
		
		
		CloseableIterator<PairedMappedFragment<SAMFragment>> reads=bamData.sortedIterator();
		
		Collection<Gene> allGenes=new TreeSet<Gene>();
		
		int counter=0;
		while(reads.hasNext()){
			PairedMappedFragment<SAMFragment> fragment=reads.next();
			Collection<Gene> genes=findGene(fragment, geneTree);
			Collection<String> states=new TreeSet<String>();
			for(Gene gene: genes){
				allGenes.add(gene);
				String state=assignState(fragment, gene);
				states.add(state);
				add(gene, state, spliced, unspliced, ambiguos);
				//write(fragment, state, exonWriter, intronWriter, ambWriter);
			}
			String consensusState=getConsensus(states);
			//add(genes, consensusState, spliced, unspliced, ambiguos);
			//write(fragment, consensusState, exonWriter, intronWriter, ambWriter);
			
			//if(states.size()>1){write(fragment, "amb", exonWriter, intronWriter, ambWriter);}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		reads.close();
		
		//write(save, spliced, unspliced);
		
		
		
		Map<Gene, Pair<Integer>> rtrn=new TreeMap<Gene, Pair<Integer>>();
		for(Gene gene: allGenes){
			Pair<Integer> pair=new Pair<Integer>(get(spliced, gene), get(unspliced,gene));
			rtrn.put(gene, pair);
		}
		return rtrn;
	}

	private String getConsensus(Collection<String> states) {
		if(states.size()>1){return "amb";}
		return states.iterator().next();
	}

	private void write(PairedMappedFragment<SAMFragment> fragment, String state, SAMFileWriter exonWriter, SAMFileWriter intronWriter, SAMFileWriter ambWriter) {
		SAMFileWriter writerToUse=null;
		
		if(state.equals("spliced")){writerToUse=exonWriter;}
		if(state.equals("unspliced")){writerToUse=intronWriter;}
		if(state.equals("amb")){writerToUse=ambWriter;}
		
		writerToUse.addAlignment(fragment.getRead1().getSamRecord());
		writerToUse.addAlignment(fragment.getRead2().getSamRecord());
		
	}

	/*private void write(String save, Map<Gene, Integer> spliced, Map<Gene, Integer> unspliced) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Gene> allGenes=new TreeSet<Gene>();
		allGenes.addAll(spliced.keySet());
		allGenes.addAll(unspliced.keySet());
		
		for(Gene gene: allGenes){
			int splicedCount=get(spliced, gene);
			int unsplicedCount=get(unspliced, gene);
			double ratio=(double)unsplicedCount/(double)(splicedCount+unsplicedCount);
			writer.write(gene.getName()+"\t"+splicedCount+"\t"+unsplicedCount+"\t"+ratio+"\n");
		}
		
		writer.close();
	}*/
	
	private int get(Map<Gene, Integer> map, Gene gene) {
		int count=0;
		if(map.containsKey(gene)){count=map.get(gene);}
		return count;
	}

	private String assignState(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		SAMFragment read1=fragment.getRead1();
		SAMFragment read2=fragment.getRead2();
		
		for(Annotation intron: gene.getIntrons()){
			if(read1.overlaps(intron) || read2.overlaps(intron)){return "unspliced";}
		}
		
		return "spliced"; //TODO Add ambiguous
		
		//if(state.equals("amb")){mapToUse=ambiguos;}
	}

	private void add(Collection<Gene> genes, String state, Map<Gene, Integer> spliced, Map<Gene, Integer> unspliced, Map<Gene, Integer> ambiguos) {
		for(Gene gene: genes){
			add(gene, state, spliced, unspliced, ambiguos);
		}
	}
	
	private void add(Gene gene, String state, Map<Gene, Integer> spliced, Map<Gene, Integer> unspliced, Map<Gene, Integer> ambiguos) {
		int count=0;
		Map<Gene, Integer> mapToUse=null;
		
		if(state.equals("spliced")){mapToUse=spliced;}
		if(state.equals("unspliced")){mapToUse=unspliced;}
		if(state.equals("amb")){mapToUse=ambiguos;}
		
		if(mapToUse.containsKey(gene)){count=mapToUse.get(gene);}
		count++;
		mapToUse.put(gene, count);
	}

	private BAMPairedFragmentCollection readsOverlappingSpliceJunction(File bamFile, Map<String, IntervalTree<Gene>> geneTree) {
		BAMPairedFragmentCollection bamData=new BAMPairedFragmentCollection(bamFile);
		File file=new File(bamFile.getAbsolutePath()+".junctionReads.bam");
		
		SAMFileHeader header=bamData.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter alignmentWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, file);
		
		
		CloseableIterator<PairedMappedFragment<SAMFragment>> reads=bamData.sortedIterator();
		
		int counter=0;
		while(reads.hasNext()){
			PairedMappedFragment<SAMFragment> fragment=reads.next();
			Collection<Gene> genes=findGene(fragment, geneTree);
			//Does fragment overlap junction
			boolean overlapsJunction=overlapsJunction(fragment, genes);
			
			if(overlapsJunction){
				alignmentWriter.addAlignment(fragment.getRead1().getSamRecord());
				alignmentWriter.addAlignment(fragment.getRead2().getSamRecord());
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		reads.close();
		alignmentWriter.close();
		
		return new BAMPairedFragmentCollection(file);
	}
	
	
	
	/*private boolean overlapsJunction(PairedMappedFragment<SAMFragment> fragment, Collection<Gene> genes) {
		for(Gene gene: genes){
			//if fragment overlaps two exons or an exon and intron
			for(Annotation intron: gene.getIntrons()){
				boolean overlaps=(fragment.contains(intron.getReferenceStartPosition()) || fragment.contains(intron.getReferenceEndPosition()));
				if(overlaps){return true;}
			}
		}
		return false;
		
	}*/
	
	private boolean overlapsJunction(PairedMappedFragment<SAMFragment> fragment, Collection<Gene> genes) {
		for(Gene gene: genes){
			//if fragment overlaps two exons or an exon and intron
			for(Annotation intron: gene.getIntrons()){
				//Overlaps only if contains 3' splice site
				int start=intron.get3PrimePosition();
				//boolean overlaps=(fragment.contains(intron.getReferenceStartPosition()) || fragment.contains(intron.getReferenceEndPosition()));
				boolean overlaps=fragment.contains(start);
				if(overlaps){return true;}
			}
		}
		return false;
		
	}

	private Map<String, IntervalTree<Gene>> makeTree(Collection<Gene> genes) {
		Map<String, IntervalTree<Gene>> tree=new TreeMap<String, IntervalTree<Gene>>();
		
		for(Gene gene: genes){
			IntervalTree<Gene> temp=new IntervalTree<Gene>();
			if(tree.containsKey(gene.getReferenceName())){
				temp=tree.get(gene.getReferenceName());
			}
			temp.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
			tree.put(gene.getReferenceName(), temp);
		}
		return tree;
	}
	
	private Collection<Gene> findGene(Annotation fragment, Map<String, IntervalTree<Gene>> geneTree) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		IntervalTree<Gene> tree=geneTree.get(fragment.getReferenceName());
		if(tree!=null){
			Iterator<Gene> genes=tree.overlappingValueIterator(fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition());
			while(genes.hasNext()){
				Gene gene=genes.next();
				if(fragment.getOrientation().equals(gene.getOrientation())){
					if(gene.getNumberOfBlocks()>1){
						rtrn.add(gene);
					}
				}
			}
		}
		return rtrn;
		
	}
	
	private static void write(String save, Map<Gene, Pair<Integer>> sample1, Map<Gene, Pair<Integer>> sample2) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Gene> allGenes=new TreeSet<Gene>();
		allGenes.addAll(sample1.keySet());
		allGenes.addAll(sample2.keySet());
		
		for(Gene gene: allGenes){
			Pair<Integer> count1=sample1.get(gene);
			Pair<Integer> count2=sample2.get(gene);
			
			if(count1==null){count1=new Pair<Integer>(0,0);}
			if(count2==null){count2=new Pair<Integer>(0,0);}
			
			double ratio1=(double)count1.getValue1()/(double)(count1.getValue1()+count1.getValue2());
			double ratio2=(double)count2.getValue1()/(double)(count2.getValue1()+count2.getValue2());
			writer.write(gene.getName()+"\t"+gene.getGenomicLength()+"\t"+count1.getValue1()+"\t"+count1.getValue2()+"\t"+count2.getValue1()+"\t"+count2.getValue2()+"\t"+ratio1+"\t"+ratio2+"\n");
		}
		
		writer.close();
	}
	
	private static void write(String save, Map<Gene, Pair<Integer>> sample1, Collection<Gene> genes) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Gene> allGenes=new TreeSet<Gene>();
		allGenes.addAll(sample1.keySet());
		
		
		for(Gene gene: genes){
			Pair<Integer> count1=sample1.get(gene);
			if(count1==null){
				count1=new Pair<Integer>(0,0);
			}
			double ratio=(double)count1.getValue1()/(double)(count1.getValue1()+count1.getValue2());
			
			writer.write(gene.getName()+"\t"+gene.getGenomicLength()+"\t"+count1.getValue1()+"\t"+count1.getValue2()+"\t"+ratio+"\n");
		}
		
		writer.close();
		
	}
	
	public Map<Gene, Pair<Integer>> getCounts(){return this.counts;}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File bam1=new File(args[0]);
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile((args[1]));
			String save=args[2];
			SplicingRatios sample1=new SplicingRatios(bam1, genes, save);
			write(save, sample1.counts, genes);
		}
		else{System.err.println(usage);}
	}
	
	


	static String usage=" args[0]=bam file 1 \n args[1]=Genes (BED file) \n args[2]=save";
	
	
}

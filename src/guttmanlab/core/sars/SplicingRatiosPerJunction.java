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
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.datastructures.Pair;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.CloseableIterator;

public class SplicingRatiosPerJunction {

	Map<Gene, Pair<Integer>> counts;
	Map<Gene, Collection<SAMFragment>> readsByJunction;
	int binSize=100000;
	
	//TODO Filter by 3' splice site
	//TODO Remove fragments that are in multiple non-overlapping genes (long inserts)
	//TODO Remove fragments where exon/intron is ambiguous across annotations
	
	public SplicingRatiosPerJunction(File bamFile, Collection<Gene> genes, String save) throws IOException{
		Collection<Gene> junctions=BEDFileIO.getJunctions(genes);
		junctions=makeUnique(junctions);
		
		BEDFileIO.writeGeneBED(junctions, save+"Juncs.bed");		
		
		System.err.println("Done writing unique junctions");
		
		Map<String, IntervalTree<Gene>> geneTree=makeTree(junctions);
		
		this.readsByJunction=new TreeMap<Gene, Collection<SAMFragment>>();
		//SAMFileHeader header=getHeader(bamFile);
		
		//Identify reads that span splice junctions
		readsOverlappingSpliceJunction(bamFile, geneTree);
		
		
		//writeReadsByJunction(save, header);
		
		Map<Gene, Pair<Integer>> countsByJunction=assignReadsByJunction();
		write(save, countsByJunction);
		writeBedgraph(save+".splicingefficiency.bedgraph", countsByJunction);
	}
	
	private void writeBedgraph(String save, Map<Gene, Pair<Integer>> countsByJunction) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<SingleInterval, Pair<Integer>> binnedVals=bin(countsByJunction, binSize);
		
		for(SingleInterval r: binnedVals.keySet()) {
			Pair<Integer> pair=binnedVals.get(r);
			double sum=pair.getValue1()+pair.getValue2();
			if(sum>100) {
				double se=(double)pair.getValue1()/sum;
				writer.write(r.toBedgraph(se)+"\n");
			}
		}
		
		
		/*for(Gene junction: countsByJunction.keySet()) {
			Pair<Integer> pair=countsByJunction.get(junction);
			double sum=pair.getValue1()+pair.getValue2();
			double splicingefficiency=(double)pair.getValue1()/sum;
			if(sum>25) {
				Annotation intron=junction.getIntrons().iterator().next();
				int ss=intron.get3PrimePosition();
				int start=ss-50;
				int end=ss+50;
				writer.write(junction.getReferenceName()+"\t"+start+"\t"+end+"\t"+splicingefficiency+"\n");
			}
		}*/
		
		writer.close();
		
	}

	private Map<SingleInterval, Pair<Integer>> bin(Map<Gene, Pair<Integer>> countsByJunction, int binSize2) {
		Map<SingleInterval, Pair<Integer>> rtrn=new TreeMap<SingleInterval, Pair<Integer>>();
		
		for(Gene g: countsByJunction.keySet()) {
			SingleInterval bin=g.getGenomicRegion().bin(binSize2);
			if(!rtrn.containsKey(bin)) {rtrn.put(bin, new Pair<Integer>(0,0));}
			Pair<Integer> pair=rtrn.get(bin);
			pair.setValue1(countsByJunction.get(g).getValue1()+pair.getValue1());
			pair.setValue2(countsByJunction.get(g).getValue2()+pair.getValue2());
		}
		
		return rtrn;
	}

	private Collection<Gene> makeUnique(Collection<Gene> junctions) {
		Map<String, IntervalTree<Gene>> trees=new TreeMap<String, IntervalTree<Gene>>();
		
		for(Gene junction: junctions) {
			//System.err.println(junction.getName());
			String chr=junction.getReferenceName();
			if(!trees.containsKey(chr)) {trees.put(chr, new IntervalTree<Gene>());}
			IntervalTree<Gene> tree=trees.get(chr);
			Annotation intron=junction.getIntrons().iterator().next();
			tree.put(intron.getReferenceStartPosition(), intron.getReferenceEndPosition(), junction);
		}
		
		//iterate through junctions and find overlappers. If ==1 add to list;
		
		
		Collection<Gene> rtrn=new TreeSet<Gene>();
		for(Gene junction: junctions) {
			//System.err.println(junction.getName());
			boolean use=getOverlappers(junction, trees);
			if(use) {rtrn.add(junction);}
		}
		
		return rtrn;
	}

	private boolean getOverlappers(Gene junction, Map<String, IntervalTree<Gene>> trees) {
		Annotation intron=junction.getIntrons().iterator().next();
		if(trees.containsKey(junction.getReferenceName())) {
			IntervalTree<Gene> tree=trees.get(junction.getReferenceName());
			Iterator<Node<Gene>> overlappers=tree.overlappers(intron.getReferenceStartPosition(), intron.getReferenceEndPosition());
			int counter=0;
			while(overlappers.hasNext()) {
				overlappers.next();
				counter++;
				if(counter>1) {return false;}
			}
			
		}
		return true;
	}

	private void write(String save, Map<Gene, Pair<Integer>> countsByJunction) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene junction: countsByJunction.keySet()) {
			Pair<Integer> pair=countsByJunction.get(junction);
			writer.write(junction.getName()+"\t"+junction.toUCSC()+"\t"+pair.getValue1()+"\t"+pair.getValue2()+"\n");
		}
		
		writer.close();
	}

	private Map<Gene, Pair<Integer>> assignReadsByJunction(String save, SAMFileHeader header) {
		Map<Gene, Pair<Integer>> rtrn=new TreeMap<Gene, Pair<Integer>>();
		for(Gene junction: this.readsByJunction.keySet()) {
			Pair<Integer> pair=assignReads(this.readsByJunction.get(junction), junction, save, header);
			rtrn.put(junction, pair);
		}
		return rtrn;
	}
	
	private Map<Gene, Pair<Integer>> assignReadsByJunction() {
		Map<Gene, Pair<Integer>> rtrn=new TreeMap<Gene, Pair<Integer>>();
		for(Gene junction: this.readsByJunction.keySet()) {
			Pair<Integer> pair=assignReads(this.readsByJunction.get(junction), junction);
			rtrn.put(junction, pair);
		}
		return rtrn;
	}

	private Pair<Integer> assignReads(Collection<SAMFragment> reads, Gene junction) {
		int exonCount=0;
		int intronCount=0;
		
		for(SAMFragment fragment: reads) {
			String state=assignState(fragment, junction);
			if(state.equals("spliced")) {exonCount++;}
			if(state.equals("unspliced")) {intronCount++;}
		}
		
		Pair<Integer> rtrn=new Pair<Integer>(exonCount, intronCount);
		
		return rtrn;
	}
	
	private Pair<Integer> assignReads(Collection<SAMFragment> reads, Gene junction, String saveBase, SAMFileHeader header) {
		String save=saveBase+"."+junction.getName();
		SAMFileWriter exonWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".exon.bam"));
		SAMFileWriter intronWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".intron.bam"));
		
		int exonCount=0;
		int intronCount=0;
		
		for(SAMFragment fragment: reads) {
			String state=assignState(fragment, junction);
			write(fragment, state, exonWriter, intronWriter);
			if(state.equals("spliced")) {exonCount++;}
			if(state.equals("unspliced")) {intronCount++;}
		}
		
		Pair<Integer> rtrn=new Pair<Integer>(exonCount, intronCount);
		
		exonWriter.close();
		intronWriter.close();
		return rtrn;
	}

	private SAMFileHeader getHeader(File bamFile) {
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMFileHeader rtrn=reader.getFileHeader();
		reader.close();
		return rtrn;
	}

	private void writeReadsByJunction(String saveBase, SAMFileHeader header) {
		for(Gene junction: this.readsByJunction.keySet()) {
			String save=saveBase+"."+junction.getName()+".bam";
			writeBam(save, this.readsByJunction.get(junction), header);
		}
		
	}

	private void writeBam(String save, Collection<SAMFragment> collection, SAMFileHeader header) {
		File file=new File(save);
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter alignmentWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, file);
		
		for(SAMFragment fragment: collection){
			alignmentWriter.addAlignment(fragment.getSamRecord());
		}
		alignmentWriter.close();
		
	}

	

	/*private Map<Gene, Pair<Integer>> assignSpliceStates(BAMPairedFragmentCollection bamData, Map<String, IntervalTree<Gene>> geneTree, String save) throws IOException {
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
	}*/
	
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
		//if(states.size()>1) {return "spliced";}
		return states.iterator().next();
	}

	private void write(SAMFragment read, String state, SAMFileWriter exonWriter, SAMFileWriter intronWriter) {
		SAMFileWriter writerToUse=null;
		
		if(state.equals("spliced")){writerToUse=exonWriter;}
		if(state.equals("unspliced")){writerToUse=intronWriter;}
		
		writerToUse.addAlignment(read.getSamRecord());
		
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
	
	
	private String assignState(SAMFragment read, Gene gene) {
		Annotation intron=gene.getIntrons().iterator().next();
		SingleInterval exon1=gene.getFirstBlock();
		SingleInterval exon2=gene.getLastBlock();
		
		if(read.overlaps(intron)) {return "unspliced";}
		else if(read.overlaps(exon1) && read.overlaps(exon2)){return "spliced";}
		else{return "amb";}
	}

	private String assignState(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		SAMFragment read1=fragment.getRead1();
		SAMFragment read2=fragment.getRead2();
		
		for(Annotation intron: gene.getIntrons()){
			if(read1.overlaps(intron) || read2.overlaps(intron)){return "unspliced";}
		}
		
		if(overlapsExon(fragment, gene)) {return "spliced";}
		
		return "amb"; //TODO Add ambiguous
		
		//if(state.equals("amb")){mapToUse=ambiguos;}
	}

	private boolean overlapsExon(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		SingleInterval exon1=gene.getFirstBlock();
		SingleInterval exon2=gene.getLastBlock();
		
		if(fragment.getRead1().overlaps(exon1) && fragment.getRead2().overlaps(exon2)) {
			return true;
		}
		if(fragment.getRead1().overlaps(exon2) && fragment.getRead2().overlaps(exon1)) {
			return true;
		}
		return false;
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

	private void readsOverlappingSpliceJunction(File bamFile, Map<String, IntervalTree<Gene>> geneTree) {
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		File file=new File("/Users/mguttman/Desktop/test.junctionReads.bam");
		
		SAMFileHeader header=reader.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter alignmentWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, file);
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment read=new SAMFragment(record);
			
			Collection<Gene> genes=findGene(read, geneTree);
			//Does fragment overlap junction
			//boolean overlapsJunction=overlapsJunction(fragment, genes);
			Collection<Gene> overlappingJunctions=overlapsJunction(read, genes);
			
			if(!overlappingJunctions.isEmpty()){
				alignmentWriter.addAlignment(record);
			}
			
			for(Gene g: overlappingJunctions) {
				if(!this.readsByJunction.containsKey(g)) {this.readsByJunction.put(g, new ArrayList<SAMFragment>());}
				Collection<SAMFragment> list=this.readsByJunction.get(g);	
				list.add(read);
			}
			
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		alignmentWriter.close();
	}
	
	/*Fragments**/
	/*
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
			//boolean overlapsJunction=overlapsJunction(fragment, genes);
			Collection<Gene> overlappingJunctions=overlapsJunction(fragment, genes);
			
			if(!overlappingJunctions.isEmpty()){
				alignmentWriter.addAlignment(fragment.getRead1().getSamRecord());
				alignmentWriter.addAlignment(fragment.getRead2().getSamRecord());
			}
			
			for(Gene g: overlappingJunctions) {
				if(!this.readsByJunction.containsKey(g)) {this.readsByJunction.put(g, new ArrayList<PairedMappedFragment<SAMFragment>>());}
				Collection<PairedMappedFragment<SAMFragment>> list=this.readsByJunction.get(g);	
				list.add(fragment);
			}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		reads.close();
		alignmentWriter.close();
		
		return new BAMPairedFragmentCollection(file);
	}*/
	
	
	private Collection<Gene> overlapsJunction(SAMFragment read, Collection<Gene> genes) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		for(Gene gene: genes){
			for(Annotation intron: gene.getIntrons()){
				//Overlaps only if contains 3' splice site
				int start=intron.get3PrimePosition();
				boolean overlaps=read.contains(start);
				if(overlaps){rtrn.add(gene);}
			}
		}
		return rtrn;
		
	}
	
	private Collection<Gene> overlapsJunction(PairedMappedFragment<SAMFragment> fragment, Collection<Gene> genes) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		for(Gene gene: genes){
			//if fragment overlaps two exons or an exon and intron
			for(Annotation intron: gene.getIntrons()){
				//Overlaps only if contains 3' splice site
				int start=intron.get3PrimePosition();
				//boolean overlaps=(fragment.contains(intron.getReferenceStartPosition()) || fragment.contains(intron.getReferenceEndPosition()));
				boolean overlaps=fragment.contains(start);
				if(overlaps){rtrn.add(gene);}
			}
		}
		return rtrn;
		
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
			SplicingRatiosPerJunction sample1=new SplicingRatiosPerJunction(bam1, genes, save);
			//write(save, sample1.counts, genes);
		}
		else{System.err.println(usage);}
	}
	
	


	static String usage=" args[0]=bam file 1 \n args[1]=Genes (BED file) \n args[2]=save";
	
	
}

package guttmanlab.core.smit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;


import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.util.CloseableIterator;

public class AssignReadsToSplicingStates {
	
	Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> splicedFragments;
	Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> unsplicedFragments;
	Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> ambigousFragments;
	private double intronLengthCutoff;
	private int splicedCount;
	private int unsplicedCount;
	private int ambiguousCount;

	Map<String, IntervalTree<Gene>> geneTree;
	private double percentileCutoff=0.99;
	private int exonCount;
	private int intronCount;
	private int exonLength;
	private int intronLength;
	
	Map<Gene, Integer> splicedCounts;
	Map<Gene, Integer> unsplicedCounts;
	Map<Gene, Integer> ambiguousCounts;
	
	public enum SpliceState {
		SPLICED('s'), UNSPLICED('u'), AMBIGUOUS('a');
		private char value;
		
		private SpliceState(char value) {
			this.value = value;
		}
		
		public String toString() {
			return "" + value;
		}
	}
	
	
	public AssignReadsToSplicingStates(BAMPairedFragmentCollection bamData, AnnotationCollection<Gene> genes, String save) throws IOException{
		geneTree=makeTree(genes);
		splicedCount=0;
		unsplicedCount=0;
		ambiguousCount=0;
		exonCount=0;
		intronCount=0;
		
		//plotPolymerasePosition(bamData, save+".polymerasePosition.bedgraph");
		
		splicedCounts=new TreeMap<Gene, Integer>();
		unsplicedCounts=new TreeMap<Gene, Integer>();
		ambiguousCounts=new TreeMap<Gene, Integer>();
		
		Collection<PairedMappedFragment<SAMFragment>> informativeFragments=new TreeSet<PairedMappedFragment<SAMFragment>>();
		Collection<PairedMappedFragment<SAMFragment>> uninformativeFragments=new TreeSet<PairedMappedFragment<SAMFragment>>();
		splicedFragments=new TreeMap<Integer, Collection<PairedMappedFragment<SAMFragment>>>();
		unsplicedFragments=new TreeMap<Integer, Collection<PairedMappedFragment<SAMFragment>>>();
		ambigousFragments=new TreeMap<Integer, Collection<PairedMappedFragment<SAMFragment>>>();
		
		//TODO Define cutoff for intron length
		//intronLengthCutoff=buildEmpiricalInsertDistribution(bamData);
		intronLengthCutoff=450;
		System.err.println("intron cutoff "+intronLengthCutoff);
		
		//Iterate through all fragments
		CloseableIterator<PairedMappedFragment<SAMFragment>> reads=bamData.sortedIterator();
		
		//Calculate the length of all introns and all exons
		getLengths(geneTree);
		
		int counter=0;
		while(reads.hasNext()){
			PairedMappedFragment<SAMFragment> fragment=reads.next();
			Collection<Gene> overlappingGenes=findGene(fragment);
			//Only analyze if overlaps a single isoform TODO: This should be fixed
			if(overlappingGenes.size()==1){
				Gene gene=overlappingGenes.iterator().next();
				classifyPolymerasePosition(fragment, gene);
				boolean overlapsSS=overlapsSS(fragment, gene);
				if(overlapsSS){
					informativeFragments.add(fragment);
					SpliceState state=assignSpliceState(fragment, gene);
				}
				else{
					//overlapsIntron(fragment);
					uninformativeFragments.add(fragment);
				}
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		reads.close();
		
		
		writeFraction(save+".count");
		
		SAMFileHeader header=bamData.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		writeSAM(new File(save+".splicedfragments.bam"), splicedFragments, header);
		writeSAM(new File(save+".ambfragments.bam"), ambigousFragments, header);
		writeSAM(new File(save+".unsplicedfragments.bam"), unsplicedFragments, header);
		writeSAM(new File(save+".informative.bam"), informativeFragments, header);
		writeSAM(new File(save+".uninformative.bam"), uninformativeFragments, header);
		
		//writeByDistance(save+".dist.bed", splicedFragments, unsplicedFragments);
		
		
		System.err.println("Exon: "+exonCount+"\nIntron: "+intronCount+"\nExon Length "+this.exonLength+"\nIntron length "+this.intronLength);
		System.err.println("Total: "+counter+"\ninformative: "+informativeFragments.size()+"\nspliced: "+splicedCount+"\nunspliced: "+unsplicedCount+"\nambiguous: "+ ambiguousCount);
	}

	private void writeFraction(String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Gene> genes=new TreeSet<Gene>();
		genes.addAll(this.splicedCounts.keySet());
		genes.addAll(this.unsplicedCounts.keySet());
		
		for(Gene gene: genes){
			int spliced=0;
			if(splicedCounts.containsKey(gene)){spliced=splicedCounts.get(gene);}
			int unspliced=0;
			if(unsplicedCounts.containsKey(gene)){unspliced=unsplicedCounts.get(gene);}
			int ambiguous=0;
			if(ambiguousCounts.containsKey(gene)){ambiguous=ambiguousCounts.get(gene);}
			writer.write(gene.getName()+"\t"+spliced+"\t"+unspliced+"\t"+ambiguous+"\n");
		}
		
		writer.close();
	}

	private void getLengths(Map<String, IntervalTree<Gene>> geneTree2) {
		for(String chr: geneTree2.keySet()){
			Iterator<Gene> iter=geneTree2.get(chr).valueIterator();
			while(iter.hasNext()){
				Gene gene=iter.next();
				int totalDistance=gene.getReferenceEndPosition()-gene.getReferenceStartPosition();
				int exonDistance=gene.size();
				int intronDistance=totalDistance-exonDistance;
				this.intronLength+=intronDistance;
				this.exonLength+=exonDistance;
			}
		}
		
	}

	private void classifyPolymerasePosition(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		boolean overlapsExon=gene.overlapsExon(new SingleInterval(gene.getReferenceName(), fragment.get3PrimePosition(), (fragment.get3PrimePosition()+1), gene.getOrientation()));
		if(overlapsExon){exonCount++;}
		else{intronCount++;}
	}

	private void plotPolymerasePosition(BAMPairedFragmentCollection bamData, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		//iterate by chromosome
		for(String chr: bamData.getReferenceCoordinateSpace().getRefSeqLengths().keySet()){
			System.err.println(chr);
			Map<Integer, Integer> positionCounter=new TreeMap<Integer, Integer>();
			SingleInterval chrRegion=new SingleInterval(chr, 0, bamData.getReferenceCoordinateSpace().getRefSeqLengths().get(chr));
			CloseableIterator<PairedMappedFragment<SAMFragment>> reads=bamData.sortedIterator(chrRegion, true);
			while(reads.hasNext()){
				PairedMappedFragment<SAMFragment> read=reads.next();
				int position=read.get3PrimePosition();
				int count=0;
				if(positionCounter.containsKey(position)){
					count=positionCounter.get(position);
				}
				count++;
				positionCounter.put(position, count);
			}
			reads.close();
			write(writer, positionCounter, chr);
		}
		
		writer.close();
	}

	private void write(FileWriter writer, Map<Integer, Integer> positionCounter, String chr) throws IOException {
		for(int position: positionCounter.keySet()){
			writer.write(chr+"\t"+position+"\t"+(position+1)+"\t"+positionCounter.get(position)+"\n");
		}
		
	}

	private void writeByDistance(String save,
			Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> splicedFragments2,
			Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> unsplicedFragments2) throws IOException {
		
		FileWriter writer=new FileWriter(save);
		
		for(Integer distance: splicedFragments2.keySet()){
			for(PairedMappedFragment<SAMFragment> fragment: splicedFragments2.get(distance)){
				writer.write(fragment.getReferenceName()+"\t"+fragment.getReferenceStartPosition()+"\t"+fragment.getReferenceEndPosition()+"\t"+"s="+distance+"\n");
			}
		}
		
		for(Integer distance: unsplicedFragments2.keySet()){
			for(PairedMappedFragment<SAMFragment> fragment: unsplicedFragments2.get(distance)){
				writer.write(fragment.getReferenceName()+"\t"+fragment.getReferenceStartPosition()+"\t"+fragment.getReferenceEndPosition()+"\t"+"u="+distance+"\n");
			}
		}
		
		writer.close();
	}

	private void writeSAM(File save, Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> splicedFragments2, SAMFileHeader header, int start, int end) {
		SAMFileWriter informativeAlignmentWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, save);
		
		for(Integer d: splicedFragments2.keySet()){
			if(d>=start && d<=end){
				for(PairedMappedFragment<SAMFragment> fragment: splicedFragments2.get(d)){
					informativeAlignmentWriter.addAlignment(fragment.getRead1().getSamRecord());
					informativeAlignmentWriter.addAlignment(fragment.getRead2().getSamRecord());
				}
			}
		}
		
		informativeAlignmentWriter.close();
	}

	private void writeSAM(File save, Collection<PairedMappedFragment<SAMFragment>> informativeFragments, SAMFileHeader header) {
		SAMFileWriter informativeAlignmentWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, save);
		
		for(PairedMappedFragment<SAMFragment> fragment: informativeFragments){
			informativeAlignmentWriter.addAlignment(fragment.getRead1().getSamRecord());
			informativeAlignmentWriter.addAlignment(fragment.getRead2().getSamRecord());
		}
		
		informativeAlignmentWriter.close();
	}

	private void writeSAM(File file, Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> fragmentsByDistance, SAMFileHeader header) {
		SAMFileWriter alignmentWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, file);
		writeSAM(alignmentWriter, fragmentsByDistance);
		alignmentWriter.close();
	}
	
	
	private void writeSAM(SAMFileWriter alignmentWriter, Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> map) {
		for(Integer distance: map.keySet()){
			Collection<PairedMappedFragment<SAMFragment>> set=map.get(distance);
			for(PairedMappedFragment<SAMFragment> pair: set){
				alignmentWriter.addAlignment(pair.getRead1().getSamRecord());
				alignmentWriter.addAlignment(pair.getRead2().getSamRecord());
			}
		}
		
		
	}
	
	private boolean onlyOverlapsFirstExon(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		SingleInterval firstExon=gene.getFirstBlock();
		if(gene.getOrientation().equals(Strand.NEGATIVE)){
			firstExon=gene.getLastBlock();
		}
		
		SingleInterval fragmentRegion=new SingleInterval(fragment.getReferenceName(), fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition());
		if(firstExon.overlaps(fragmentRegion)){
			//TODO Does it overlaps elsewhere too?
			Iterator<SingleInterval> exons=gene.getBlocks();
			while(exons.hasNext()){
				SingleInterval exon= exons.next();
				if(exon!=firstExon){
					if(exon.overlaps(fragmentRegion)){return false;}
				}
			}
			return true;
			
		}
		return false;
	}


	private void writeSAM(SAMFileWriter informativeAlignmentWriter, PairedMappedFragment<SAMFragment> fragment) {
		informativeAlignmentWriter.addAlignment(fragment.getRead1().getSamRecord());
		informativeAlignmentWriter.addAlignment(fragment.getRead2().getSamRecord());
	}


	private double buildEmpiricalInsertDistribution(BAMPairedFragmentCollection bamData) {
		System.err.println("Started building");
		//Build an EmpiricalDistribution from read fragments that are within an exon or within an intron
		//Iterate through all fragments
		CloseableIterator<PairedMappedFragment<SAMFragment>> iter=bamData.sortedIterator();
		List<Double> fragmentLength=new ArrayList<Double>();
		int counter=0;
		while(iter.hasNext()){
			PairedMappedFragment<SAMFragment> fragment=iter.next();
			Collection<Gene> overlappingGenes=findGene(fragment);
			
			//Only consider single isoform genes
			if(overlappingGenes.size()==1){
				Gene gene=overlappingGenes.iterator().next();
				//retain fragment if it is entirely within 1 intron or entirely within 1 exon
				//int numExons=gene.getNumberOfOverlappingBlocks(fragment);
				
				int hasOneOverlappingExon=gene.hasOneOverlappingExon(fragment);
				int hasOneOverlappingIntron=gene.hasOneOverlappingIntron(fragment);
				
				if(hasOneOverlappingIntron + hasOneOverlappingExon ==1){
					fragmentLength.add((double)fragment.fragmentLength());
				}
				
			}
			counter++;
			if(counter%100000 ==0){System.err.println("Building insert distribution "+counter);}
		}
		
		double rtrn=getPercentile(fragmentLength, percentileCutoff);
		
		
		iter.close();
		return rtrn;
	}
	
	/*private EmpiricalDistribution buildEmpiricalInsertDistribution(BAMPairedFragmentCollection bamData) {
		System.err.println("Started building");
		//Build an EmpiricalDistribution from read fragments that are within an exon or within an intron
		//Iterate through all fragments
		CloseableIterator<PairedMappedFragment<SAMFragment>> iter=bamData.sortedIterator();
		List<Integer> fragmentLength=new ArrayList<Integer>();
		int counter=0;
		while(iter.hasNext()){
			PairedMappedFragment<SAMFragment> fragment=iter.next();
			Collection<Gene> overlappingGenes=findGene(fragment);
			
			//Only consider single isoform genes
			if(overlappingGenes.size()==1){
				Gene gene=overlappingGenes.iterator().next();
				//retain fragment if it is entirely within 1 intron or entirely within 1 exon
				//int numExons=gene.getNumberOfOverlappingBlocks(fragment);
				
				int hasOneOverlappingExon=gene.hasOneOverlappingExon(fragment);
				int hasOneOverlappingIntron=gene.hasOneOverlappingIntron(fragment);
				
				if(hasOneOverlappingIntron + hasOneOverlappingExon ==1){
					fragmentLength.add(fragment.fragmentLength());
				}
				
			}
			counter++;
			if(counter%100000 ==0){System.err.println("Building insert distribution "+counter);}
		}
		
		double[] values=getValues(fragmentLength);
		EmpiricalDistribution rtrn=new EmpiricalDistribution(fragmentLength.size());
		rtrn.load(values);
		
		iter.close();
		return rtrn;
	}*/

	
	private double getPercentile(List<Double> fragmentLength, double percentile) {
		Collections.sort(fragmentLength);
		return Statistics.quantile(fragmentLength, percentile);
	}

	private double[] getValues(List<Integer> fragmentLength) {
		double[] rtrn=new double[fragmentLength.size()];
		
		int i=0;
		for(Integer val: fragmentLength){
			rtrn[i]=val;
			i++;
		}
		
		return rtrn;
	}


	private boolean overlaps3PrimeSS(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		if(gene.getNumberOfBlocks()>1){
			SingleInterval firstExon=getFirstExon(gene);
			//iterate through blocks and get 3'SS
			Iterator<SingleInterval> exons=gene.getBlocks();
			//Test if fragment contains it
			while(exons.hasNext()){
				SingleInterval exon=exons.next();
				if(exon!=firstExon){
					int ss=exon.get5PrimePosition();
					if(overlap(fragment, ss)){return true;}
				}
			}
		}
		
		return false;
	}
	
	private boolean overlapsSS(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		if(gene.getNumberOfBlocks()>1){
			SingleInterval firstExon=getFirstExon(gene);
			SingleInterval lastExon=getLastExon(gene);
			//iterate through blocks and get 3'SS
			Iterator<SingleInterval> exons=gene.getBlocks();
			//Test if fragment contains it
			while(exons.hasNext()){
				SingleInterval exon=exons.next();
				if(exon!=firstExon){
					int ss1=exon.get5PrimePosition();
					if(overlap(fragment, ss1)){return true;}
				}
				if(exon!=lastExon){
					int ss2=exon.get3PrimePosition();
					if(overlap(fragment, ss2)){return true;}
				}
			}
		}
		
		return false;
	}

	private ArrayList<SpliceSite> get3PrimeSS(Gene gene){
		ArrayList<SpliceSite> rtrn=new ArrayList<SpliceSite>();
		
		if(gene.getNumberOfBlocks()>1){
			Collection<Annotation> introns=gene.getIntrons();
			for(Annotation intron: introns){
				SpliceSite ss=new SpliceSite(intron.get3PrimePosition(), intron);
				rtrn.add(ss);
			}
			
			
			/*SingleInterval firstExon=getFirstExon(gene);
			//iterate through blocks and get 3'SS
			Iterator<SingleInterval> exons=gene.getBlocks();
			//Test if fragment contains it
			while(exons.hasNext()){
				SingleInterval exon=exons.next();
				if(exon!=firstExon){
					int ss=exon.get5PrimePosition();
					rtrn.add(ss);
				}
			}*/
		}
		return rtrn;
	}
	

	private SingleInterval getFirstExon(Gene gene) {
		SingleInterval firstExon=gene.getFirstBlock();
		if(gene.getOrientation().equals(Strand.NEGATIVE)){
			firstExon=gene.getLastBlock();
		}
		return firstExon;
	}
	
	private SingleInterval getLastExon(Gene gene) {
		SingleInterval lastExon=gene.getLastBlock();
		if(gene.getOrientation().equals(Strand.NEGATIVE)){
			lastExon=gene.getFirstBlock();
		}
		return lastExon;
	}


	public AssignReadsToSplicingStates(List<PairedMappedFragment<SAMFragment>> fragments, AnnotationCollection<Gene> genes, String save) throws IOException{
		
		geneTree=makeTree(genes);
		
		splicedFragments=new TreeMap<Integer, Collection<PairedMappedFragment<SAMFragment>>>();
		unsplicedFragments=new TreeMap<Integer, Collection<PairedMappedFragment<SAMFragment>>>();
		
		//Iterate through all fragments
		Iterator<PairedMappedFragment<SAMFragment>> iter=fragments.iterator();
		
		int counter=0;
		while(iter.hasNext()){
			PairedMappedFragment<SAMFragment> fragment=iter.next();
			
			Collection<Gene> overlappingGenes=findGene(fragment);
			
			SpliceState state=SpliceState.AMBIGUOUS;
			
			if(overlappingGenes.size()==1){
				Gene gene=overlappingGenes.iterator().next();
				
				//Only consider fragments overlapping an exon, the rest are ambigous
				boolean overlapsExon=gene.overlapsExon(new SingleInterval(fragment.getReferenceName(), fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition()));
				
				if(overlapsExon){
					//get 3' position
					int position=distanceToSpliceSite(fragment, gene);
					
					//assign to splice state
					state=assignSpliceState(position, fragment, gene);
					
					fragment.setName("d="+position+",s="+state);
				}
				
				
				
				
			}
			//if(state!=SpliceState.AMBIGUOUS){System.err.println(counter+" "+fragments.size()+" "+state);}
			counter++;
			if(counter%10000 ==0){System.err.println(counter+" "+fragments.size());}
		}
		writeFraction(save+".count", splicedFragments, unsplicedFragments);	
		write(save+".splicedfragments.bed", splicedFragments);
		write(save+".unsplicedfragments.bed", unsplicedFragments);
	}
	
	private void write(String save, List<Annotation> fragments) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Annotation fragment: fragments){
			writer.write(fragment.toBED()+"\n");
		}
		
		writer.close();
	}
	
	private void write(String save, Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> fragments) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Integer d: fragments.keySet()){
			for(Annotation fragment: fragments.get(d)){
				writer.write(fragment.toBED()+"\n");
			}
		}
		
		writer.close();
	}

	public Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> getSplicedFragments(){return this.splicedFragments;}
	
	public Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> getUnsplicedFragments(){return this.unsplicedFragments;}
	
	
	/*private SpliceState assignSpliceState(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		int position=distanceToSpliceSite(fragment, gene);
		if(position<0){
			//add(position, fragment, this.ambigousFragments);
			return SpliceState.AMBIGUOUS;
		}
		
		//Trim fragment back to last exon
		
		
		if(position>0){
		
			//if read spans exon-exon junction --> spliced
			if(fragment.getRead1().getNumberOfBlocks()>1 || fragment.getRead2().getNumberOfBlocks()>1){
				add(position, fragment, this.splicedFragments);	
				return SpliceState.SPLICED;
			}
			
			
			
			//if read spans exon-intron junction --> unspliced
			else if(spansExonIntron(gene, fragment)){
				add(position, fragment, this.unsplicedFragments);	
				return SpliceState.UNSPLICED;
			}
			
			//Reads are in different introns
			else if(differentIntrons(gene, fragment)){
				add(position, fragment, this.unsplicedFragments);	
				return SpliceState.UNSPLICED;
			}
			
			//Both reads are in different exons --> Check insert sizes
			else if(differentExons(gene, fragment)){
				//TODO Check insert distribution
				int unsplicedSize=fragment.fragmentLength();
				double p=1-this.insertDistribution.cumulativeProbability(unsplicedSize);
				if(p>0.01){
					add(position, fragment, this.unsplicedFragments);	
					return SpliceState.UNSPLICED;
				}
				else{
					add(position, fragment, this.splicedFragments);	
					return SpliceState.SPLICED;
				}
			}
			
			//Check fragment length of spliced and unspliced and decide based on probability
			else{
				//TODO This could be implemented better to do a Likelihood Ratio Test
				//Compute sizes for spliced and unspliced fragments
				//Look up sizes in insertDistribution
				int unsplicedSize=fragment.fragmentLength();
				double p=1-this.insertDistribution.cumulativeProbability(unsplicedSize);
				if(p<0.01){
					add(position, fragment, this.splicedFragments);	
					return SpliceState.SPLICED;
				}
				else{
					add(position, fragment, this.unsplicedFragments);	
					return SpliceState.UNSPLICED;
				}
			}
			
			
			
		}
		add(position, fragment, this.ambigousFragments);	
		return SpliceState.AMBIGUOUS;
	}*/
	
	private SpliceState assignSpliceState(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		ArrayList<SpliceSite> spliceSites=get3PrimeSS(gene);
		Collection<SpliceState> spliceStates=new ArrayList<SpliceState>();
		
		for(SpliceSite spliceSite: spliceSites){
			//for each 3'SS that fragment overlaps
			//Is this splice site informative? Is distance between exons well seperated for spliced and unspliced?
			boolean informativeDistance=informativeDistance(spliceSite);
			if(overlap(fragment, spliceSite)){
				
				//get distance of 3' position to 3'SS
				int distance=distance(fragment, spliceSite);
				//check if fragment is spliced or unspliced BEFORE 3'SS (exon-exon or exon-intron)
				
				
				
				//trim read to 3'SS -- is remaining read in exon or intron?
				boolean spliced=isSpliced(fragment, spliceSite);
				boolean isUnspliced=isUnspliced(fragment, spliceSite);
				
				if(isUnspliced){
					unsplicedCount++;
					add(this.unsplicedFragments, distance, fragment);
					spliceStates.add(SpliceState.UNSPLICED);
				}
				else if(!isUnspliced && informativeDistance){
					splicedCount++;
					add(this.splicedFragments, distance, fragment);
					spliceStates.add(SpliceState.SPLICED);
				}
				else if(spliced){
					splicedCount++;
					add(this.splicedFragments, distance, fragment);
					spliceStates.add(SpliceState.SPLICED);
				}
				else{
					this.ambiguousCount++;
					add(this.ambigousFragments, distance, fragment);
					spliceStates.add(SpliceState.AMBIGUOUS);
				}
				
				/*if(spliced && informativeDistance){
					splicedCount++;
					add(this.splicedFragments, distance, fragment);
				}
				else if(!spliced){
					unsplicedCount++;
					add(this.unsplicedFragments, distance, fragment);
				}
				else{
					this.ambiguousCount++;
					add(this.ambigousFragments, distance, fragment);
				}*/
			}
		}
		return consensusSpliceState(gene, spliceStates);
	}
	
	private SpliceState consensusSpliceState(Gene gene, Collection<SpliceState> spliceStates) {
		//If has unspliced then unspliced
		boolean unspliced=false;
		boolean spliced=false;
		for(SpliceState state: spliceStates){
			if(state.equals(SpliceState.UNSPLICED)){unspliced=true;}
			if(state.equals(SpliceState.SPLICED)){spliced=true;}
		}
		if(unspliced){
			int count=0;
			if(this.unsplicedCounts.containsKey(gene)){count=this.unsplicedCounts.get(gene);}
			count++;
			this.unsplicedCounts.put(gene, count);
			return SpliceState.UNSPLICED;
		}
		if(spliced){
			int count=0;
			if(this.splicedCounts.containsKey(gene)){count=this.splicedCounts.get(gene);}
			count++;
			this.splicedCounts.put(gene, count);
			return SpliceState.SPLICED;
		}
		
		int count=0;
		if(this.ambiguousCounts.containsKey(gene)){count=this.ambiguousCounts.get(gene);}
		count++;
		this.ambiguousCounts.put(gene, count);
		
		return SpliceState.AMBIGUOUS;
	}

	private boolean isUnspliced(PairedMappedFragment<SAMFragment> fragment, SpliceSite spliceSite) {
		//Check if read actually overlaps intron
		boolean intronOverlap=fragment.getRead1().overlaps(spliceSite.getIntron()) || fragment.getRead2().overlaps(spliceSite.getIntron());
		
		if(intronOverlap){return true;}
		
		return false;
	}

	private boolean informativeDistance(SpliceSite spliceSite) {
		int intronLength=spliceSite.getIntron().getGenomicLength();
		return intronLength>intronLengthCutoff;
	}

	private void add(Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> map, int distance,PairedMappedFragment<SAMFragment> fragment) {
		Collection<PairedMappedFragment<SAMFragment>> list=new TreeSet<PairedMappedFragment<SAMFragment>>();
		if(map.containsKey(distance)){
			list=map.get(distance);
		}
		list.add(fragment);
		map.put(distance, list);
	}


	private boolean isSpliced(PairedMappedFragment<SAMFragment> fragment, SpliceSite spliceSite) {
		/*boolean intronOverlap=fragment.getRead1().overlaps(spliceSite.getIntron()) || fragment.getRead2().overlaps(spliceSite.getIntron());
		
		if(intronOverlap){return false;}
		
		return true;*/
		
		//Check if the read is actually spliced
		if(isUnspliced(fragment, spliceSite)){return false;}
		
		boolean spliced=spliced(fragment.getRead1(), spliceSite.getIntron()) || spliced(fragment.getRead2(), spliceSite.getIntron());
		return spliced;
	}


	


	private boolean spliced(SAMFragment read, Annotation intron1) {
		boolean spliced=false;
		//Check if introns overlap
		Iterator<Annotation> introns=read.getIntrons().iterator();
		while(introns.hasNext()){
			Annotation intron=introns.next();
			if(intron.overlaps(intron1)){spliced=true;}
		}
		return spliced;
	}

	private int distance(PairedMappedFragment<SAMFragment> fragment, SpliceSite spliceSite) {
		if(fragment.getOrientation().equals(Strand.NEGATIVE)){return spliceSite.getSpliceSite()-fragment.getReferenceStartPosition();}
		return fragment.getReferenceEndPosition()-spliceSite.getSpliceSite();
	}


	private boolean overlap(PairedMappedFragment<SAMFragment> fragment, SpliceSite ss){
		return overlap(fragment, ss.getSpliceSite());
	}
	
	private boolean overlap(PairedMappedFragment<SAMFragment> fragment, int ss) {
		if(fragment.getReferenceStartPosition()<ss && fragment.getReferenceEndPosition()>ss){return true;}
		return false;
	}


	private boolean spansExonIntron(Gene gene, PairedMappedFragment<SAMFragment> fragment) {
		int spliceSite=gene.get3PrimeExon(fragment.get3PrimePosition()).get5PrimePosition();
		
		boolean overlapsExon1=gene.getOverlappingBlocks(fragment.getRead1()).hasNext();
		Iterator<SingleInterval> overlappingIntrons1=gene.getOverlappingIntrons(fragment.getRead1());
		boolean overlapsIntron1=false;
		
		
		while(overlappingIntrons1.hasNext()){
			SingleInterval intron=overlappingIntrons1.next();
			//is the overlap before the splice site? If so, make overlaps1 true
			if(gene.getOrientation().equals(Strand.POSITIVE) && intron.getReferenceStartPosition()<spliceSite){overlapsIntron1=true;}
			if(gene.getOrientation().equals(Strand.NEGATIVE) && intron.getReferenceEndPosition()>spliceSite){overlapsIntron1=true;}
		}
		
		
		
		boolean overlapsExon2=gene.getOverlappingBlocks(fragment.getRead2()).hasNext();
		Iterator<SingleInterval> overlappingIntrons2=gene.getOverlappingIntrons(fragment.getRead2());
		
		boolean overlapsIntron2=false;
		
		while(overlappingIntrons2.hasNext()){
			SingleInterval intron=overlappingIntrons2.next();
			//is the overlap before the splice site? If so, make overlaps1 true
			if(gene.getOrientation().equals(Strand.POSITIVE) && intron.getReferenceStartPosition()<spliceSite){overlapsIntron2=true;}
			if(gene.getOrientation().equals(Strand.NEGATIVE) && intron.getReferenceEndPosition()>spliceSite){overlapsIntron2=true;}
		}
		
		if(overlapsExon1 && overlapsIntron1){return true;}
		
		if(overlapsExon2 && overlapsIntron2){return true;}
		
		boolean overlapsExon=overlapsExon1 || overlapsExon2;
		boolean overlapsIntron=overlapsIntron1 || overlapsIntron2;
		
		if(overlapsExon && overlapsIntron){return true;}
		
		return false;
	}


	private boolean differentExons(Gene gene, PairedMappedFragment<SAMFragment> fragment) {
		Iterator<SingleInterval> iter1=gene.getOverlappingBlocks(fragment.getRead1());
		Iterator<SingleInterval> iter2=gene.getOverlappingBlocks(fragment.getRead2());
		
		TreeSet<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		while(iter1.hasNext()){rtrn.add(iter1.next());}
		while(iter2.hasNext()){rtrn.add(iter2.next());}
		return rtrn.size()>1;
	}
	
	private boolean differentIntrons(Gene gene, PairedMappedFragment<SAMFragment> fragment) {
		Iterator<SingleInterval> iter1=gene.getOverlappingIntrons(fragment.getRead1());
		Iterator<SingleInterval> iter2=gene.getOverlappingIntrons(fragment.getRead2());
		
		TreeSet<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		while(iter1.hasNext()){rtrn.add(iter1.next());}
		while(iter2.hasNext()){rtrn.add(iter2.next());}
		return rtrn.size()>1;
	}


	private boolean intronOverlapBefore3SS(Gene gene, PairedMappedFragment<SAMFragment> fragment) {
		int spliceSite=gene.get3PrimeExon(fragment.get3PrimePosition()).get5PrimePosition();
		for(SingleInterval intron: gene.getIntronSet()){
			if(precedes(intron, spliceSite)){
				if(fragment.getRead1().overlaps(intron) || fragment.getRead2().overlaps(intron)){return true;}
			
			}
		}
		return false;
	}


	private boolean precedes(SingleInterval intron, int spliceSite) {
		if(intron.getOrientation().equals(Strand.POSITIVE)){
			return intron.getReferenceStartPosition()<spliceSite;
		}
		return intron.getReferenceEndPosition()>spliceSite;
	}


	//Old
	/*private SpliceState assignSpliceState(int position, PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		if(position<0){return SpliceState.AMBIGUOUS;}
		
		
		//The fragment has to overlap the 3' splice site to be considered
		//if(gene.inIntron(position)){return SpliceState.AMBIGUOUS;}
		
		if(position>0){
			
			//If it has intron in it --> unspliced
			//if fragment is spliced --> spliced
			if(fragment.getRead1().getNumberOfBlocks()>1 || fragment.getRead2().getNumberOfBlocks()>1){
				add(position, fragment, this.splicedFragments);	
				return SpliceState.SPLICED;
			}
			//Both reads are in different exons --> Check insert sizes
			
			
			//both reads are within same exon or same intron --> ambiguous
			
			
			//if fragment is not spliced --> unspliced
			else {
				
				//1 intron only or 1 exon only --> ambiguous
				SpliceState state=overlapsIntron(gene, fragment);
				//IF overlaps intron --> unspliced (or completely in an exon)
				if (state.equals(SpliceState.UNSPLICED)){add(position, fragment, this.unsplicedFragments);}
				return state;
				
			}
		}
		return SpliceState.AMBIGUOUS;
	}*/
	
	private SpliceState assignSpliceState(int position, PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		if(position<0){return SpliceState.AMBIGUOUS;}
		
		/*If 3' position is past exon, then trim back to exon for consideration
		Annotation trimmedFragment=trimFragment(fragment);*/
		
		//The fragment has to overlap the 3' splice site to be considered
		//if(gene.inIntron(position)){return SpliceState.AMBIGUOUS;}
		
		if(position>0){
			
			//If it has intron in it --> unspliced
			//if fragment is spliced --> spliced
			if(fragment.getNumberOfBlocks()>1){
				add(position, fragment, this.splicedFragments);	
				return SpliceState.SPLICED;
			}
			
			//if fragment is not spliced --> unspliced
			else if(fragment.getNumberOfBlocks()==1){
				//1 intron only or 1 exon only --> ambiguous
				SpliceState state=overlapsIntron(gene, fragment);
				//IF overlaps intron --> unspliced (or completely in an exon)
				if (state.equals(SpliceState.UNSPLICED)){add(position, fragment, this.unsplicedFragments);}
				return state;
				
			}
		}
		return SpliceState.AMBIGUOUS;
	}

	
	private boolean containedWithinIntron(Gene gene, Annotation fragment){
		for(SingleInterval intron: gene.getIntronSet()){
			if(fragment.getReferenceStartPosition()>=intron.getReferenceStartPosition() && fragment.getReferenceEndPosition()<=intron.getReferenceEndPosition()){return true;}
		}
		return false;
	}
	
	private boolean containedWithinExon(Gene gene, Annotation fragment){
		for(Annotation exon: gene.getExonSet()){
			if(fragment.getReferenceStartPosition()>=exon.getReferenceStartPosition() && fragment.getReferenceEndPosition()<=exon.getReferenceEndPosition()){return true;}
		}
		return false;
	}
	
	/**
	 * Tests if both reads are contained within the same exon
	 * @param gene
	 * @param fragment
	 * @return
	 */
	private boolean withinExon(Gene gene, PairedMappedFragment<SAMFragment> fragment){
		for(Annotation exon: gene.getExonSet()){
			boolean fullyContained=fragment.getReferenceStartPosition()>=exon.getReferenceStartPosition() && fragment.getReferenceEndPosition()<=exon.getReferenceEndPosition();
			if(fullyContained){return true;}
		}
		return false;
	}
	
	private boolean withinIntron(Gene gene, PairedMappedFragment<SAMFragment> fragment){
		for(SingleInterval intron: gene.getIntronSet()){
			boolean fullyContained=fragment.getReferenceStartPosition()>=intron.getReferenceStartPosition() && fragment.getReferenceEndPosition()<=intron.getReferenceEndPosition();
			if(fullyContained){return true;}
		}
		return false;
	}
	
	private SpliceState overlapsIntron(Gene gene, Annotation fragment){
		for(SingleInterval intron: gene.getIntronSet()){
			boolean endsIn=fragment.getReferenceEndPosition()>=intron.getReferenceStartPosition() && fragment.getReferenceEndPosition()<=intron.getReferenceEndPosition();
			boolean startsIn=fragment.getReferenceStartPosition()>=intron.getReferenceStartPosition() && fragment.getReferenceStartPosition()<=intron.getReferenceEndPosition();
			boolean middleOf=intron.getReferenceStartPosition()>=fragment.getReferenceStartPosition() && intron.getReferenceEndPosition()<=fragment.getReferenceEndPosition();
			boolean fullyContained=fragment.getReferenceStartPosition()>=intron.getReferenceStartPosition() && fragment.getReferenceEndPosition()<=intron.getReferenceEndPosition();
			if(fullyContained){return SpliceState.AMBIGUOUS;}
			if(endsIn || startsIn || middleOf){return SpliceState.UNSPLICED;}
		}
		return SpliceState.AMBIGUOUS;
	}
	
	/*private boolean completelyContainedInAnExon(Gene gene, Annotation fragment) {
		int exonCount=gene.getNumberOfOverlappingBlocks(fragment);
		System.err.println("Exon count: "+exonCount);
		if(exonCount>1){return false;}
		SingleInterval exon=gene.getOverlappingBlocks(fragment).next();
		return fragment.getReferenceStartPosition()>=exon.getReferenceStartPosition() && fragment.getReferenceEndPosition()<=exon.getReferenceEndPosition();
	}*/

	private void add(int position, PairedMappedFragment<SAMFragment> fragment, Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> fragments) {
		Collection<PairedMappedFragment<SAMFragment>> temp=new TreeSet<PairedMappedFragment<SAMFragment>>();
		if(fragments.containsKey(position)){
			temp=fragments.get(position);
		}
		temp.add(fragment);
		fragments.put(position, temp);
	}




	private void writeFraction(String save, Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> splicedFragments, Map<Integer, Collection<PairedMappedFragment<SAMFragment>>> unsplicedFragments) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Distance\tSpliced\tUnspliced\t%Spliced\tSpliced (adjusted)\tUnspliced (adjusted)\tratio (adjusted)\n");
		
		Collection<Integer> positions=new TreeSet<Integer>();
		positions.addAll(splicedFragments.keySet());
		positions.addAll(unsplicedFragments.keySet());
		
		for(Integer position: positions){
			Collection<PairedMappedFragment<SAMFragment>> spliced=splicedFragments.get(position);
			Collection<PairedMappedFragment<SAMFragment>> unspliced=unsplicedFragments.get(position);
			int splicedCount=0;
			int unsplicedCount=0;
			if(spliced!=null){splicedCount=spliced.size();}
			if(unspliced!=null){unsplicedCount=unspliced.size();}
			
			double ratio=(double)splicedCount/(double)(splicedCount+unsplicedCount);
			
			double exonRate=(double)this.exonCount/(double)this.exonLength;
			double intronRate=(double) this.intronCount/(double)this.intronLength;
			
			double splicedAdjusted=splicedCount;
			double unsplicedAdjusted=unsplicedCount*(exonRate/intronRate);
			double ratioAdjusted=splicedAdjusted/(splicedAdjusted+unsplicedAdjusted);
			
			writer.write(position+"\t"+splicedCount+"\t"+unsplicedCount+"\t"+ratio+"\t"+splicedAdjusted+"\t"+unsplicedAdjusted+"\t"+ratioAdjusted+"\n");
		}
		
		writer.close();
	}

	private int distanceToSpliceSite(Annotation fragment, Gene gene) {
		return gene.distanceFrom5PrimeSpliceSite(fragment);
		
		//int positionOfPolymerase=fragment.get3PrimePosition();
		//return gene.distanceFrom5PrimeSpliceSite(positionOfPolymerase);
		
		/*int positionOfPolymerase=fragment.get3PrimePosition();
		int distance=-1;
		Collection<Annotation> exons=gene.getBlockSet();
		for(Annotation exon: exons){
			if(positionOfPolymerase>exon.getReferenceStartPosition() && positionOfPolymerase<exon.getReferenceEndPosition()){
				distance=positionOfPolymerase-exon.getReferenceStartPosition();
				if(gene.getOrientation().equals(Strand.NEGATIVE)){
					distance=exon.getReferenceEndPosition()-positionOfPolymerase;
					}
				}
			}
		return distance;*/
		
	}
		
	
	private Map<String, IntervalTree<Gene>> makeTree(AnnotationCollection<Gene> genes) {
		Map<String, IntervalTree<Gene>> tree=new TreeMap<String, IntervalTree<Gene>>();
		
		CloseableIterator<Gene>iter=genes.sortedIterator();
		while(iter.hasNext()){
			Gene gene=iter.next();
			IntervalTree<Gene> temp=new IntervalTree<Gene>();
			if(tree.containsKey(gene.getReferenceName())){
				temp=tree.get(gene.getReferenceName());
			}
			temp.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
			tree.put(gene.getReferenceName(), temp);
		}
		return tree;
	}


	private Collection<Gene> findGene(Annotation fragment) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		IntervalTree<Gene> tree=geneTree.get(fragment.getReferenceName());
		if(tree!=null){
			Iterator<Gene> genes=tree.overlappingValueIterator(fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition());
			while(genes.hasNext()){
				Gene gene=genes.next();
				if(fragment.getOrientation().equals(gene.getOrientation())){
					rtrn.add(gene);
				}
			}
		}
		return rtrn;
		
	}


	private SpliceState assignRead(Annotation fragment, Gene gene) {
		//If fragment is entirely in intron then ambiguous
		for(SingleInterval intron: gene.getIntronSet()){
			if(fragment.getReferenceStartPosition()>=intron.getReferenceStartPosition() && fragment.getReferenceEndPosition()<=intron.getReferenceEndPosition()){
				return SpliceState.AMBIGUOUS;
			}
		}
		
		//If fragment is entirely a single exon then ambiguous
		for(Annotation exon: gene.getExonSet()){
			if(fragment.getReferenceStartPosition()>=exon.getReferenceStartPosition() && fragment.getReferenceEndPosition()<=exon.getReferenceEndPosition()){
				return SpliceState.AMBIGUOUS;
			}
		}
		
		
		//If fragment overlaps any intron --> Unspliced
		for(SingleInterval intron: gene.getIntronSet()){
			if(fragment.overlaps(intron)){return SpliceState.UNSPLICED;}
		}
		
		//Else is spliced --> Spliced
		return SpliceState.SPLICED;
	}

	public static class SpliceSite{
		int spliceSite;
		Annotation intron;
		
		public SpliceSite(int spliceSite, Annotation intron){
			this.spliceSite=spliceSite;
			this.intron=intron;
		}
		
		public Annotation getIntron(){return intron;}
		public int getSpliceSite(){return spliceSite;}
	}
	

	

	private static List<Annotation> parse(String file) throws IOException {
		List<Annotation> rtrn=new ArrayList<Annotation>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			Annotation annotation=BEDFileIO.parse(nextLine);
			rtrn.add(annotation);
			counter++;
			if(counter%1000000 ==0){System.err.println(counter+" "+nextLine);}
		}
		reader.close();
		return rtrn;
	}
	
	
	
	private static void fractionSpliced(AssignReadsToSplicingStates sample1, AssignReadsToSplicingStates sample2, AnnotationCollection<Gene> genes, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Iterator<Gene> geneIter=genes.sortedIterator();
		while(geneIter.hasNext()){
			Gene gene=geneIter.next();
			
			int spliced1=get(sample1.splicedCounts, gene);
			int unspliced1=get(sample1.unsplicedCounts, gene);
			int amb1=get(sample1.ambiguousCounts, gene);
			
			int spliced2=get(sample2.splicedCounts, gene);
			int unspliced2=get(sample2.unsplicedCounts, gene);
			int amb2=get(sample2.ambiguousCounts, gene);
			
			double ratio1=(double)unspliced1/(double)(spliced1+unspliced1+amb1);
			double ratio2=(double)unspliced2/(double)(spliced2+unspliced2+amb2);
			
			writer.write(gene.getName()+"\t"+spliced1+"\t"+unspliced1+"\t"+amb1+"\t"+spliced2+"\t"+unspliced2+"\t"+amb2+"\t"+ratio1+"\t"+ratio2+"\n");
		}
		
		writer.close();
	}

	private static int get(Map<Gene, Integer> map, Gene gene) {
		int count=0;
		if(map.containsKey(gene)){count=map.get(gene);}
		return count;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			BAMPairedFragmentCollection data=new BAMPairedFragmentCollection(new File(args[0]));
			BAMPairedFragmentCollection data2=new BAMPairedFragmentCollection(new File(args[1]));
			
			AnnotationCollection<Gene> genes=BEDFileIO.loadFromFile((args[2]));
			String save=args[3];
			
			AssignReadsToSplicingStates sample1=new AssignReadsToSplicingStates(data, genes, save+".s1");
			AssignReadsToSplicingStates sample2=new AssignReadsToSplicingStates(data2, genes, save+".s2");
			
			fractionSpliced(sample1, sample2, genes, save);
					
		}
		else{
			System.err.println(usage);
		}
	}

	static String usage=" args[0]=Fragments BAM \n args[1]= Fragments2 \n args[2]=Genes \n args[3]=save";
	
}

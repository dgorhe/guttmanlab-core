package guttmanlab.core.sharp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.ScanStat;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class AssignReads {
	
	static final String exon="exon";
	static final String intron="intron";
	static final String ambiguous="amb";
	static final String intergenic="inter";
	
	int binSize=100;
	
	public AssignReads(File bam, Map<String, IntervalTree<Annotation>> genes, String save) throws IOException {
		//go through reads and assign to exons, introns, or ambiguous
		assignPairs(bam, genes, save);
	}
	
	private void significantWindows(File bam, int binSize, Map<String, Integer> geneSizes, String save) throws IOException {
		Map<String, Integer> geneCounts=new TreeMap<String, Integer>();
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator reads=inputReader.iterator();
		
		Map<SingleInterval, Integer> positionCount=new TreeMap<SingleInterval, Integer>();
		Map<SingleInterval, Collection<String>> binToGene=new TreeMap<SingleInterval, Collection<String>>();
		
		int totalCount=0;
		
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				Collection<String> overlappingGenes=getGenes(read);
				Collection<SingleInterval> allBins=SAMFragment.allBins(read, binSize);
				for(SingleInterval bin: allBins) {
					int score=0;
					if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
					score++;
					positionCount.put(bin, score);
					add(overlappingGenes, geneCounts);
					Collection<String> list=new TreeSet<String>();
					if(binToGene.containsKey(bin)) {list=binToGene.get(bin);}
					list.addAll(overlappingGenes);
					binToGene.put(bin, list);
				}
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount+" "+read.getReferenceName());}
		}
				
		reads.close();
		inputReader.close();
		
		
		
		FileWriter writer=new FileWriter(save);
		for(SingleInterval bin: positionCount.keySet()) {
			Collection<String> overlappingGenes=binToGene.get(bin);
			int score=positionCount.get(bin);
			
			double minEnrich;
			double maxP=1.0;
			for(String g: overlappingGenes) {
				int geneScore=geneCounts.get(g);
				int size=geneSizes.get(g);
				double expected=(double)geneScore/(double)size;
				double observed=(double)score/(double)bin.size();
				double enrichment=observed/expected;
				double pval=ScanStat.getPValue(score, expected, bin.size(), size);
				if(score>20) {
					if(pval<0.01){writer.write(bin.toBedgraph(enrichment)+"\n");}
				}
			}
		}
		writer.close();
	}

	private Map<String, Integer> getIntronSizes(Map<String, IntervalTree<Annotation>> genes) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(String chr: genes.keySet()) {
			Iterator<Annotation> iter=genes.get(chr).valueIterator();
			while(iter.hasNext()) {
				Annotation gene=iter.next();
				int size=gene.getGenomicLength()-gene.size();
				if(size>0) {
					rtrn.put(gene.getName(), size);
				}
			}
		}
		
		
		return rtrn;
	}
	
	private Map<String, Integer> getExonSizes(Map<String, IntervalTree<Annotation>> genes) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(String chr: genes.keySet()) {
			Iterator<Annotation> iter=genes.get(chr).valueIterator();
			while(iter.hasNext()) {
				Annotation gene=iter.next();
				int size=gene.size();
				if(size>0) {
					rtrn.put(gene.getName(), size);
				}
			}
		}
		
		
		return rtrn;
	}

	private Collection<String> getGenes(SAMRecord read) {
		String geneNames=(String)read.getAttribute("XL");
		Collection<String> rtrn=new TreeSet<String>();
		String[] tokens=geneNames.split(",");
		for(int i=0; i<tokens.length; i++) {rtrn.add(tokens[i]);}
		return rtrn;
	}

	
	
	

	
	
	
	private void add(Collection<String> overlappingGenes, Map<String, Integer> geneCounts) {
		for(String gene: overlappingGenes) {
			int score=0;
			if(geneCounts.containsKey(gene)) {
				score=geneCounts.get(gene);
			}
			score+=1;
			geneCounts.put(gene, score);
		}		
	}

	/*private void assign(File bam, Map<String, IntervalTree<Annotation>> genes, String save) {
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator reads=inputReader.iterator();
		
		SAMFileHeader header=inputReader.getFileHeader();
		
		SAMFileWriter exonWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".exon.bam"));
		SAMFileWriter intronWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".intron.bam"));
		SAMFileWriter ambiguousWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".ambiguous.bam"));
		SAMFileWriter interWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".intergenic.bam"));
		
		int totalCount=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			String state=assign(read, genes);
			write(read, state, exonWriter, intronWriter, ambiguousWriter, interWriter);
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		exonWriter.close();
		intronWriter.close();
		ambiguousWriter.close();
		interWriter.close();
	}*/

	
	private void assignPairs(File bam, Map<String, IntervalTree<Annotation>> genes, String save) {
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator reads=inputReader.iterator();
		
		SAMFileHeader header=inputReader.getFileHeader();
		
		SAMFileWriter exonWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".exon.bam"));
		SAMFileWriter intronWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".intron.bam"));
		SAMFileWriter ambiguousWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".ambiguous.bam"));
		SAMFileWriter interWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".intergenic.bam"));
		
		
		Map<String, Pair<SAMRecord>> pairedReads=new TreeMap<String, Pair<SAMRecord>>();
		
		
		int totalCount=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			
			String name=read.getReadName();
			Pair<SAMRecord> pair=new Pair<SAMRecord>();
			if(pairedReads.containsKey(name)) {pair=pairedReads.get(name);}
			pair=update(pair, read);
			if(complete(pair)) {
				pairedReads.remove(name);
				Pair<String> state=assign(pair, genes);
				write(pair, state, exonWriter, intronWriter, ambiguousWriter, interWriter);
			}
			else {
				pairedReads.put(name, pair);
			}
			
			
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		exonWriter.close();
		intronWriter.close();
		ambiguousWriter.close();
		interWriter.close();
	}
	
	
	

	private boolean complete(Pair<SAMRecord> pair) {
		return pair.isComplete();
	}

	private Pair<SAMRecord> update(Pair<SAMRecord> pair, SAMRecord read) {
		if(read.getSecondOfPairFlag()) {pair.setValue2(read);}
		if(read.getFirstOfPairFlag()) {pair.setValue1(read);}
		return pair;
	}
	
	
	public static Pair<String> assign(Pair<SAMRecord> pair, Map<String, IntervalTree<Annotation>> genes) {
		Collection<Annotation> overlappingGenes1=getOverlappers(pair.getValue1(), genes);
		Collection<Annotation> overlappingGenes2=getOverlappers(pair.getValue2(), genes);
		
		
		Collection<String> exonNames=new TreeSet<String>();
		Collection<String> intronNames=new TreeSet<String>();
		
		//Only consider annotations that both overlap
		Collection<Annotation> overlappingGenes=intersect(overlappingGenes1, overlappingGenes2);
		if(overlappingGenes.isEmpty()) {
			
			if(!overlappingGenes1.isEmpty() || !overlappingGenes2.isEmpty()) {
				return new Pair<String>(ambiguous, "");
			}
			else {
				return new Pair<String>(intergenic, "");
			}
		}
		
		SAMFragment read1=new SAMFragment(pair.getValue1());
		SAMFragment read2=new SAMFragment(pair.getValue2());
		
		Collection<String> states=new TreeSet<String>();
		for(Annotation gene: overlappingGenes) {
			
			String state1=assign(read1, gene);
			String state2=assign(read2, gene);
			
			if(state1.equals(exon) && state2.equals(exon)) {
				states.add(exon);
				exonNames.add(gene.getName());
			}
			if(state1.equals(intron) || state2.equals(intron)) {
				states.add(intron);
				intronNames.add(gene.getName());
			}
			
		}
		
		String state=consensusState(states);
		String geneNames=getConsensusGenes(state, exonNames, intronNames);
		
		Pair<String> rtrn=new Pair<String>();
		rtrn.setValue1(state);
		rtrn.setValue2(geneNames);
		return rtrn;
	}
	
	
	private static String getConsensusGenes(String state, Collection<String> exonNames, Collection<String> intronNames) {
		String rtrn="";
		
		if(state.equals(exon)) {rtrn=toString(exonNames);}
		else if(state.equals(intron)) {rtrn=toString(intronNames);}
		
		return rtrn;
	}

	private static String toString(Collection<String> exonNames) {
		String rtrn="";
		
		for(String n: exonNames) {
			rtrn+=n+",";
		}
		
		rtrn=rtrn.substring(0, rtrn.length()-1);
		
		return rtrn;
	}

	private static String assign(SAMFragment read, Annotation gene) {
		//if both overlap exons --> exon
		if(overlapsExons(gene, read)) {return exon;}
		
		else if(read.isSpliced()) {return ambiguous;}
		
		//if both overlap introns --> intron
		else if(overlapsIntrons(gene, read)){return intron;}
		
		
		else if(overlapsExonIntronJunction(gene, read)) {
			if(withinNofExon(gene, read, 5)) {return ambiguous;}
			return intron;
		}
		
		
		else if(overlapsIntronAtAll(gene, read)) {return intron;}
		
		//if one overlaps exon and one intron --> intron
		else if(overlapsExonAtAll(gene, read)) {return exon;}
		
		return ambiguous;
	}
	
	private static boolean withinNofExon(Annotation gene, SAMFragment read, int i) {
		SingleInterval t=new SingleInterval(read.getMateReferenceName(), read.getReferenceStartPosition(), read.getReferenceEndPosition()-i, read.getOrientation());
		SingleInterval t2=new SingleInterval(read.getMateReferenceName(), read.getReferenceStartPosition()+i, read.getReferenceEndPosition(), read.getOrientation());
		
		if(gene.fullyContained(t)) {return true;}
		if(gene.fullyContained(t2)) {return true;}
		
		return false;
	}

	private String assign(SAMRecord read, Map<String, IntervalTree<Annotation>> genes) {
		Collection<Annotation> overlappingGenes=getOverlappers(read, genes);
		
		if(overlappingGenes.isEmpty()) {return intergenic;}
		
		SAMFragment read1=new SAMFragment(read);
		
		Collection<String> states=new TreeSet<String>();
		for(Annotation gene: overlappingGenes) {
			String s=assign(read1, gene);
			states.add(s);
		}
		
		String state=consensusState(states);
		return state;
	}

	private static boolean overlapsExonAtAll(Annotation gene, SAMFragment read2) {
		return gene.overlaps(read2);
	}

	private static Collection<Annotation> intersect(Collection<Annotation> overlappingGenes1, Collection<Annotation> overlappingGenes2) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		for(Annotation gene: overlappingGenes1) {
			if(overlappingGenes2.contains(gene)) {rtrn.add(gene);}
		}
		
		return rtrn;
	}

	/*private String assign(SAMRecord read, Map<String, IntervalTree<Annotation>> genes) {
		Collection<Annotation> overlappingGenes=getOverlappers(read, genes);
		if(overlappingGenes.isEmpty()) {return intergenic;}
		
		SAMFragment f=new SAMFragment(read);
		
		Collection<String> states=new TreeSet<String>();
		
		for(Annotation gene: overlappingGenes) {
			//overlaps exons --> exon
			if(overlapsExons(gene, f)) {states.add(exon);}
			
			//overlaps intron --> intron
			else if(overlapsIntrons(gene, f)) {states.add(intron);}
			
			//overlaps intron+exon --> intron
			else if(overlapsExonIntronJunction(gene, f)) {states.add(intron);}
		}
		
		String state=consensusState(states);
		return state;
	}*/
	
	private static boolean overlapsExonIntronJunction(Annotation gene, SAMFragment f) {
		Annotation genomic=gene.getSingleInterval();
		return genomic.fullyContained(f);
	}

	private static boolean overlapsIntrons(Annotation gene, SAMFragment f) {
		Collection<Annotation> introns=gene.getIntrons();
		for(Annotation intron: introns) {
			if(intron.fullyContained(f)) {return true;}
		}
		return false;
	}
	
	private static boolean overlapsIntronAtAll(Annotation gene, SAMFragment f) {
		Collection<Annotation> introns=gene.getIntrons();
		for(Annotation intron: introns) {
			if(intron.overlaps(f)) {return true;}
		}
		return false;
	}

	
	
	private static String consensusState(Collection<String> states) {
		if(states.size()==1) {return states.iterator().next();}
		
		if(states.contains(exon)) {return exon;}
		
		return ambiguous;
	}

	private static boolean overlapsExons(Annotation g, SAMFragment f) {
		return g.fullyContained(f);
	}
	
	
	
	private static Collection<Annotation> getOverlappers(SAMRecord read, Map<String, IntervalTree<Annotation>> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		if(genes.containsKey(read.getReferenceName())) {
			IntervalTree<? extends Annotation> tree=genes.get(read.getReferenceName());
			Iterator<? extends Annotation> iter=tree.overlappingValueIterator(read.getAlignmentStart(), read.getAlignmentEnd());
			while(iter.hasNext()) {
				Annotation g=iter.next();
				boolean overlapsStrand=overlapsStrand(g, read);
				if(overlapsStrand) {rtrn.add(g);}
			}
		}
		return rtrn;
	}
	
	
	private static boolean overlapsStrand(Annotation g, SAMRecord read) {
		Strand readStrand=getStrand(read);
		if(g.getOrientation().equals(readStrand)) {return true;}
		return false;
	}
	
	private static Strand getStrand(SAMRecord read) {
		Strand s=Strand.UNKNOWN;
		if(!read.getFirstOfPairFlag()) {
			if(read.getReadNegativeStrandFlag()) {s=Strand.NEGATIVE;}
			else{s=Strand.POSITIVE;}
		}
		else {
			if(read.getReadNegativeStrandFlag()) {s=Strand.POSITIVE;}
			else{s=Strand.NEGATIVE;}
		}
		return s;
	}
	

	private void write(SAMRecord read, Pair<String> statePair, SAMFileWriter exonWriter, SAMFileWriter intronWriter, SAMFileWriter ambiguousWriter, SAMFileWriter interWriter) {
		String state=statePair.getValue1();
		//read.setAttribute("XT", state);
		//read.setAttribute("XL", statePair.getValue2());
		if(state.equals(exon)) {exonWriter.addAlignment(read);}
		else if(state.equals(intron)) {intronWriter.addAlignment(read);}
		else if(state.equals(ambiguous)) {ambiguousWriter.addAlignment(read);}
		else if(state.equals(intergenic)) {interWriter.addAlignment(read);}
	}
	
	
	private void write(Pair<SAMRecord> reads, Pair<String> state, SAMFileWriter exonWriter, SAMFileWriter intronWriter, SAMFileWriter ambiguousWriter, SAMFileWriter interWriter) {
		write(reads.getValue1(), state, exonWriter, intronWriter, ambiguousWriter, interWriter);
		write(reads.getValue2(), state, exonWriter, intronWriter, ambiguousWriter, interWriter);
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
		File bam=new File(args[0]);
		Map<String, IntervalTree<Annotation>> genes=BEDFileIO.loadTreeAnnotation(args[1]);
		String save=args[2];
		new AssignReads(bam, genes, save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=gene file (BED format) \n args[2]=output name";
}

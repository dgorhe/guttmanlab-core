package guttmanlab.core.sharp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.barcodeidentification.PeakCalling;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class ScoreFullFragment {
	private int totalCount;
	
	public ScoreFullFragment(File bam, int binSize, Map<String, IntervalTree<Gene>> genes, String save) throws IOException {
		preprocessPairs(bam, save);
		//Collection<Annotation> binScores=scoreFragment(bam, binSize, genes);
		//BEDFileIO.writeBED(binScores, save);
	}
	
	public ScoreFullFragment(File bam, String save) throws IOException {
		preprocessPairs(bam, save);
		//Collection<Annotation> binScores=scoreFragment(bam, binSize, genes);
		//BEDFileIO.writeBED(binScores, save);
	}
	
	
	private void preprocessPairs(File bam, String save) throws IOException {
		Map<String, Pair<SingleInterval>> pairs=new TreeMap<String, Pair<SingleInterval>>();
		FileWriter writer=new FileWriter(save);
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getProperPairFlag()){
				SingleInterval r=new SingleInterval(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
				r.setOrientation(getStrand(read));
				
				Pair<SingleInterval> pair=new Pair<SingleInterval>();
				if(pairs.containsKey(read.getReadName())) {
					pair=pairs.remove(read.getReadName());
					pair.setValue2(r);
					write(pair, writer);
				}
				else {
					pair.setValue1(r);
					pairs.put(read.getReadName(), pair);
				}
				
			}
			totalCount++;
			if(totalCount%1000000 ==0) {System.err.println(totalCount);}
		}	
		
		reads.close();
		inputReader.close();
		writer.close();
	}


	private void write(Pair<SingleInterval> pair, FileWriter writer) throws IOException {
		if(pair.getValue1().getReferenceName().equals(pair.getValue2().getReferenceName())) {
			int min=Math.min(pair.getValue1().getReferenceStartPosition(), pair.getValue2().getReferenceStartPosition());
			int max=Math.max(pair.getValue1().getReferenceEndPosition(), pair.getValue2().getReferenceEndPosition());
			writer.write(pair.getValue1().getReferenceName()+"\t"+min+"\t"+max+"\t"+pair.getValue1().getOrientation()+"\n");
		}
	}


	private Collection<Annotation> scoreFragment(File bam, int binSize, Map<String, IntervalTree<Gene>> genes){
		Collection<Annotation> rtrn=new ArrayList<Annotation>();
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<SingleInterval, Double> positionCount=new TreeMap<SingleInterval, Double>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()==255 && read.getProperPairFlag() && read.getFirstOfPairFlag()){
				
				
				//TODO Min to max (is max correct)
				//TODO find overlapping annotations
				//TODO intersect
				
				//If read1 and read2 in exons --> extend just exons
				//else extend fully
				
				/*if(overlapsExons(read, genes)) {
					//extend just over gene
				}
				
				int start=Math.min(read.getAlignmentStart(), read.getMateAlignmentStart());
				int end=Math.max(read.getAlignmentEnd(), read.getMateAlignmentStart()+read.getReadLength());
				
				SingleInterval fragment=new SingleInterval(read.getReferenceName(), start, end);
				
				System.out.println(fragment.toBED());*/
				
				Collection<Annotation> fragments=overlapsExons(read, genes);
				if(fragments.size()>1) {
					rtrn.addAll(fragments);
					//System.err.println(read.getReadName()+" "+fragments.size());
				}
				else {rtrn.addAll(fragments);}
				
				/*Collection<SingleInterval> bins=bin(fragment, binSize);
				
				for(SingleInterval bin: bins) {
					double score=0;
					if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
					score+=(1.0/bins.size());
					positionCount.put(bin, score);
				}*/
				
			}
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
		this.totalCount=totalCount;
		reads.close();
		inputReader.close();
	
		return rtrn;
	}
	
	private Collection<Annotation> overlapsExons(SAMRecord read, Map<String, IntervalTree<Gene>> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		if(!genes.containsKey(read.getReferenceName())) {
			return rtrn;
		}
		IntervalTree<Gene> tree=genes.get(read.getReferenceName());
		Iterator<Gene> iter1=tree.overlappingValueIterator(read.getAlignmentStart(), read.getAlignmentEnd());
		Iterator<Gene> iter2=tree.overlappingValueIterator(read.getMateAlignmentStart(), read.getMateAlignmentStart()+read.getReadLength());
		
		Collection<Gene> genes1=new TreeSet<Gene>();
		Collection<Gene> genes2=new TreeSet<Gene>();
		
		while(iter1.hasNext()) {
			Gene g=iter1.next();
			boolean overlapsExon=overlapsExons(g, read);
			boolean overlapsStrand=overlapsStrand(g, read);
			if(overlapsExon && overlapsStrand) {genes1.add(g);}
		}
		
		int min=Math.min(read.getAlignmentStart(), read.getMateAlignmentStart());
		int max=Math.max(read.getAlignmentEnd(), read.getMateAlignmentStart()+read.getReadLength());
		
		
		SingleInterval read2=new SingleInterval(read.getReferenceName(), read.getMateAlignmentStart(), read.getMateAlignmentStart()+read.getReadLength());
		read2.setOrientation(getStrand(read));
		while(iter2.hasNext()) {
			Gene g=iter2.next();
			boolean overlapsExon=overlapsExons(g, read2);
			boolean overlapsStrand=overlapsStrand(g, read2);
			if(overlapsExon && overlapsStrand) {genes2.add(g);}
		}
		
		boolean hasExon=false;
		for(Gene g: genes1) {
			if(genes2.contains(g)) {
				Annotation r=g.trim(min, max);
				r.setName(read.getReadName());
				//System.err.println(g.toBED()+" "+min+" "+max);
				rtrn.add(r);
				hasExon=true;
			}
		}
		
		/*if(!hasExon) {
			SingleInterval fragment=new SingleInterval(read.getReferenceName(), min, max);
			fragment.setOrientation(getStrand(read));
			rtrn.add(fragment);
		}*/
		
		return rtrn;
	}
	
	private boolean overlapsStrand(Gene g, SingleInterval read) {
		return read.getOrientation().equals(g.getOrientation());
	}


	private boolean overlapsExons(Gene g, SingleInterval read) {
		//return read.overlaps(g);
		return g.fullyContained(read);
	}

	
	private boolean overlapsStrand(Gene g, SAMRecord read) {
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


	private boolean overlapsExons(Gene g, SAMRecord read) {
		SAMFragment f=new SAMFragment(read);
		//return f.overlaps(g);
		return g.fullyContained(f);
		
	}
	

	private Collection<SingleInterval> bin(SingleInterval fragment, int binSize) {
		Collection<SingleInterval> rtrn=fragment.getBins(binSize);
		//System.err.println(fragment.toUCSC()+" "+fragment.getLength()+" "+rtrn);
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException {
		if(args.length>1) {
			File bam=new File(args[0]);
			String save=args[1];
			new ScoreFullFragment(bam, save);
		}

		else {System.err.println(usage);}
	}
	
	
	static String usage=" args[0]=bam file \n args[1]=save";
	
}

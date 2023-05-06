package guttmanlab.core.splicing.nanopore;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SMITAnalysis {
	
	int minDistanceFromEnd=50;
	int binSize;

	public SMITAnalysis(File bam, Map<String, IntervalTree<Gene>> geneTree, String save, int windowSize) throws IOException {
		this.binSize=windowSize;
		
		//iterate through each read
		//find overlapping gene
		//intersect gene by start and end of read
		Map<Gene, Map<Integer, Pair<Double>>> spliceScores=getSplicingScores(bam, geneTree);
		write(spliceScores, save);
		
	}
	
	

	private void write(Map<Gene, Map<Integer, Pair<Double>>> spliceScores, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: spliceScores.keySet()) {
			Map<Integer, Pair<Double>> map=spliceScores.get(gene);
			//Map<Integer, Pair<Double>> binnedMap=bin(map, binSize);
			for(Integer distance: map.keySet()) {
				Pair<Double> vals=map.get(distance);
				double spliceRatio=1-(vals.getValue1()/vals.getValue2());
				double runningAvg=1-getRunningAvg(map, distance, binSize);
				writer.write(gene.getName()+"\t"+(distance)+"\t"+spliceRatio+"\t"+runningAvg+"\t"+vals.getValue1()+"\t"+vals.getValue2()+"\n");
			}
		}
		
		writer.close();
	}

	private double getRunningAvg(Map<Integer, Pair<Double>> spliceScores, Integer distance, int binSize2) {
		double num=0;
		double denom=0;
		for(int i=distance; i<=distance+binSize2; i++) {
			if(spliceScores.containsKey(i)) {
				Pair<Double> pair=spliceScores.get(i);
				num+=pair.getValue1();
				denom+=pair.getValue2();
			}
		}
		
		return num/denom;
	}

	/*private Map<Integer, Pair<Double>> bin(Map<Integer, Pair<Double>> map, int binSize2) {
		Map<Integer, Pair<Double>> rtrn=new TreeMap<Integer, Pair<Double>>();
		
		for(int distance: map.keySet()) {
			int binned=distance/binSize2;
			Pair<Double> pair=map.get(distance);
			if(rtrn.containsKey(binned)) {
				Pair<Double> val=rtrn.get(binned);
				pair.setValue1(pair.getValue1()+val.getValue1());
				pair.setValue2(val.getValue2()+pair.getValue2());
			}
			rtrn.put(binned, pair);
		}
		
		return rtrn;
	}*/

	/*private void printFilteredReads(File bam, Map<String, IntervalTree<Gene>> geneTree, String save) {
		SAMFileReader reader=new SAMFileReader(bam);
		reader.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileHeader header=reader.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(header, false, new File(save));
		
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(!record.getNotPrimaryAlignmentFlag()) {
				Collection<Gene> genes=overlapsGene(record, geneTree);
				if(!genes.isEmpty()) {
					writer1.addAlignment(record);
				}
			}
			
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		
		reads.close();
		reader.close();
		writer1.close();
	}*/
	
	private Map<Gene, Map<Integer, Pair<Double>>> getSplicingScores(File bam, Map<String, IntervalTree<Gene>> geneTree, Collection<SingleInterval> excludeList) {
		SAMFileReader reader=new SAMFileReader(bam);
		reader.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator reads=reader.iterator();
		
		Map<String, Integer> geneCounts=new TreeMap<String, Integer>();
		Map<SingleInterval, Integer> pol2Counts=new TreeMap<SingleInterval, Integer>();
		
		//Distance by sum and total
		Map<Gene, Map<Integer, Pair<Double>>> fullSpliceScore=new TreeMap<Gene, Map<Integer, Pair<Double>>>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(!record.getNotPrimaryAlignmentFlag()) {
				Map<Gene, Map<Integer, Integer>> spliceScore=scoreDistanceFrom3SS(record, geneTree);
				
				Collection<SingleInterval> pol2Pos=getPos(record, spliceScore.keySet());
				if(!overlaps(pol2Pos, excludeList)) {
					update(geneCounts, spliceScore);
					add(spliceScore, fullSpliceScore);
					updatePol2(spliceScore, record, pol2Counts);
				}
				
			}
			
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		
		//print(geneCounts);
		print(pol2Counts);
		
		reads.close();
		reader.close();
		
		return fullSpliceScore;
		
	}

	private boolean overlaps(Collection<SingleInterval> pol2Pos, Collection<SingleInterval> excludeList) {
		for(SingleInterval r: excludeList) {
			for(SingleInterval r2: pol2Pos) {
				if(r.overlaps(r2)) {
					System.err.println(r.toUCSC()+" "+r2.toUCSC());
					return true;
				}
			}
		}
		return false;
	}



	private Collection<SingleInterval> getPos(SAMRecord record, Set<Gene> keySet) {
		Collection<SingleInterval> positions=new TreeSet<SingleInterval>();
		for(Gene g: keySet) {
			SingleInterval pos=getPol2Position(g, record);
			positions.add(pos);
		}
		return positions;
	}



	private Map<Gene, Map<Integer, Pair<Double>>> getSplicingScores(File bam, Map<String, IntervalTree<Gene>> geneTree) {
		SAMFileReader reader=new SAMFileReader(bam);
		reader.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator reads=reader.iterator();
		
		Map<String, Integer> geneCounts=new TreeMap<String, Integer>();
		Map<SingleInterval, Integer> pol2Counts=new TreeMap<SingleInterval, Integer>();
		
		//Distance by sum and total
		Map<Gene, Map<Integer, Pair<Double>>> fullSpliceScore=new TreeMap<Gene, Map<Integer, Pair<Double>>>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(!record.getNotPrimaryAlignmentFlag()) {
				Map<Gene, Map<Integer, Integer>> spliceScore=scoreDistanceFrom3SS(record, geneTree);
				update(geneCounts, spliceScore);
				add(spliceScore, fullSpliceScore);
				updatePol2(spliceScore, record, pol2Counts);
				
			}
			
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		
		//print(geneCounts);
		print(pol2Counts);
		
		reads.close();
		reader.close();
		
		return fullSpliceScore;
		
	}
	
	private void print(Map<SingleInterval, Integer> pol2Counts) {
		for(SingleInterval r: pol2Counts.keySet()) {
			System.out.println(r.toBedgraph(pol2Counts.get(r)));
		}
		
	}

	private void updatePol2(Map<Gene, Map<Integer, Integer>> spliceScore, SAMRecord record, Map<SingleInterval, Integer> pol2Counts) {
		Collection<SingleInterval> positions=new TreeSet<SingleInterval>();
		for(Gene g: spliceScore.keySet()) {
			SingleInterval pos=getPol2Position(g, record);
			positions.add(pos);
		}
		
		for(SingleInterval pos: positions) {
			int count=0;
			if(pol2Counts.containsKey(pos)) {count=pol2Counts.get(pos);}
			count++;
			pol2Counts.put(pos, count);
		}
		
	}

	private SingleInterval getPol2Position(Gene g, SAMRecord record) {
		if(g.getOrientation().equals(Strand.POSITIVE)) {return new SingleInterval(record.getReferenceName(),record.getAlignmentEnd(), record.getAlignmentEnd()+1);}
		return new SingleInterval(record.getReferenceName(),record.getAlignmentStart(), record.getAlignmentStart()+1);
	}

	/*private void print(Map<String, Integer> geneCounts) {
		for(String gene: geneCounts.keySet()) {
			int count=geneCounts.get(gene);
			System.out.println(gene+"\t"+count);
		}
	}*/

	private void update(Map<String, Integer> geneCounts, Map<Gene, Map<Integer, Integer>> spliceScore) {
		for(Gene g: spliceScore.keySet()) {
			int count=0;
			if(geneCounts.containsKey(g.getName())) {
				count=geneCounts.get(g.getName());
			}
			count++;
			geneCounts.put(g.getName(), count);
		}
	}

	private void add(Map<Gene, Map<Integer, Integer>> spliceScore, Map<Gene, Map<Integer, Pair<Double>>> fullSpliceScore) {
		for(Gene g: spliceScore.keySet()) {
			if(!fullSpliceScore.containsKey(g)) {
				fullSpliceScore.put(g, new TreeMap<Integer, Pair<Double>>());
			}
			
			Map<Integer, Pair<Double>> score=fullSpliceScore.get(g);
			Map<Integer, Integer> temp=spliceScore.get(g);
			
			for(Integer distance: temp.keySet()) {
				if(!score.containsKey(distance)) {score.put(distance, new Pair<Double>(0.0, 0.0));}
				Pair<Double> vals=score.get(distance);
				int ss=temp.get(distance);
				vals.setValue1(vals.getValue1()+ss);
				vals.setValue2(vals.getValue2()+1);
			}
		}
		
	}

	private Collection<Gene> overlapsGene(SAMRecord record, Map<String, IntervalTree<Gene>> geneTree) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		if(geneTree.containsKey(record.getReferenceName())) {
			IntervalTree<Gene> tree=geneTree.get(record.getReferenceName());
			Iterator<Gene> iter=tree.overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
			while(iter.hasNext()) {
				Gene gene=iter.next();
				if(readIsContainedWithinGene(record, gene)) {
					Collection<SingleInterval> overlappingJunctions=getOverlappingJunctions(record, gene);
					if(!overlappingJunctions.isEmpty()) {
							rtrn.add(gene);
					}
				}
			}
		}
		return rtrn;
	}
	

	private Map<Gene, Map<Integer, Integer>> scoreDistanceFrom3SS(SAMRecord record, Map<String, IntervalTree<Gene>> geneTree) {
		Map<Gene, Map<Integer, Integer>> rtrn=new TreeMap<Gene, Map<Integer, Integer>>();
		if(geneTree.containsKey(record.getReferenceName())) {
			IntervalTree<Gene> tree=geneTree.get(record.getReferenceName());
			Iterator<Gene> iter=tree.overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
			while(iter.hasNext()) {
				Gene gene=iter.next();
				if(readIsContainedWithinGene(record, gene)) {
					Collection<SingleInterval> overlappingJunctions=getOverlappingJunctions(record, gene);
					if(!overlappingJunctions.isEmpty()) {
						SAMFragment frag=new SAMFragment(record);
						Annotation a=frag.getAnnotation();
						
						
						//for each overlapping junction, determine if intron is contained within read or not
						Map<SingleInterval, Integer> spliced=isJunctionSpliced(record, overlappingJunctions);
						
						Map<Integer, Integer> map=getDistanceToSplicingState(spliced, gene, record);
						
						/*String name="";
						for(Integer distance: map.keySet()) {
							name+="d="+distance+" s="+map.get(distance)+",";
						}
						a.setName(name);
						System.out.println(a.toBED());*/
						
						rtrn.put(gene, map);
					}
				}
				
			}
		}
		return rtrn;
	}

	private Map<Integer, Integer> getDistanceToSplicingState(Map<SingleInterval, Integer> spliced, Gene gene, SAMRecord record) {
		Map<Integer, Integer> rtrn=new TreeMap<Integer, Integer>();
		int pol2Position=getPol2Position(record, gene);
		
		for(SingleInterval intron: spliced.keySet()) {
			int distance=getDistanceTo3SS(intron, pol2Position);
			int spliceScore=spliced.get(intron);
			rtrn.put(distance, spliceScore);
		}
		
		return rtrn;
	}

	private int getDistanceTo3SS(SingleInterval intron, int pol2Position) {
		if(intron.getOrientation().equals(Strand.POSITIVE)) {return pol2Position-intron.getReferenceEndPosition();}
		return intron.getReferenceStartPosition()-pol2Position;
	}

	private boolean hasUnspliced(Map<SingleInterval, Integer> spliced) {
		for(SingleInterval intron: spliced.keySet()) {
			int score=spliced.get(intron);
			if(score>0) {return true;}
		}
		return false;
	}

	private Map<SingleInterval, Integer> isJunctionSpliced(SAMRecord record, Collection<SingleInterval> overlappingJunctions) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		SAMFragment frag=new SAMFragment(record);
		for(SingleInterval intron: overlappingJunctions) {
			int score=0;
			SingleInterval intronAS=SingleInterval.antisense(intron);
			if(frag.overlaps(intron) || frag.overlaps(intronAS)) {score=1;}
			rtrn.put(intron, score);
		}
		
		return rtrn;
	}

	private int getPol2Position(SAMRecord record, Gene gene) {
		if(gene.getOrientation().equals(Strand.POSITIVE)) {return record.getAlignmentEnd();}
		return record.getAlignmentStart();
	}

	private Collection<SingleInterval> getOverlappingJunctions(SAMRecord record, Gene gene) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		Collection<SingleInterval> introns=gene.getIntronSet();
		//Iterator<SingleInterval> iter=gene.getBlocks();
		//while(iter.hasNext()) {
		for(SingleInterval intron: introns) {
			//SingleInterval exon=iter.next();
			int junctionPosition=intron.getReferenceEndPosition();
			if(gene.getOrientation().equals(Strand.NEGATIVE)) {junctionPosition=intron.getReferenceStartPosition();}
			if(record.getAlignmentStart()<junctionPosition && record.getAlignmentEnd()>junctionPosition) {rtrn.add(intron);}
		}
		
		return rtrn;
	}

	/*private boolean readSpansJunction(SAMRecord record, Gene gene) {
		Collection<Integer> junctionPositions=getJunctions(gene);
		
		for(Integer junctionPosition: junctionPositions) {
			if(record.getAlignmentStart()<junctionPosition && record.getAlignmentEnd()>junctionPosition) {return true;}
		}
		
		return false;
	}*/

	private Collection<Integer> getJunctions(Gene gene) {
		Collection<Integer> rtrn=new TreeSet<Integer>();
		Iterator<SingleInterval> iter= gene.getBlocks();
		while(iter.hasNext()) {
			SingleInterval exon=iter.next();
			if(gene.getOrientation().equals(Strand.POSITIVE)) {rtrn.add(exon.getReferenceStartPosition());}
			if(gene.getOrientation().equals(Strand.NEGATIVE)) {rtrn.add(exon.getReferenceEndPosition());}
		}
		return rtrn;
	}

	private boolean readIsContainedWithinGene(SAMRecord record, Gene gene) {
		int geneStart=gene.getReferenceStartPosition();
		int geneEnd=gene.getReferenceEndPosition();
		
		if(gene.getOrientation().equals(Strand.POSITIVE)) {geneEnd=geneEnd-minDistanceFromEnd;}
		if(gene.getOrientation().equals(Strand.NEGATIVE)) {geneStart=geneStart+minDistanceFromEnd;}
		
		if(record.getAlignmentStart()>=geneStart && record.getAlignmentEnd()<=geneEnd) {return true;}
		return false;
	}

	public static void main(String[] args) throws IOException {
		File file=new File(args[0]);
		Map<String, IntervalTree<Gene>> genes=BEDFileIO.loadTree(args[1]);
		String save=args[2];
		int windowSize=Integer.parseInt(args[3]);
		new SMITAnalysis(file, genes, save, windowSize);
	}
	
}

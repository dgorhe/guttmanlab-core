package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.ScanStat;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;


public class ChIPUtils {

	static double alpha=0.01;
	static double threshold=5.0;
	static double minCount=10;
	
	
	public static void scoresInRegion(File bam, Map<String, IntervalTree<SingleInterval>> regions, String save) throws IOException {
		
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Double> scores=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			Iterator<SingleInterval> iter=get(record, regions);
			while(iter.hasNext()) {
				SingleInterval r=iter.next();
				double score=0;
				if(scores.containsKey(r)) {score=scores.get(r);}
				score++;
				scores.put(r, score);
			}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		
		reader.close();
		reads.close();
		BEDFileIO.writeUCSCScore(scores, save);
	}
	
	
	private static Iterator<SingleInterval> get(SAMRecord record, Map<String, IntervalTree<SingleInterval>> tree) {
		if(tree.containsKey(record.getReferenceName())) {
			return tree.get(record.getReferenceName()).overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
		}
		return new ArrayList<SingleInterval>().iterator();
	}


	public static void getReadsPerWindow(File bam, int resolution, String save, Collection<SingleInterval> exclude) throws IOException {
		
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Double> scores=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment f=new SAMFragment(record);
			SingleInterval binned=f.getSingleInterval().bin(resolution);
			
			if(!overlaps(binned, exclude)) {
				double score=0;
				if(scores.containsKey(binned)) {
					score=scores.get(binned);
				}
				score++;
				scores.put(binned, score);
			}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		double median=median(scores, "chrX");
		scores=divideBy(scores, median);
		
		reader.close();
		reads.close();
		BEDFileIO.writeBEDGraph(scores, save);
	}
	
	
public static void getReadsPerWindow(File bam, int resolution, String save) throws IOException {
		
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Double> scores=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment f=new SAMFragment(record);
			SingleInterval binned=f.getSingleInterval().bin(resolution);
			
			double score=0;
			if(scores.containsKey(binned)) {
				score=scores.get(binned);
			}
			score++;
			scores.put(binned, score);
			
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		//double median=median(scores, "chrX");
		scores=divideBy(scores, counter);
		
		reader.close();
		reads.close();
		BEDFileIO.writeBEDGraph(scores, save);
	}


public static Map<SingleInterval, Double> getReadsPerWindow(File bam, int resolution) throws IOException {
	
	SAMFileReader reader=new SAMFileReader(bam);
	
	SAMRecordIterator reads=reader.iterator();
	
	Map<SingleInterval, Double> scores=new TreeMap<SingleInterval, Double>();
	
	int counter=0;
	while(reads.hasNext()){
		SAMRecord record=reads.next();
		
		SAMFragment f=new SAMFragment(record);
		SingleInterval binned=f.getSingleInterval().bin(resolution);
		
		double score=0;
		if(scores.containsKey(binned)) {
			score=scores.get(binned);
		}
		score++;
		scores.put(binned, score);
		
			
		counter++;
		if(counter%1000000 ==0){System.err.println(counter);}
	}
	
	//double median=median(scores, "chrX");
	scores=divideBy(scores, counter);
	
	reader.close();
	reads.close();
	return scores;
}

	public static Map<SingleInterval, Integer> getWindowCounts(File bam, int resolution) throws IOException {
		
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment f=new SAMFragment(record);
			SingleInterval binned=f.getSingleInterval().bin(resolution);
			//binned.setOrientation(f.getOrientation());
			
			int score=0;
			if(scores.containsKey(binned)) {
				score=scores.get(binned);
			}
			score++;
			scores.put(binned, score);
			
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		
		
		reader.close();
		reads.close();
		return scores;
	}
	
	
	public static Map<SingleInterval, Integer> getAllWindowCounts(File bam, int resolution) throws IOException {
		
		
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		
		Map<String, IntervalTree<SingleInterval>> tree=makeAllBins(CoordinateSpace.getGenomeLengths(reader.getFileHeader()), resolution);
		
		
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			updateAll(record, scores, tree);
			
			/*SAMFragment f=new SAMFragment(record);
			SingleInterval binned=f.getSingleInterval().bin(resolution);
			binned.setOrientation(f.getOrientation());
			
			int score=0;
			if(scores.containsKey(binned)) {
				score=scores.get(binned);
			}
			score++;
			scores.put(binned, score);*/
			
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		
		
		reader.close();
		reads.close();
		return scores;
	}
	
	
	private static Map<String, IntervalTree<SingleInterval>> makeAllBins(Map<String, Integer> genomeLengths, int resolution) {
		Map<String, IntervalTree<SingleInterval>> rtrn=new TreeMap<String, IntervalTree<SingleInterval>>();
		for(String chr: genomeLengths.keySet()) {
			int start=0;
			int end=genomeLengths.get(chr);
			IntervalTree<SingleInterval> tree=makeBins(chr, start, end, resolution);
			rtrn.put(chr, tree);
		}
		return rtrn;
	}


	private static IntervalTree<SingleInterval> makeBins(String chr, int start, int end, int resolution) {
		IntervalTree<SingleInterval> rtrn=new IntervalTree<SingleInterval>();
		
		for(int i=start; i<end; i++) {
			SingleInterval r=new SingleInterval(chr, i, i+resolution);
			rtrn.put(i,  i+resolution, r);
		}
		
		return rtrn;
	}


	private static void updateAll(SAMRecord record, Map<SingleInterval, Integer> scores, Map<String, IntervalTree<SingleInterval>> tree) {
		String chr=record.getReferenceName();
		int start=record.getAlignmentStart();
		int end=record.getAlignmentEnd();
		
		if(tree.containsKey(chr)) {
			Iterator<SingleInterval> iter=tree.get(chr).overlappingValueIterator(start, end);
			while(iter.hasNext()) {
				SingleInterval region=iter.next();
				update(region, scores);
			}
			
		}
		
	}


	private static void update(SingleInterval region, Map<SingleInterval, Integer> scores) {
		int score=0;
		if(scores.containsKey(region)) {score=scores.get(region);}
		score++;
		scores.put(region, score);
	}


	public static Collection<SingleInterval> getRegions(Map<String, IntervalTree<SingleInterval>> tree, SAMRecord record) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		if(tree.containsKey(record.getReferenceName())) {
			Iterator<SingleInterval> iter=tree.get(record.getReferenceName()).overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
			while(iter.hasNext()) {
				rtrn.add(iter.next());
			}
		}
		
		return rtrn;
	}
	
	
	public static Collection<SingleInterval> getRegions(Map<String, IntervalTree<SingleInterval>> tree, SingleInterval record) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		if(tree.containsKey(record.getReferenceName())) {
			Iterator<SingleInterval> iter=tree.get(record.getReferenceName()).overlappingValueIterator(record.getReferenceStartPosition(), record.getReferenceEndPosition());
			while(iter.hasNext()) {
				rtrn.add(iter.next());
			}
		}
		
		return rtrn;
	}
	
	public static Map<SingleInterval, Integer> score(File file, Collection<SingleInterval> regions) {
		Map<String, IntervalTree<SingleInterval>> tree=makeTree(regions);
		
		
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			Collection<SingleInterval> overlappers=getRegions(tree, record);
			for(SingleInterval o: overlappers) {
				increment(counts, o);
			}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		reader.close();
		reads.close();
		
		return counts;
	}
	
	private static void increment(Map<SingleInterval, Integer> counts, SingleInterval o) {
		int count=0;
		if(counts.containsKey(o)){
			count=counts.get(o);
		}
		count++;
		counts.put(o,  count);
	}
	
	public static Map<String, IntervalTree<SingleInterval>> makeTree(Collection<SingleInterval> regions) {
		Map<String, IntervalTree<SingleInterval>> tree=new TreeMap<String, IntervalTree<SingleInterval>>();
		
		for(SingleInterval r: regions) {
			if(!tree.containsKey(r.getReferenceName())) {
				tree.put(r.getReferenceName(), new IntervalTree<SingleInterval>());
			}
			tree.get(r.getReferenceName()).put(r.getReferenceStartPosition(), r.getReferenceEndPosition(), r);
		}
		
		
		return tree;
	}

	private static int totalCount(File bam) {
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		return counter;
	}

	
	
	public static void CountByChr(File bam, String save, Collection<SingleInterval> exclude) throws IOException {
		
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<String, Double> counts=new TreeMap<String, Double>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment f=new SAMFragment(record);
			SingleInterval binned=f.getSingleInterval().bin(10000);
			
			if(!overlaps(binned, exclude)) {
				double score=0;
				if(counts.containsKey(binned.getReferenceName())) {
					score=counts.get(binned.getReferenceName());
				}
				score++;
				counts.put(binned.getReferenceName(), score);
			}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		/*double median=median(scores);
		scores=divideBy(scores, median);*/
		
		reader.close();
		reads.close();
		BEDFileIO.write(counts, save);
	}
	
	private static boolean overlaps(SingleInterval binned, Collection<SingleInterval> exclude) {
		for(SingleInterval e: exclude) {
			if(binned.overlaps(e)){return true;}
		}
		return false;
	}


	private static double median(Map<SingleInterval, Double> scores) {
		List<Double> list=new ArrayList<Double>();
		list.addAll(scores.values());
		return Statistics.quantile(list, 0.5);
	}
	
	public static Map<SingleInterval, Double> medianNormalize(Map<SingleInterval, Double> scores) {
		double median=median(scores);
		
		Map<SingleInterval, Double> rtrn=divideBy(scores, median);
		
		return rtrn;
	}
	
	
	private static double average(Map<SingleInterval, Double> scores) {
		List<Double> list=new ArrayList<Double>();
		list.addAll(scores.values());
		return Statistics.mean(list);
	}
	
	private static double median(Map<SingleInterval, Double> scores, String chr) {
		List<Double> list=new ArrayList<Double>();
		
		for(SingleInterval r: scores.keySet()) {
			if(r.getReferenceName().equalsIgnoreCase(chr)) {list.add(scores.get(r));}
		}
		
		//list.addAll(scores.values());
		return Statistics.quantile(list, 0.5);
	}


	private static Map<SingleInterval, Double> divideBy(Map<SingleInterval, Double> scores, double denom) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval r: scores.keySet()) {
			double score=scores.get(r);
			double norm=(score/denom);
			rtrn.put(r, norm);
		}
		
		return rtrn;
	}


	public static double countReadsInRegion(File bam, SingleInterval region) {
		
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.query(region.getReferenceName(), region.getReferenceStartPosition(), region.getReferenceEndPosition(), true);
		
		//SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		
		reader.close();
		reads.close();
		return counter;
	}
	
	
	private static double numberOfReads(File bam) {
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		return counter;
	}
	
	
	

	private static Iterator<SingleInterval> getOverlaps(SAMRecord record, Map<String, IntervalTree<SingleInterval>> windowTree) {
		if(windowTree.containsKey(record.getReferenceName())) {
			IntervalTree<SingleInterval> tree=windowTree.get(record.getReferenceName());
			return tree.overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
		}
		return new ArrayList<SingleInterval>().iterator();
	}
	
	
	private static void normalize(Map<SingleInterval, Double> vals1, Map<SingleInterval, Double> vals2, String save) throws IOException {
		Map<SingleInterval, Double> norm=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval r: vals1.keySet()) {
			double num=vals1.get(r);
			double denom=0;
			if(vals2.containsKey(r)) {
				denom=Math.max(vals2.get(r), denom);
			}
			norm.put(r, (num+1)/(denom+1));
		}
		
		BEDFileIO.writeBEDGraph(norm, save);
	}
	
	public static void significantWindows(File bam1, File bam2, int windowSize, String save) throws IOException {
		//FileWriter writer=new FileWriter(save+".bed");
		FileWriter writerBEDGraph=new FileWriter(save+".bedgraph");
		
		Map<SingleInterval, Integer> counts1=getWindowCounts(bam1, windowSize);
		Map<SingleInterval, Integer> counts2=getWindowCounts(bam2, windowSize);
		int sampleTotal=totalCount(bam1);
		int inputTotal=totalCount(bam2);
		
		double inputMedian=scoreMedian(counts2);
		double sampleMedian=scoreMedian(counts1);
		
		for(SingleInterval r: counts1.keySet()) {
			double binScore=counts1.get(r);
			double geneScore=sampleMedian;
			double enrichment1=binScore/geneScore;
	
			double window2=0;
			if(counts2.containsKey(r)) {
				window2=counts2.get(r);
			}
			
			window2=Math.max(window2, inputMedian);
			
			double inputWindowScore=window2;
			
				
			double enrichment2=getEnrichment(binScore, inputWindowScore, sampleTotal, inputTotal);
				
			double overallEnrichment=Math.min(enrichment1, enrichment2);
			
			writerBEDGraph.write(r.toBedgraph(overallEnrichment)+"\n");
		}
		
		writerBEDGraph.close();
		
	}
	
	private static int scoreMedian(Map<SingleInterval, Integer> binScores) {
		List<Integer> vals=new ArrayList<Integer>();
		vals.addAll(binScores.values());
		
		
		int median=Statistics.quantileInt(vals, 0.5);
		
		return median;
	}
	
	
	/*public static void significantWindows(File bam1, File bam2, int windowSize, String save, Map<SingleInterval, Collection<Gene>> binsToGenes, Map<Gene, Integer> geneScores) throws IOException {
		FileWriter writer=new FileWriter(save+".bed");
		FileWriter writerBEDGraph=new FileWriter(save+".bedgraph");
		
		Map<SingleInterval, Integer> counts1=getWindowCounts(bam1, windowSize);
		Map<SingleInterval, Integer> counts2=getWindowCounts(bam2, windowSize);
		int total1=totalCount(bam1);
		int total2=totalCount(bam2);
		
		for(SingleInterval r: counts1.keySet()) {
			int window1=counts1.get(r);
			int window2=0;
			if(counts2.containsKey(r)) {
				window2=counts2.get(r);
			}
			
			window2=getGeneScores(r, binsToGenes, geneScores, window2);
			
			double enrichment=getEnrichment(window1, window2, total1, total2);
			writerBEDGraph.write(r.toBedgraph(enrichment)+"\n");
			double p=getPValue(window1, window2, total1, total2);
			if(p<0.001 && window1+window2>minCount && enrichment>threshold) {
				writer.write(r.toBedgraph(p)+"\n");
			}
		}
		writer.close();
		writerBEDGraph.close();
		
	}*/
	
	
	
	public static void significantWindows(File bam1, File bam2, int windowSize, String save, Map<String, IntervalTree<Gene>> genes, Map<Gene, Integer> geneScores, boolean rewrite) throws IOException {
		FileWriter writer=new FileWriter(save+".bed");
		FileWriter writerBEDGraph=new FileWriter(save+".bedgraph");
		
		Map<SingleInterval, Integer> counts1=getWindowCounts(bam1, windowSize);
		Map<SingleInterval, Integer> counts2=getWindowCounts(bam2, windowSize);
		int total1=totalCount(bam1);
		int total2=totalCount(bam2);
		
		for(SingleInterval r: counts1.keySet()) {
			int window1=counts1.get(r);
			int window2=0;
			if(counts2.containsKey(r)) {
				window2=counts2.get(r);
			}
			
			window2=getGeneScores(r, genes, geneScores, window2);
			//System.out.println(r.toBedgraph(window2));
			
			double enrichment=getEnrichment(window1, window2, total1, total2);
			writerBEDGraph.write(r.toBedgraph(enrichment)+"\n");
			double p=getPValue(window1, window2, total1, total2);
			if(p<0.001 && window1+window2>minCount && enrichment>threshold) {
				writer.write(r.toShortBED()+"\n");
			}
		}
		writer.close();
		writerBEDGraph.close();
		
	}
	
	
	public static void significantWindowsWithControls(File bam1, File bam2, Collection<SingleInterval> regions, String save) throws IOException {
		FileWriter writer=new FileWriter(save+".bed");
		
		Map<SingleInterval, Integer> counts1=score(bam1, regions);
		Map<SingleInterval, Integer> counts2=score(bam2, regions);
		int total1=totalCount(bam1);
		int total2=totalCount(bam2);
		double totalLength=genomeLength(bam1);
		double lambda=total1/totalLength;
		
		for(SingleInterval r: counts1.keySet()) {
			int window1=counts1.get(r);
			int window2=0;
			if(counts2.containsKey(r)) {
				window2=counts2.get(r);
			}
			
			
			double scanP=ScanStat.getPValue(window1, lambda, r.size(), totalLength);
			double p=getPValue(window1, window2, total1, total2);
			double enrich1=getEnrichment(window1, window2, total1, total2);
			double enrich2=((double)window1/(double)r.size())/lambda;
			if(p<alpha && scanP<alpha) {
				writer.write(r.toShortBED("E1="+enrich1+" E2="+enrich2)+"\n");
			}
		}
		writer.close();
	}
	
	
	
	
	/*private static int getGeneScores(SingleInterval r, Map<SingleInterval, Collection<Gene>> binsToGenes, Map<Gene, Integer> geneScores, int window2) {
		int rtrn=window2;
		
		if(binsToGenes.containsKey(r)) {
			Collection<Gene> genes=binsToGenes.get(r);
			
			for(Gene g: genes) {
				if(geneScores.containsKey(g)) {
					double score=(double)r.size()* ((double)geneScores.get(g)/(double)g.size());
					int score2=new Double(score).intValue()+1;
					System.out.println(r.toUCSC()+"\t"+ g.toUCSC()+" "+ window2+" "+score+" "+score2);
					rtrn=Math.max(score2,  rtrn);
				}
			}
		}
		
		return rtrn;
		
	}*/
	
	private static double genomeLength(File bam1) {
		SAMFileReader reader=new SAMFileReader(bam1);
		double rtrn=CoordinateSpace.getGenomeLength(reader.getFileHeader());
		reader.close();
		return rtrn;
	}


	private static int getGeneScores(SingleInterval r, Map<String, IntervalTree<Gene>> geneTree, Map<Gene, Integer> geneScores, int window2) {
		int rtrn=window2;
		
		Collection<Gene> genes=getOverlappingGenes(r, geneTree);
		
		for(Gene g: genes) {
			if(geneScores.containsKey(g)) {
				double score=(double)r.size()* ((double)geneScores.get(g)/(double)g.size());
				int score2=new Double(score).intValue()+1;
				//System.out.println(r.toUCSC()+"\t"+ g.toUCSC()+" "+ window2+" "+score+" "+score2);
				rtrn=Math.max(score2,  rtrn);
			}
		}
		
		
		return rtrn;
		
	}


	private static Collection<Gene> getOverlappingGenes(SingleInterval read, Map<String, IntervalTree<Gene>> geneTree) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		if(geneTree.containsKey(read.getReferenceName())) {
			IntervalTree<Gene> tree=geneTree.get(read.getReferenceName());
			Iterator<Gene> iter=tree.overlappingValueIterator(read.getReferenceStartPosition(), read.getReferenceEndPosition());
			while(iter.hasNext()) {
				Gene g=iter.next();
				//rtrn.add(g);
				boolean overlapsExon=overlapsExons(g, read);
				boolean overlapsStrand=overlapsStrand(g, read);
				if(overlapsExon && overlapsStrand) {rtrn.add(g);}
			}
		}
		return rtrn;
	}
		
	private static boolean overlapsStrand(Gene g, SingleInterval read) {
		return read.getOrientation().equals(g.getOrientation());
	}


	private static boolean overlapsExons(Gene g, SingleInterval read) {
		//return read.overlaps(g);
		return g.fullyContained(read);
	}


	private static double getEnrichment(double window1, double window2, int total1, int total2) {
		double num=((double)window1+1)/(double)total1;
		double denom=((double)window2+1)/(double)total2;
		return num/denom;
	}


	private static double getPValue(int window_sampleCount, int window_inputCount, int total_sampleCount, int total_inputCount) {
		//Return binomial p
		int n=window_inputCount+window_sampleCount;
		if(n>0){
			double p=(double)total_sampleCount/(double)(total_inputCount+total_sampleCount);
			Binomial b=new Binomial(n, p);
			return 1-b.cdf(window_sampleCount);
		}
		return 1.0;
	}

	
	private static void normalizeBedgraph(Map<SingleInterval, Double> vals1, String save) throws IOException {
		double median=median(vals1, "chrX");
		
		
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: vals1.keySet()) {
			double score=vals1.get(r);
			double norm=score/median;
			writer.write(r.toBedgraph(norm)+"\n");
		}
		
		writer.close();
	}
	
	private static Collection<SingleInterval> collapse(Collection<SingleInterval> sigRegions) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		Iterator<SingleInterval> iter=sigRegions.iterator();
		
		ArrayList<SingleInterval> list=new ArrayList<SingleInterval>();
		SingleInterval current=null;
		while(iter.hasNext()) {
			SingleInterval next=iter.next();
			if(overlaps(current, next)) {
				list.add(next);
				current=next;
			}
			else {
				SingleInterval c=collapseList(list);
				rtrn.add(c);
				list=new ArrayList<SingleInterval>();
				list.add(next);
				current=next;
			}
		}
		
		return rtrn;
	}
	
	private static boolean overlaps(SingleInterval current, SingleInterval next) {
		if(current==null) {return true;}
		
		if(next.getReferenceStartPosition()==(current.getReferenceEndPosition())) {
			if(next.getReferenceName().equals(current.getReferenceName())) {return true;}
		}
		
		return false;
	}

	private static SingleInterval collapseList(ArrayList<SingleInterval> list) {
		SingleInterval first=list.get(0);
		
		SingleInterval last=list.get(list.size()-1);
		
		return new SingleInterval(first.getReferenceName(), first.getReferenceStartPosition(), last.getReferenceEndPosition());
	}

	private static void write(String save, Collection<SingleInterval> collapsed) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: collapsed) {writer.write(r.toShortBED()+"\n");}
		
		writer.close();
	}
	
	
	//Take the merged BAM file and randomly sample a new BAM of 
	private static void enrichment(File bam1, File bam2, int binSize) throws IOException {
		Map<SingleInterval, Integer> counts1=getWindowCounts(bam1, binSize);
		Map<SingleInterval, Integer> counts2=getWindowCounts(bam2, binSize);
		//double total1=totalCount(bam1);
		//double total2=totalCount(bam2);
		
		Collection<SingleInterval> allRegions=new TreeSet<SingleInterval>();
		allRegions.addAll(counts1.keySet());
		allRegions.addAll(counts2.keySet());
		
		for(SingleInterval r: allRegions) {
			double window1=0;
			if(counts1.containsKey(r)) {
				window1=counts1.get(r);
			}
			double window2=0;
			if(counts2.containsKey(r)) {
				window2=counts2.get(r);
			}
			
			//double enrich1=getEnrichment(window1, window2, total1, total2);
			double enrich1=(window1)-(window2);
			System.out.println(r.toUCSC()+"\t"+enrich1);	
		}
	}
	
	
	private static void smooth(File bam1, int binSize, String save) throws IOException {
		Map<SingleInterval, Integer> counts1=getWindowCounts(bam1, binSize);
		//extend up and down by adding and subtracting start and end positions from counts
		writeSmoothed(counts1, save);
	}
	

	private static void writeSmoothed(Map<SingleInterval, Integer> counts1, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: counts1.keySet()) {
			writer.write(r.getMidPoint().toBedgraph(counts1.get(r))+"\n");
		}
		
		writer.close();
	}

	
	private static void windowCounts(File bam1, int binSize, String save) throws IOException {
		Map<SingleInterval, Integer> counts1=getWindowCounts(bam1, binSize);
		
		write(save, counts1);
		
	}

	private static void write(String save, Map<SingleInterval, Integer> counts1) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: counts1.keySet()) {
			writer.write(r.toUCSC()+"\t"+counts1.get(r)+"\n");
		}
		
		writer.close();
	}

	
	private static void mergeFiles(File[] files, String save) throws IOException {
		Map<String, Integer>[] maps=new Map[files.length];
		Collection<String> allStrings=new TreeSet<String>();
		FileWriter writer=new FileWriter(save);
		writer.write("position");
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i]);
			maps[i]=getScores(files[i]);
			allStrings.addAll(maps[i].keySet());
			writer.write("\t"+files[i].getName());
		}
		writer.write("\n");
		
		
		for(String key: allStrings) {
			writer.write(key);
			for(int i=0; i<maps.length; i++) {
				double score=get(key, maps[i]);
				writer.write("\t"+score);
			}
			writer.write("\n");
		}
		
		writer.close();
	}
	

	private static Map<String, Integer> getScores(File file) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		List<String> lines=BEDFileIO.loadLines(file.getAbsolutePath());
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			String key=tokens[0];
			int score=Integer.parseInt(tokens[1]);
			rtrn.put(key, score);
		}
		return rtrn;
	}


	private static double get(String key, Map<String, Integer> map) {
		if(map.containsKey(key)) {return map.get(key);}
		return 0;
	}


	public static void main (String[] args) throws IOException {
		if(args.length>1) {
			File bam1=new File(args[0]);
			String save=args[1];
			int binSize=Integer.parseInt(args[2]);
			windowCounts(bam1, binSize, save);
			
			/*File[] files=new File(args[0]).listFiles();
			String save=args[1];
			mergeFiles(files, save);*/
			
		}
		else {System.err.println(usage);}
		
	}

	

	static String usage=" args[0]=sample bam \n args[1]=save \n args[2]=binSize";

	

	

	
	
}

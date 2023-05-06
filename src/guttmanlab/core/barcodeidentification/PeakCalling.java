package guttmanlab.core.barcodeidentification;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.PoissonDistribution;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.ScanStat;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PeakCalling {

	//private Map<Annotation, Double> geneScores;
	private Map<Annotation, Integer> geneCounts;
	double alpha=0.01;
	double threshold=5.0;
	double minCount=10;
	int numPerm=10;
	int extension=10;
	private int totalCount;
	int binSize;
	int numberOfWindows;
	
	Map<Annotation, Integer> binCounts;
	Map<Annotation, Double> binPvals;
	Map<Annotation, Double> binEnrichment;
	Map<String, IntervalTree<Annotation>> geneTree;
	

	
	//public Map<Annotation, Double> getGeneScores(){return this.geneScores;}
	public int getTotalCount() {return totalCount;}
	
	public void setTotalCount(int total) {
		this.totalCount=total;
	}
	
	
	public PeakCalling(File bam, int binSize, Map<String, IntervalTree<Annotation>> genes) throws IOException {
		this.geneTree=genes;
		this.binSize=binSize;
		
		this.binCounts=score(bam, binSize, geneTree);
	}
	
	private void score() {
		this.binPvals=new TreeMap<Annotation, Double>();
		this.binEnrichment=new TreeMap<Annotation, Double>();
		
		for(Annotation bin: binCounts.keySet()) {
			Collection<Annotation> overlappingGenes=this.getGenes(bin, geneTree);
			int score=binCounts.get(bin);
			Annotation g=getMaxGene(geneCounts, overlappingGenes);
			boolean use=false;
			double expected=1;
			if(g!=null && geneCounts.containsKey(g)) {
				int geneScore=geneCounts.get(g);
				double lambda=(double)geneScore/(double)g.size();
				expected=lambda;
				use=true;
			}
			else if(isIntergenic(bin, geneTree)) {
				//System.out.println(bin.toBED());
				g=extend(bin, extension);
				double lambda=getScore(bin, extension, binCounts)/(double)binSize;
				expected=lambda;
				use=true;
			}
			if(use) {
				double observed=(double)score/(double)binSize;
				double enrichment=observed/expected;
				double pval=ScanStat.getPValue(score, expected, binSize, g.size());
				binEnrichment.put(bin, enrichment);
				binPvals.put(bin, pval);
			}
			//else {System.out.println(bin.toBED());}
		}
	}
	
	public Map<Annotation, Double> sampleGeneScores() {
		Map<Annotation, Double> rtrn=new TreeMap<Annotation, Double>();
		
		int counter=0;
		for(Annotation g: geneCounts.keySet()) {
			int total=geneCounts.get(g);
			double lambda=(double)total/(double)g.size();
			lambda=binSize*lambda;
			int numberBins=g.size()/binSize;
			double perm=getPerms(lambda, numberBins, numPerm);
			rtrn.put(g, perm);
			counter++;
			if(counter%1000==0) {System.err.println("sampling "+counter+" "+geneCounts.size());}
		}
		
		return rtrn;
	}
	
	private double getPerms(double lambda, int numberBins, int numPerm2) {
		double sum=0;
		PoissonDistribution dist=new PoissonDistribution(lambda);
		
		
		for(int i=0; i<numPerm2; i++) {
			sum+=max(dist, numberBins);
		}
		
		return sum/(double)numPerm2;
	}
	
	public double getPerms(double lambda) {
		double sum=0;
		PoissonDistribution dist=new PoissonDistribution(lambda);
		int numberBins=2*extension;
		
		for(int i=0; i<numPerm; i++) {
			sum+=max(dist, numberBins);
		}
		
		return sum/(double)numPerm;
	}
	
	
	private double max(PoissonDistribution dist, int numberBins) {
		
		int max=0;
		for(int i=0; i<numberBins; i++) {
			int val=dist.sample();
			max=Math.max(val, max);
		}
		return max;
	}
	
	private Annotation extend(Annotation bin, int extension2) {
		int newStart=bin.getReferenceStartPosition()+(-extension2*bin.size());
		int newEnd=bin.getReferenceStartPosition()+(extension2*bin.size())+bin.size();
		
		Annotation rtrn=new SingleInterval(bin.getReferenceName(), newStart, newEnd);
		rtrn.setOrientation(bin.getOrientation());
		return rtrn;
	}
	
	private double getScore(Annotation bin, int extension2, Map<Annotation, Integer> binCounts2) {
		int sum=0;
		int count=0;
		for(int start=-extension2; start<extension2; start++) {
			int newStart=bin.getReferenceStartPosition()+(start*bin.size());
			int newEnd=newStart+bin.size();
			Annotation r=new SingleInterval(bin.getReferenceName(), newStart, newEnd);
			r.setOrientation(bin.getOrientation());
			sum+=get(r, binCounts2);
			count++;
		}
		
		double expected= (double)sum/(double)count;
		//System.out.println(bin.getReferenceName()+"\t"+bin.getReferenceStartPosition()+(-extension2*bin.size())+"\t"+bin.getReferenceStartPosition()+(extension2*bin.size())+"\t"+expected);
		return expected;
		
	}
	
	
	public double getScore(Annotation bin) {
		int sum=0;
		int count=0;
		for(int start=-extension; start<extension; start++) {
			int newStart=bin.getReferenceStartPosition()+(start*bin.size());
			int newEnd=newStart+bin.size();
			Annotation r=new SingleInterval(bin.getReferenceName(), newStart, newEnd);
			r.setOrientation(bin.getOrientation());
			sum+=get(r, binCounts);
			count++;
		}
		
		
		int newStart=bin.getReferenceStartPosition()+(-extension*bin.size());
		int newEnd=bin.getReferenceStartPosition()+(extension*bin.size())+binSize;
		double expected= (double)sum/(double)count;
		//System.out.println(bin.getReferenceName()+"\t"+newStart+"\t"+newEnd+"\t"+sum);
		return expected;
		
	}
	
	private int get(Annotation r, Map<Annotation, Integer> binCounts2) {
		if(binCounts2.containsKey(r)) {return binCounts2.get(r);}
		return 0;
	}
	private void bin(Annotation gene, int binSize, IntervalTree<Annotation> tree) {
		Collection<Annotation> windows=gene.getSplicedWindows(binSize);
		for(Annotation bin: windows) {
			tree.put(bin.getReferenceStartPosition(), bin.getReferenceEndPosition(), bin);
		}
	}
	
	/*public int getMaxGeneScore(Annotation bin) {
		Annotation g=getMaxGene(bin);
		return getMaxGeneScore(bin, g);
		
		
	}*/
	
	/*public int getMaxGeneScore(Annotation bin, Annotation g) {
		int rtrn=0;
		if(g!=null && geneScores.containsKey(g)) {
			double score=geneScores.get(g);
			rtrn=(int)Math.ceil(score);	
		}
		else {
			double scaled=getPerms(getScore(bin, extension, this.binCounts), 2*extension, numPerm);
			rtrn=(int)Math.ceil(scaled);
		}
			
		return rtrn;
	}*/
	
	public Annotation getMaxGene(Annotation bin) {
		Collection<Annotation> overlappingGenes=this.getGenes(bin, geneTree);
		Annotation g=getMaxGene(geneCounts, overlappingGenes);
		return g;	
	}
	
	private Annotation getMaxGene(Map<Annotation, Integer> geneScores2, Collection<Annotation> overlappingGenes) {
		double max=0;
		Annotation maxGene=null;
		
		for(Annotation g: overlappingGenes) {
			if(geneScores2.containsKey(g)) {
				double score=(double)geneScores2.get(g)/(double)g.size();
				if(score>max) {
					max=score;
					maxGene=g;
				}
			}
			//else {System.err.println(g);}
		}
		
		return maxGene;
	}
	
	
	private int get(Map<Annotation, Integer> binCounts2, Annotation bin) {
		int score=0;
		if(binCounts2.containsKey(bin)) {score=binCounts2.get(bin);}
		return score;
	}
	
	
	public Map<Annotation, Double> getPvalues(){
		if(this.binPvals==null || this.binPvals.isEmpty()) {
			score();
			System.err.println("scoring");
		}
		
		return this.binPvals;
	}
	
	public Map<Annotation, Double> getEnrichments(){
		if(this.binEnrichment==null || this.binEnrichment.isEmpty()) {
			score();
			System.err.println("scoring");
		}
		return this.binEnrichment;
	}
	
	public Map<Annotation, Integer> getCounts(){return this.binCounts;}
	
	
	
	
	private TreeMap<Annotation, Integer> score(File bam, Map<String, IntervalTree<Annotation>> geneBins){
		//TODO For introns, we should estimate gene coverage over unique regions
		//TODO but we should score windows across entire genomic interval (including exons)
		this.geneCounts=new TreeMap<Annotation, Integer>();
		
		System.err.println("Made gene bins");
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<Annotation, Integer> positionCount=new TreeMap<Annotation, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				Collection<Annotation> allBins=getGenes(read, geneBins);
				
				for(Annotation bin: allBins) {
					int score=0;
					if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
					score++;
					positionCount.put(bin, score);
				}
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount+" "+read.getReferenceName());}
		}
				
		this.totalCount=totalCount;
		reads.close();
		inputReader.close();
		//this.numberOfWindows=positionCount.keySet().size();
		this.numberOfWindows=numWindows(geneBins);
		
		System.out.println("Number of windows over chromosomes:"+this.numberOfWindows+" total windows:"+positionCount.keySet().size());
		
		return positionCount;
	}

	private int numWindows(Map<String, IntervalTree<Annotation>> geneBins) {
		int sum=0;
		for(String chr: geneBins.keySet()) {
			IntervalTree<Annotation> tree=geneBins.get(chr);
			sum+=tree.size();
		}
		return sum;
	}

	private TreeMap<Annotation, Integer> score(File bam, int binSize, Map<String, IntervalTree<Annotation>> genes){
		//TODO For introns, we should estimate gene coverage over unique regions
		//TODO but we should score windows across entire genomic interval (including exons)
		
		
		
		this.geneCounts=new TreeMap<Annotation, Integer>();
		
		Map<String, IntervalTree<Annotation>> geneBins=getBins(genes, binSize);
		
		System.err.println("Made gene bins");
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<Annotation, Integer> positionCount=new TreeMap<Annotation, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				Collection<Annotation> overlappingGenes=getGenes(read, genes);
				Collection<Annotation> allBins=new TreeSet<Annotation>();
				if(overlappingGenes.size()>0) {
					Collection<Annotation> additionalBins=this.getGenes(read, geneBins);
					allBins.addAll(additionalBins);
				}
				allBins.addAll(SAMFragment.allBins(read, binSize));
				for(Annotation bin: allBins) {
					int score=0;
					if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
					score++;
					positionCount.put(bin, score);
					add(overlappingGenes);
				}
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount+" "+read.getReferenceName());}
		}
				
		this.totalCount=totalCount;
		reads.close();
		inputReader.close();
		//this.numberOfWindows=positionCount.keySet().size();
		this.numberOfWindows=numWindows(positionCount, genes.keySet());
		
		System.out.println("Number of windows over chromosomes:"+this.numberOfWindows+" total windows:"+positionCount.keySet().size());
		
		return positionCount;
	}
	
	
	
	
	
	private int numWindows(TreeMap<Annotation, Integer> positionCount, Set<String> keySet) {
		int count=0; 
		
		for(Annotation r: positionCount.keySet()) {
			if(keySet.contains(r.getReferenceName())) {count++;}
		}
		
		return count;
	}

	private Map<String, IntervalTree<Annotation>> getBins(Map<String, IntervalTree<Annotation>> genes, int binSize) {
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		
		
		for(String chr: genes.keySet()) {
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			Iterator<Annotation> iter=genes.get(chr).valueIterator();
			while(iter.hasNext()) {
				Annotation g=iter.next();
				if(!g.getName().contains("intron")) {
					bin(g, binSize, tree);
				}
			}
			rtrn.put(chr, tree);
		}
		
		return rtrn;
	}

	private Collection<Annotation> getGenes(Annotation read, Map<String, IntervalTree<Annotation>> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		if(genes.containsKey(read.getReferenceName())) {
			IntervalTree<? extends Annotation> tree=genes.get(read.getReferenceName());
			Iterator<? extends Annotation> iter=tree.overlappingValueIterator(read.getReferenceStartPosition(), read.getReferenceEndPosition());
			while(iter.hasNext()) {
				Annotation g=iter.next();
				boolean overlapsExon=overlapsExons(g, read);
				boolean overlapsStrand=overlapsStrand(g,read);
				//boolean overlapsIntron=overlapsIntron(g, read);
				//boolean overlaps=overlapsIntron || overlapsExon;
				if(overlapsExon && overlapsStrand) {rtrn.add(g);}
				//if(overlapsIntron && overlapsStrand) {System.out.println(read.toBED());}
			}
		}
		return rtrn;
	}
	
	
	private boolean isIntergenic(Annotation read, Map<String, IntervalTree<Annotation>> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		if(genes.containsKey(read.getReferenceName())) {
			IntervalTree<? extends Annotation> tree=genes.get(read.getReferenceName());
			Iterator<? extends Annotation> iter=tree.overlappingValueIterator(read.getReferenceStartPosition(), read.getReferenceEndPosition());
			while(iter.hasNext()) {
				Annotation g=iter.next();
				boolean overlapsStrand=overlapsStrand(g,read);
				//boolean overlapsIntron=overlapsIntron(g, read);
				//boolean overlaps=overlapsIntron || overlapsExon;
				if(overlapsStrand) {rtrn.add(g);}
				//if(overlapsIntron && overlapsStrand) {System.out.println(read.toBED());}
			}
		}
		return rtrn.isEmpty();
	}
	
	private boolean overlapsIntron(Annotation gene, Annotation read) {
		if(gene.getName().contains("intron")) {
			return read.overlaps(gene);
		}
		return false;
	}


	private boolean overlapsStrand(Annotation g, Annotation read) {
		return read.getOrientation().equals(g.getOrientation());
	}


	private boolean overlapsExons(Annotation g, Annotation read) {
		//return read.overlaps(g);
		return g.fullyContained(read); //TODO Either all exon or all intron --> deal with spanning seperately
	}


	private Collection<Annotation> getGenes(SAMRecord read, Map<String, IntervalTree<Annotation>> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		if(genes.containsKey(read.getReferenceName())) {
			IntervalTree<? extends Annotation> tree=genes.get(read.getReferenceName());
			Iterator<? extends Annotation> iter=tree.overlappingValueIterator(read.getAlignmentStart(), read.getAlignmentEnd());
			while(iter.hasNext()) {
				Annotation g=iter.next();
				//rtrn.add(g);
				boolean overlaps=overlapsExons(g, read);
				//boolean overlapsExon=overlapsExons(g, read);
				boolean overlapsStrand=overlapsStrand(g, read);
				//boolean overlapsIntron=overlapsIntron(g, read);
				/*if(print && overlapsStrand && !overlapsExon && overlapsIntron) {
					SAMFragment f=new SAMFragment(read);
					System.out.println(f.toBED());
				}*/
				if(overlaps && overlapsStrand) {rtrn.add(g);}
			}
		}
		return rtrn;
	}
	
	
	private boolean overlapsStrand(Annotation g, SAMRecord read) {
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


	private boolean overlapsExons(Annotation g, SAMRecord read) {
		SAMFragment f=new SAMFragment(read);
		//return f.overlaps(g);
		return g.fullyContained(f);
		
	}
	
	private void add(Collection<Annotation> overlappingGenes) {
		for(Annotation gene: overlappingGenes) {
			int score=0;
			if(this.geneCounts.containsKey(gene)) {
				score=geneCounts.get(gene);
			}
			score+=1;
			this.geneCounts.put(gene, score);
		}		
	}


	public Map<Annotation, Integer> getSignificantBins() {
		Map<Annotation, Integer> rtrn=new TreeMap<Annotation, Integer>();
		
		for(Annotation a: this.getEnrichments().keySet()) {
			double enrich=getEnrich(this.getEnrichments(), a);
			int count=get(this.getCounts(), a);
			double p=getP(this.getPvalues(),a);
			/*if(count>0) {
				System.out.println(a.toBED(enrich));
			}*/
			
			if(p<alpha && enrich>threshold && count>minCount) {
				rtrn.put(a, count);
			}
		}
		
		return rtrn;
	}
	
	private double getEnrich(Map<Annotation, Double> enrichments, Annotation a) {
		double rtrn=0;
		if(enrichments.containsKey(a)) {rtrn=enrichments.get(a);}
		return rtrn;
	}
	
	private double getP(Map<Annotation, Double> enrichments, Annotation a) {
		double rtrn=1;
		if(enrichments.containsKey(a)) {rtrn=enrichments.get(a);}
		return rtrn;
	}

	public int getNumberOfWindows() {
		return numberOfWindows;
	}
	

}

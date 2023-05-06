package guttmanlab.core.sharp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.PoissonDistribution;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.ScanStat;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class CLAPAnalysis_AnnotateReads {

	//TODO Use paired interval, second read pileups, or reads
	boolean useFullPair=true; //TODO make this a variable
	
	public CLAPAnalysis_AnnotateReads(File sampleBAM, File inputBAM, Map<String, IntervalTree<Annotation>> genes, int binSize, String save) throws IOException {
		//go through reads and assign to exons, introns, or ambiguous
		SampleValues sample=new SampleValues(sampleBAM, genes, binSize); // make bin counts and bin to Gene maps
		SampleValues input=new SampleValues(inputBAM, genes, binSize);
		
		//Differential scores (sample versus input)
		differentialScores(sample, input, save);
	}

	private void differentialScores(SampleValues sample, SampleValues input, String save) throws IOException {
		//go through all exon, intron, inter -> compute enrichment, pvals sample to input
		FileWriter writer=new FileWriter(save);
		
		for(Annotation bin: sample.exonBinCounts.keySet()) {
			writeScores(writer, bin, sample, input, AssignReads.exon);
		}
		
		for(Annotation bin: sample.intronBinCounts.keySet()) {
			writeScores(writer, bin, sample, input, AssignReads.intron);
		}
		
		for(Annotation bin: sample.intergenicBinCounts.keySet()) {
			writeScores(writer, bin, sample, input, AssignReads.intergenic);
		}
		
		writer.close();
		
	}
	
	private void writeScores(FileWriter writer, Annotation bin, SampleValues sample, SampleValues input, String type) throws IOException {
		int binCount=sample.getBinCount(bin, type);
		int inputCount=Math.max(input.getBinCount(bin, type),1);
		String names="NA";
		double localEnrichment=Double.NaN;
		double localP=Double.NaN;
		
		if(!type.equals(AssignReads.intergenic)) {
			Collection<String> genes=getGenes(bin, sample, input, type);
			inputCount=input.getMaxScore(bin, genes, type);
			localEnrichment=localEnrichment(binCount, bin, sample, type);
			localP=localPVal(binCount, bin, sample, type);
			names=toString(genes);
		}
		
		double diffEnrichment=getEnrichment(binCount, inputCount, sample.getTotalReads(type), input.getTotalReads(type));
		double diffP=getPValue(binCount, inputCount, sample.getTotalReads(type), input.getTotalReads(type));
	
		writer.write(bin.toBlocks()+"\t"+bin.getNumberOfBlocks()+"\t"+bin.getOrientation()+"\t"+names+"\t"+type+"\t"+binCount+"\t"+inputCount+"\t"+localEnrichment+"\t"+diffEnrichment+"\t"+localP+"\t"+diffP+"\n");
		
	}


	private static String toString(Collection<String> exonNames) {
		String rtrn="";
		for(String n: exonNames) {
			rtrn+=n+",";
		}
		rtrn=rtrn.substring(0, rtrn.length()-1);
		return rtrn;
	}
	
	private double localPVal(int score, Annotation bin, SampleValues sample, String type) {
		double rtrn=-1;
		
		Collection<String> genes=sample.getGenes(bin, type);
		for(String gene: genes) {
			int size=sample.getGeneSize(gene, type);
			if(size>0) {
				int geneScore=sample.getGeneScore(gene, type);
				double expected=(double)geneScore/(double)size;
				double pval=ScanStat.getPValue(score, expected, bin.size(), size);
				rtrn=Math.max(rtrn, pval);
			}
		}
		return rtrn;
	}


	private double localEnrichment(int score, Annotation bin, SampleValues sample, String type) {
		Collection<String> genes=sample.getGenes(bin, type);
		
		double rtrn=Double.MAX_VALUE;
		for(String gene: genes) {
			int size=sample.getGeneSize(gene, type);
			if(size>0) {
				int geneScore=sample.getGeneScore(gene, type);
				double expected=(double)geneScore/(double)size;
				double observed=(double)score/(double)bin.size();
				double enrichment=observed/expected;
				rtrn=Math.min(rtrn, enrichment);
			}
		}
		return rtrn;
	}

	
	private Collection<String> getGenes(Annotation bin, SampleValues sample, SampleValues input, String type) {
		if(type.equals(AssignReads.exon)) {return getExonGenes(bin, sample, input);}
		else if(type.equals(AssignReads.intron)) {return getIntronGenes(bin, sample, input);}
		return new TreeSet<String>();
	}

	private Collection<String> getExonGenes(Annotation bin, SampleValues sample, SampleValues input) {
		Collection<String> rtrn=new TreeSet<String>();
		rtrn.addAll(sample.getExonGenes(bin));
		rtrn.addAll(input.getExonGenes(bin));
		return rtrn;
	}
	
	private Collection<String> getIntronGenes(Annotation bin, SampleValues sample, SampleValues input) {
		Collection<String> rtrn=new TreeSet<String>();
		rtrn.addAll(sample.getIntronGenes(bin));
		rtrn.addAll(input.getIntronGenes(bin));
		return rtrn;
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

	
	private static double getEnrichment(double window1, double window2, int total1, int total2) {
		double num=((double)window1)/(double)total1;
		double denom=((double)window2)/(double)total2;
		return num/denom;
	}

	
	
	
	private class SampleValues{
		
		Map<Annotation, Integer> exonBinCounts;
		Map<Annotation, Integer> intronBinCounts;
		Map<Annotation, Integer> intergenicBinCounts;
		
		Map<Annotation, Collection<String>> exonBinToGene;
		Map<Annotation, Collection<String>> IntronBinToGene;
		
		Map<String, Integer> exonGeneCount;
		Map<String, Integer> intronGeneCount;
		
		Map<String, Double> exonSampledCount;
		Map<String, Double> intronSampledCount;
		
		int totalReads;
		int exonTotal;
		int intronTotal;
		int intergenicTotal;
		
		Map<String, Integer> intronSizes;
		Map<String, Integer> exonSizes;
		Map<String, Annotation> genesByName;
		
		boolean sampled;
		
		int numPerm=10;
		
		
		public SampleValues(File bam, Map<String, IntervalTree<Annotation>> genes, int binSize) {
			sampled=false;
			getGeneSizes(genes, binSize);
			assignPairs(bam, genes, binSize); // make bin counts and bin to Gene maps
		}
		
		
		public int getTotalReads(String type) {
			if(type.equals(AssignReads.exon)) {return exonTotal;}
			if(type.equals(AssignReads.intron)) {return intronTotal;}
			if(type.equals(AssignReads.intergenic)) {return intergenicTotal;}
			return getTotalReads();
		}

		public int getBinCount(Annotation bin, String type) {
			int score=0;
			if(type.equals(AssignReads.exon)) {
				if(this.exonBinCounts.containsKey(bin)) {score=this.exonBinCounts.get(bin);}
			}
			else if(type.equals(AssignReads.intron)) {
				if(this.intronBinCounts.containsKey(bin)) {score=this.intronBinCounts.get(bin);}
			}
			else if(type.equals(AssignReads.intergenic)) {
				if(this.intergenicBinCounts.containsKey(bin)) {score=this.intergenicBinCounts.get(bin);}
			}
			return score;
		}

		public int getGeneScore(String gene, String type) {
			int score=0;
			
			if(type.equals(AssignReads.exon)) {
				score=this.exonGeneCount.get(gene);
			}
			else if(type.equals(AssignReads.intron)) {
				score=this.intronGeneCount.get(gene);
			}
			
			return score;
		}

		public int getGeneSize(String gene, String type) {
			if(type.equals(AssignReads.exon)) {return this.exonSizes.get(gene);}
			else if(type.equals(AssignReads.intron)) {return this.intronSizes.get(gene);}
			return 0;
		}

		public Collection<String> getGenes(Annotation bin, String type) {
			if(type.equals(AssignReads.exon)) {return getExonGenes(bin);}
			else if(type.equals(AssignReads.intron)) {return getIntronGenes(bin);}
			return new TreeSet<String>();
		}

		private Collection<String> getExonGenes(Annotation bin) {
			if(this.exonBinToGene.containsKey(bin)) {return this.exonBinToGene.get(bin);}
			return new TreeSet<String>();
		}

		private Collection<String> getIntronGenes(Annotation bin) {
			if(this.IntronBinToGene.containsKey(bin)) {return this.IntronBinToGene.get(bin);}
			return new TreeSet<String>();
		}
		
		private void getGeneSizes(Map<String, IntervalTree<Annotation>> genes, int binSize) {
			this.exonSizes=new TreeMap<String, Integer>();
			this.intronSizes=new TreeMap<String, Integer>();
			this.genesByName=new TreeMap<String, Annotation>();
			
			for(String chr: genes.keySet()) {
				Iterator<Annotation> iter=genes.get(chr).valueIterator();
				while(iter.hasNext()) {
					Annotation gene=iter.next();
					int size=gene.getGenomicLength()-gene.size();
					if(size>0) {
						intronSizes.put(gene.getName(), size);
					}
					exonSizes.put(gene.getName(), gene.size());
					genesByName.put(gene.getName(), gene);
				}
			}
			System.err.println("done with genes");
		}
	
	
	

		public int getMaxScore(Annotation bin, Collection<String> geneNames, String type) {
			if(!sampled) {
				sample(bin.size());
				sampled=true;
			}
			if(type.equals(AssignReads.exon)) {return getMaxExonScore(bin, geneNames);}
			else if(type.equals(AssignReads.intron)) {return getMaxIntronScore(bin, geneNames);}
			return 1;
		}
		
		private void sample(int binSize) {
			this.exonSampledCount=new TreeMap<String, Double>();
			this.intronSampledCount=new TreeMap<String, Double>();
			
			for(String gene: this.exonGeneCount.keySet()) {
				double lambda=((double)exonGeneCount.get(gene)/(double)exonSizes.get(gene))*binSize;
				int numberBins=exonSizes.get(gene)/binSize;
				double perm=getPerms(lambda, numberBins, numPerm);
				exonSampledCount.put(gene, perm);
			}
			
			
			for(String gene: this.intronGeneCount.keySet()) {
				double lambda=((double)intronGeneCount.get(gene)/(double)intronSizes.get(gene))*binSize;
				int numberBins=intronSizes.get(gene)/binSize;
				double perm=getPerms(lambda, numberBins, numPerm);
				intronSampledCount.put(gene, perm);
			}
			
		}
		
		
		private double getPerms(double lambda, int numberBins, int numPerm2) {
			double sum=0;
			PoissonDistribution dist=new PoissonDistribution(lambda);
			
			for(int i=0; i<numPerm2; i++) {
				sum+=max(dist, numberBins);
			}
			return sum/(double)numPerm2;
		}
		
		
		private double max(PoissonDistribution dist, int numberBins) {
			
			int max=0;
			for(int i=0; i<numberBins; i++) {
				int val=dist.sample();
				max=Math.max(val, max);
			}
			return max;
		}


		public int getMaxExonScore(Annotation bin, Collection<String> geneNames) {
			//get counts over bin OR average across genes
			double count1=1;
			if(exonBinCounts.containsKey(bin)) {
				count1=exonBinCounts.get(bin);
			}
			
			double count2=0;
			for(String gene: geneNames) {
				if(exonGeneCount.containsKey(gene)) {
					double score=((double)exonGeneCount.get(gene)/(double)exonSizes.get(gene))*bin.size();
					count2=Math.max(count2, score);
				}
			}
			
			double count3=0;
			for(String gene: geneNames) {
				if(exonSampledCount.containsKey(gene)) {
					double score=exonSampledCount.get(gene);
					count3=Math.max(count3, score);
				}
			}
			
			count2=Math.max(count2, count3);
			
			return (int)Math.ceil(Math.max(count1, count2));
		}
		
		public int getMaxIntronScore(Annotation bin, Collection<String> geneNames) {
			//get counts over bin OR average across genes
			double count1=1;
			if(intronBinCounts.containsKey(bin)) {
				count1=intronBinCounts.get(bin);
			}
			
			double count2=0;
			for(String gene: geneNames) {
				if(intronGeneCount.containsKey(gene)) {
					double score=((double)intronGeneCount.get(gene)/(double)intronSizes.get(gene))*bin.size();
					count2=Math.max(count2, score);
				}
			}
			
			double count3=0;
			for(String gene: geneNames) {
				if(intronSampledCount.containsKey(gene)) {
					double score=intronSampledCount.get(gene);
					count3=Math.max(count3, score);
				}
			}
			
			count2=Math.max(count2, count3);
			
			return (int)Math.ceil(Math.max(count1, count2));
		}

		public int getTotalReads() {
			return totalReads;
		}
		
		
		private void assignPairs(File bam, Map<String, IntervalTree<Annotation>> genes, int binSize) {
			this.exonBinCounts=new TreeMap<Annotation, Integer>();
			this.intronBinCounts=new TreeMap<Annotation, Integer>();
			this.intergenicBinCounts=new TreeMap<Annotation, Integer>();
			this.exonBinToGene=new TreeMap<Annotation, Collection<String>>();
			this.IntronBinToGene=new TreeMap<Annotation, Collection<String>>();
			
			this.totalReads=0;
			this.exonGeneCount=new TreeMap<String, Integer>();
			this.intronGeneCount=new TreeMap<String, Integer>();
			
			SAMFileReader inputReader= new SAMFileReader(bam);
			SAMRecordIterator reads=inputReader.iterator();
			
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
					Pair<String> state=AssignReads.assign(pair, genes);
					add(pair, state, binSize);
					incrementTotal(state);
				}
				else {
					pairedReads.put(name, pair);
				}
				
				
				totalCount++;
				if(totalCount%1000000 ==0){System.err.println(totalCount);}
			}
					
			this.totalReads=totalCount;
			reads.close();
			inputReader.close();
		}
		
		
		private void incrementTotal(Pair<String> pairState) {
			String state=pairState.getValue1();
			if(state.equals(AssignReads.exon)) {exonTotal++;}
			else if(state.equals(AssignReads.intron)) {intronTotal++;}
			else if(state.equals(AssignReads.intergenic)) {intergenicTotal++;}
		}

		private Pair<SAMRecord> update(Pair<SAMRecord> pair, SAMRecord read) {
			if(read.getSecondOfPairFlag()) {pair.setValue2(read);}
			if(read.getFirstOfPairFlag()) {pair.setValue1(read);}
			return pair;
		}
		
		private boolean complete(Pair<SAMRecord> pair) {
			return pair.isComplete();
		}

		private void add(Pair<SAMRecord> pair, Pair<String> pairState, int binSize) {
			String state=pairState.getValue1();
			Collection<String> geneList=parse(pairState.getValue2());
			
			Collection<Annotation> allBins=bin(pair, binSize, useFullPair); //TODO this is where we should do the work to count as reads, paired fragments, or seond read starts
			
			if(state.equals(AssignReads.exon)) {
				Collection<Annotation> genes=getGenes(geneList);
				allBins=bin(pair, binSize, genes);
				increment(this.exonBinCounts, allBins);
				addGenes(this.exonBinToGene, allBins, geneList);
				incrementGene(this.exonGeneCount, geneList);
			}
			else if(state.equals(AssignReads.intron)) {
				increment(this.intronBinCounts, allBins);
				addGenes(this.IntronBinToGene, allBins, geneList);
				incrementGene(this.intronGeneCount, geneList); 
			}
			else if(state.equals(AssignReads.intergenic)) {
				increment(this.intergenicBinCounts, allBins);
			}
			
		}

		private Collection<Annotation> getGenes(Collection<String> geneList) {
			Collection<Annotation> rtrn=new TreeSet<Annotation>();
			for(String name: geneList) {
				Annotation gene=this.genesByName.get(name);
				rtrn.add(gene);
			}
			return rtrn;
		}


		private Collection<Annotation> bin(Pair<SAMRecord> pair, int binSize, Collection<Annotation> genes) {
			Collection<Annotation> rtrn=new TreeSet<Annotation>();
			SAMFragment read1=new SAMFragment(pair.getValue1());
			SAMFragment read2=new SAMFragment(pair.getValue2());
			
			for(Annotation gene: genes) {
				rtrn.addAll(bin(read1, gene, binSize));
				rtrn.addAll(bin(read2, gene, binSize));
			}
			
			return rtrn;
		}

		private Collection<? extends Annotation> bin(SAMFragment read1, Annotation gene, int binSize) {
			int relativeStart=gene.getRelativePositionFrom5PrimeOfFeature(read1.getReferenceStartPosition());
			int relativeEnd=gene.getRelativePositionFrom5PrimeOfFeature(read1.getReferenceEndPosition());
			
			int startIndex=relativeStart/binSize;
			int endIndex=relativeEnd/binSize;
			
			Collection<Annotation> rtrn=new TreeSet<Annotation>();
			
			for(int i=startIndex; i<=endIndex; i++) {
				int newStart=i*binSize;
				int newEnd=newStart+binSize;
				
				Annotation newInterval=gene.trimByRelativePositions(newStart, newEnd);
				//newInterval.setName(gene.getName());
				//newInterval.setOrientation(gene.getOrientation());
				rtrn.add(newInterval);
				//System.out.println(read1.toUCSC()+"\t"+newInterval.toUCSC());	
			}
			
			return rtrn;
		}


		private void incrementGene(Map<String, Integer> exonGeneCount2, Collection<String> geneList) {
			for(String gene: geneList) {
				int count=0;
				if(exonGeneCount2.containsKey(gene)) {count=exonGeneCount2.get(gene);}
				count++;
				exonGeneCount2.put(gene, count);
			}
		}

		/*private Collection<Annotation> bin(Pair<SAMRecord> pair, int binSize) {
			Collection<Annotation> rtrn=new TreeSet<Annotation>();
			rtrn.addAll(SAMFragment.allBins(pair.getValue1(), binSize));
			rtrn.addAll(SAMFragment.allBins(pair.getValue2(), binSize));
			return rtrn;
		}*/
		
		private Collection<Annotation> bin(Pair<SAMRecord> pair, int binSize, boolean useFullPair) {
			Collection<Annotation> rtrn=new TreeSet<Annotation>();
			if(useFullPair) {
				SingleInterval fragment=merge(pair);
				rtrn.addAll(fragment.allBins(binSize));
			}
			
			else {
				rtrn.addAll(SAMFragment.allBins(pair.getValue1(), binSize));
				rtrn.addAll(SAMFragment.allBins(pair.getValue2(), binSize));
			}
			
			return rtrn;
		}

		private SingleInterval merge(Pair<SAMRecord> pair) {
			SAMFragment i1=new SAMFragment(pair.getValue1());
			SAMFragment i2=new SAMFragment(pair.getValue2());
			
			SingleInterval rtrn=new SingleInterval(i1.getReferenceName(), Math.min(i1.getReferenceStartPosition(), i2.getReferenceStartPosition()), Math.max(i1.getReferenceEndPosition(), i2.getReferenceEndPosition()));
			rtrn.setOrientation(i1.getOrientation());
			return rtrn;
		}


		private void addGenes(Map<Annotation, Collection<String>> exonBinToGene2, Collection<Annotation> allBins, Collection<String> geneList) {
			for(Annotation bin: allBins) {
				Collection<String> list=new TreeSet<String>();
				if(exonBinToGene2.containsKey(bin)) {list=exonBinToGene2.get(bin);}
				list.addAll(geneList);
				exonBinToGene2.put(bin, list);
			}
		}

		private void increment(Map<Annotation, Integer> exonBinCounts2, Collection<Annotation> allBins) {
			for(Annotation bin: allBins) {
				int count=0;
				if(exonBinCounts2.containsKey(bin)) {
					count=exonBinCounts2.get(bin);
				}
				count++;
				exonBinCounts2.put(bin, count);
			}
		}

		private Collection<String> parse(String geneNames) {
			Collection<String> rtrn=new TreeSet<String>();
			String[] tokens=geneNames.split(",");
			for(int i=0; i<tokens.length; i++) {rtrn.add(tokens[i]);}
			return rtrn;
		}
		
	}


	public static void main(String[] args) throws IOException {
		if(args.length>4) {
			File sample=new File(args[0]);
			File input=new File(args[1]);
			Map<String, IntervalTree<Annotation>> genes=BEDFileIO.loadTreeAnnotation(args[2]);
			int binSize=Integer.parseInt(args[3]);
			String save=args[4];
			new CLAPAnalysis_AnnotateReads(sample, input, genes, binSize, save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=sample bam \n args[1]=input bam \n args[2]=genes \n args[3]=bin size \n args[4]=save";
}

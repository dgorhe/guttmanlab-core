package guttmanlab.core.splicing.speckle;

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

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.EnumeratedDistribution;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SplicingRatiosPerIsoform {
	
	/*private Map<SingleInterval, Integer> intronBinCounts;
	int binResolution=10;*/

	public SplicingRatiosPerIsoform(File bamFile, Collection<Gene> genes, String save, boolean debug) throws IOException{
		Map<String, IntervalTree<Gene>> geneTree=makeTree(genes);
		
		Map<Gene, Collection<Gene>> genesToJunctions=makeGenesToJunctions(genes);
				
		Map<String, IntervalTree<Gene>> junctionTree=makeTree(genesToJunctions);
		
		
		
		Map<Gene, Pair<Integer>> countsByGene=readsOverlappingSpliceJunction(bamFile, geneTree, save, debug);
		Map<Gene, Pair<Integer>> countsByJunction=readsOverlappingSpliceJunction(bamFile, junctionTree, save, debug);
		
		write(save, countsByGene, countsByJunction, genesToJunctions);
	}
	
	
	
	public SplicingRatiosPerIsoform(File bamFile, GTFToJunctions gtf, String save, boolean debug) throws IOException{
		Map<String, IntervalTree<Gene>> geneTree=gtf.getJunctionTree();
		Map<Gene, Pair<Integer>> countsByGene=readsOverlappingSpliceJunction(bamFile, geneTree, save, debug);
		//writeIntron(intronBinCounts, save+".intron.bedgraph");
		write(save, countsByGene, gtf);
	}
	
	private void writeIntron(Map<SingleInterval, Integer> intronBinCounts2, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: intronBinCounts2.keySet()) {
			int count=intronBinCounts2.get(r);
			writer.write(r.toBedgraph(count)+"\n");
		}
		
		writer.close();
	}



	private void write(String save, Map<Gene, Pair<Integer>> countsByGene, GTFToJunctions gtf) throws IOException {
		FileWriter writer=new FileWriter(save);
		FileWriter writerBG=new FileWriter(save+".bedgraph");
		FileWriter writerLowBG=new FileWriter(save+".low.bedgraph");
		FileWriter writerHighBG=new FileWriter(save+".high.bedgraph");
		
		double[] vals=new double[101];
		double interval=1.0/((double)vals.length-1);
		
		for(int i=0; i<vals.length; i++) {
			vals[i]=i*interval;
		}
		
		for(Gene gene: countsByGene.keySet()) {
			Pair<Integer> pair=countsByGene.get(gene);
			
			
			int sum=pair.getValue1()+pair.getValue2();
			double r=(double)pair.getValue1()/(double)sum;
			double[] probs=getProbabilities(sum, pair.getValue1(), vals);
			if(Statistics.max(probs)>0) {
				
				//EnumeratedDistribution<Double> dist=makeDist(probs, vals);
				SingleInterval mid=gene.getIntrons().iterator().next().getSingleInterval().getMidPoint();
				mid=new SingleInterval(mid.getReferenceName(), mid.getReferenceStartPosition()-50, mid.getReferenceEndPosition()+50);
				double valAtMax=getMax(vals, probs);
				Pair<Double> range=getRange(vals, probs, 0.05);
				int intronLength=gene.getGenomicLength()-gene.size();
				writer.write(gene.getName()+"\t"+gene.toUCSC()+"\t"+gene.size()+"\t"+intronLength+"\t"+pair.getValue1()+"\t"+pair.getValue2()+"\t"+(valAtMax*sum)+"\t"+(range.getValue1()*sum)+"\t"+(range.getValue2()*sum));
				/*for(int i=0; i<numPerm; i++) {
					writer.write("\t"+dist.sample());
				}*/
				writer.write("\n");
				
				writerBG.write(mid.toBedgraph(r)+"\n");
				writerLowBG.write(mid.toBedgraph(range.getValue1())+"\n");
				writerHighBG.write(mid.toBedgraph(range.getValue2())+"\n");
			}
		}
		writer.close();
		writerBG.close();
		writerLowBG.close();
		writerHighBG.close();
	}

	

	/*private List<Double> getValues(SingleInterval intron, Map<SingleInterval, Integer> intronBinCounts2) {
		List<Double> rtrn=new ArrayList<Double>();
		
		Collection<SingleInterval> bins=intron.allBins(binResolution);
		for(SingleInterval bin: bins) {
			if(intron.fullyContained(bin)) {
				System.out.println(bin.toShortBED());
				if(intronBinCounts2.containsKey(bin)) {rtrn.add((double)intronBinCounts2.get(bin));}
				else {rtrn.add(0.0);}
			}
		}
		
		return rtrn;
	}*/



	private EnumeratedDistribution<Double> makeDist(double[] probs, double[] vals) {
		List<org.apache.commons.math3.util.Pair<Double, Double>> list=new ArrayList<org.apache.commons.math3.util.Pair<Double, Double>>();
		
		for(int i=0; i<probs.length; i++) {
			org.apache.commons.math3.util.Pair<Double, Double> pair=new org.apache.commons.math3.util.Pair<Double, Double>(vals[i], probs[i]);
			list.add(pair);
		}
		
		EnumeratedDistribution<Double> rtrn=new EnumeratedDistribution<Double>(list);
		return rtrn;
	}



	private Pair<Double> getRange(double[] vals, double[] probs, double ci) {
		double max=Statistics.max(probs);
				
		int minIndex=getMinIndex(probs, max, ci);
		int maxIndex=getMaxIndex(probs, max, ci);
		
		Pair<Double> rtrn=new Pair<Double>();
		
		rtrn.setValue1(vals[minIndex]);
		rtrn.setValue2(vals[maxIndex]);
		
		return rtrn;
	}

	private int getMaxIndex(double[] probs, double max, double ci) {
		for(int i=probs.length-1; i>=0; i--) {
			double rp=probs[i]/max;
			if(rp>ci) {
				return i;
			}
		}
		return probs.length-1;
	}


	private int getMinIndex(double[] probs, double max, double ci) {
		for(int i=0; i<probs.length; i++) {
			double rp=probs[i]/max;
			if(rp>ci) {
				return i;
			}
		}
		return 0;
	}



	private double getMax(double[] vals, double[] probs) {
		double max=0;
		int index=0;
		for(int i=0; i<probs.length; i++) {
			if(probs[i]>max) {
				max=probs[i];
				index=i;
			}
		}
		return vals[index];
	}



	private double[] getProbabilities(int count, int observedSpliced, double[] testFraction) {
		double[] rtrn=new double[testFraction.length];
		for(int i=0; i<testFraction.length; i++) {
			rtrn[i]=sampleProbability(testFraction[i], count, observedSpliced);
		}
		
		
		
		return rtrn;
	}



	private double sampleProbability(double p, int count, int observedSpliced) {
		
		//TODO binomial
		//int n=count;
		int k=observedSpliced;
		BinomialDistribution dist=new BinomialDistribution(count, p);
		return dist.probability(k);
		
		/*double[] randomFract=new double[numPerm];
		for(int i=0; i<numPerm; i++) {
			double r=sample(p, count);
			randomFract[i]=r;
		}
		
		EmpiricalDistribution dist=new EmpiricalDistribution();
		dist.load(randomFract);
		
		return dist.density(observedRatio);*/
	}



	private double sample(double p, int count) {
		double spliced=0;
		for(int i=0; i<count; i++) {
			double rand=Math.random();
			if(rand<p) {spliced++;}
		}
		double frac=spliced/(double)count;
		return frac;
	}



	private Map<Gene, Collection<Gene>> makeGenesToJunctions(Collection<Gene> genes) {
		Map<Gene, Collection<Gene>> rtrn=new TreeMap<Gene, Collection<Gene>>();
		
		for(Gene gene: genes) {
			rtrn.put(gene, gene.getJunctions());
		}
		
		return rtrn;
	}





	private void write(String save, Map<Gene, Pair<Integer>> countsByGene, Map<Gene, Pair<Integer>> countsByJunction, Map<Gene, Collection<Gene>> genesToJunctions) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: countsByGene.keySet()) {
			Pair<Integer> pair=countsByGene.get(gene);
			double ratio=getExpectedRatio(gene);
			double r2=(double)pair.getValue1()/(double)(pair.getValue1()+pair.getValue2());
			writer.write(gene.getName()+"\t"+gene.toUCSC()+"\t"+ratio+"\t"+pair.getValue1()+"\t"+pair.getValue2());
			
			
			Collection<Gene> junctions=genesToJunctions.get(gene);
			for(Gene junction: junctions) {
				Pair<Integer> jScore=new Pair<Integer>(0,0);
				if(countsByJunction.containsKey(junction)) {
					jScore=countsByJunction.get(junction);
				}
					int v1=jScore.getValue1();
					int v2=jScore.getValue2();
					double ratio2=(double)v1/(double)(v1+v2);
					writer.write("\t"+junction.toUCSC()+"\t"+v1+"\t"+v2);
				}
				writer.write("\n");
			}
		writer.close();
		
	}





	private Map<String, IntervalTree<Gene>> makeTree(Map<Gene, Collection<Gene>> genesToJunctions) {
		Collection<Gene> junctions=new TreeSet<Gene>();
		
		for(Gene g: genesToJunctions.keySet()) {
			junctions.addAll(genesToJunctions.get(g));
		}
		
		return makeTree(junctions);
	}





	private void write(String save, Map<Gene, Pair<Integer>> countsByJunction) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene junction: countsByJunction.keySet()) {
			Pair<Integer> pair=countsByJunction.get(junction);
			double ratio=getExpectedRatio(junction);
			writer.write(junction.getName()+"\t"+junction.toUCSC()+"\t"+ratio+"\t"+pair.getValue1()+"\t"+pair.getValue2()+"\n");
		}
		
		writer.close();
	}

	private double getExpectedRatio(Gene g) throws IOException {
		double ratio=(double)g.size()/(double)g.getGenomicLength();
		return ratio;
	}

	
	
	private String assignState(SAMFragment read, Gene gene) {
		//If read is spliced, check if splice junction matches gene junction, if so --> spliced
		if(read.isSpliced()) {
			if(junctionMatch(read, gene)) {return "spliced";}
			else {return "amb";}
		}
		
		for(Annotation intron: gene.getIntrons()) {
			if(read.overlaps(intron)) {
				//TODO Trim read by 3 bases
				return "unspliced";
			}
		}
		//return "spliced";
		return "amb";
	}

	private boolean junctionMatch(SAMFragment read, Gene gene) {
		if(read.getOrientation().equals(gene.getOrientation())) {
			Collection<SingleInterval> readJunctions=getJunctions(read);
			Collection<SingleInterval> geneJunctions=getJunctions(gene);
			for(SingleInterval readJunction: readJunctions) {
				if(geneJunctions.contains(readJunction)){return true;}
			}	
		}
		return false;
	}





	private Map<Gene, Pair<Integer>> readsOverlappingSpliceJunction(File bamFile, Map<String, IntervalTree<Gene>> geneTree, String save, boolean debug) {
		//this.intronBinCounts=new TreeMap<SingleInterval, Integer>();
		Map<Gene, Pair<Integer>> rtrn=new TreeMap<Gene, Pair<Integer>>();
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileWriter exonWriter=null;
		SAMFileWriter intronWriter=null;
		SAMFileWriter ambWriter=null;
		
		if(debug) {
			exonWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save+".exon.bam"));
			intronWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save+".intron.bam"));
			ambWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save+".amb.bam"));
		}
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment read=new SAMFragment(record);
			
			Collection<Gene> genes=findGene(read, geneTree);
			
			Collection<Gene> overlappingGenes=overlapsGene(read, genes);
			
			
			/*if(!read.isSpliced()) {
				add(read.getSingleInterval().allBins(binResolution), intronBinCounts);
			}*/
			
			for(Gene g: overlappingGenes) {
				Pair<Integer> pair=assign(g, read);
				if(debug) {
					if(pair.getValue1()==1) {exonWriter.addAlignment(record);}
					else if(pair.getValue2()==1) {intronWriter.addAlignment(record);}
					else {ambWriter.addAlignment(record);}
				}
				
				if(rtrn.containsKey(g)) {
					Pair<Integer> old=rtrn.get(g);
					int val1=pair.getValue1()+old.getValue1();
					int val2=pair.getValue2()+old.getValue2();
					/*if(pair.getValue2()==1) {
						add(read.getSingleInterval().allBins(binResolution), intronBinCounts);
					}*/
					
					pair=new Pair<Integer>(val1, val2);
				}
				
				rtrn.put(g, pair);
			}
			
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		if(debug) {
			exonWriter.close();
			intronWriter.close();
			ambWriter.close();
		}
		reader.close();
		reads.close();
		return rtrn;
	}
	
	
	
	private void add(Collection<SingleInterval> allBins, Map<SingleInterval, Integer> intronBinCounts2) {
		for(SingleInterval bin: allBins) {
			int count=0;
			if(intronBinCounts2.containsKey(bin)) {count=intronBinCounts2.get(bin);}
			count++;
			intronBinCounts2.put(bin, count);
		}
	}



	private Pair<Integer> assign(Gene g, SAMFragment read) {
		
		int exonCount=0;
		int intronCount=0;
		
		String state=assignState(read, g);
		if(state.equals("spliced")) {
			exonCount++;
		}
		if(state.equals("unspliced")) {
			intronCount++;
		}
		Pair<Integer> rtrn=new Pair<Integer>(exonCount, intronCount);
		
		return rtrn;
		
		
	}

	
	
	private Collection<Gene> overlapsGene(SAMFragment read, Collection<Gene> genes) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		for(Gene gene: genes){
			if(overlaps(gene, read)) {rtrn.add(gene);}
		}
		return rtrn;
		
	}
	
	private boolean overlaps(Gene gene, SAMFragment read) {
		if(read.getOrientation().equals(gene.getOrientation())) {
			
			Collection<SingleInterval> readJunctions=getJunctions(read);
			Collection<SingleInterval> geneJunctions=getJunctions(gene);
			for(SingleInterval readJunction: readJunctions) {
				if(geneJunctions.contains(readJunction)){return true;}
			}
			
			if(read.getReferenceName().equals(gene.getReferenceName()) && read.getReferenceStartPosition()>=gene.getReferenceStartPosition() && read.getReferenceEndPosition()<=gene.getReferenceEndPosition()) {
				return true;
			}
		}
		return false;
	}


	private Collection<SingleInterval> getJunctions(Annotation gene) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		Collection<Annotation> blocks=gene.getIntrons();
		
		for(Annotation intron: blocks) {
			SingleInterval intronSI=new SingleInterval(intron.getReferenceName(), intron.getReferenceStartPosition(), intron.getReferenceEndPosition());
			intronSI.setOrientation(gene.getOrientation());
			rtrn.add(intronSI);
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
	
	
	
	private static void writeTable(File[] bams, Map<Gene, Pair<Integer>>[] scores, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Gene> allGenes=new TreeSet<Gene>();
		
		writer.write("Gene Name \t Coordinates");
		for(int i=0; i<scores.length; i++) {
			allGenes.addAll(scores[i].keySet());
			writer.write("\t"+bams[i].getName());
		}
		writer.write("\n");
		
		for(Gene gene: allGenes) {
			writer.write(gene.getName()+"\t"+gene.toUCSC());
			for(int i=0; i<scores.length; i++) {
				double ratio=getRatio(scores[i], gene);
				writer.write("\t"+ratio);
			}
			writer.write("\n");
		}
		
		
		writer.close();
	}



	
	private static double getRatio(Map<Gene, Pair<Integer>> map, Gene gene) {
		double ratio=0;
		if(map.containsKey(gene)) {
			Pair<Integer> pair=map.get(gene);
			double sum=pair.getValue1()+pair.getValue2();
			ratio=(double)pair.getValue1()/sum;
		}
		return ratio;
	}

	
	private static List<File> getBAMs(String string) {
		List<File> rtrn=new ArrayList<File>();
		
		File[] files=new File(string).listFiles();
		for(int i=0; i<files.length; i++) {
			if(files[i].getName().endsWith(".bam")) {rtrn.add(files[i]);}
	
		}
		
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File bam=new File(args[0]);
			GTFToJunctions gtf=new GTFToJunctions(new File(args[1]));
			String save=args[2];
			SplicingRatiosPerIsoform sample1=new SplicingRatiosPerIsoform(bam, gtf, save, false);
		}
		else{System.err.println(usage);}
	}
	
	


	

	//TODO filter junctions that are less than 1% of max junctions in gene






	static String usage=" args[0]=bam file (time point) \n args[1]=GTF \n args[3]=save";
	
	
}

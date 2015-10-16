package guttmanlab.core.examples;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.MappedFragment;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AbstractAnnotationCollection;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMFragmentCollectionFactory;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.math.ScanStat;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

import jsc.distributions.Binomial;
import net.sf.samtools.util.CloseableIterator;

public class ComputeLocalizedEnrichment {
	
	static int windowSize=100;
	
	/*public ComputeLocalizedEnrichment(File bamFile, File geneFile, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		//convert BAM to cDNA
		AbstractAnnotationCollection<? extends MappedFragment> bam=BAMFragmentCollectionFactory.createFromBam(bamFile, false);
		AnnotationCollection<Gene> genes=BEDFileIO.loadFromFile(geneFile, bam.getReferenceCoordinateSpace());
		
		CloseableIterator<Gene> iter=genes.sortedIterator();
		
		//scan cDNA
		while(iter.hasNext()){
			Gene gene=iter.next();
			int readsPerGene=bam.numOverlappers(gene, true);
			double expected=(double)readsPerGene/(double)gene.size();
			CloseableIterator<? extends PopulatedWindow<? extends MappedFragment>> windows=bam.getPopulatedWindows(gene, windowSize, 1);
			while(windows.hasNext()){
				PopulatedWindow w=windows.next();
				double observed=(double)w.getNumberOfAnnotationsInWindow()/(double)windowSize;
				double score=observed/expected;
				writer.write(w.getReferenceName()+"\t"+w.getReferenceStartPosition()+"\t"+w.getReferenceEndPosition()+"\t"+score+"\n");
			}
		}
		writer.close();
	}*/
	
	public ComputeLocalizedEnrichment(File sampleFile, File inputFile, File geneFile, String save) throws IOException{
		
	}
	
	public static void localEnrichments(File sampleFile, File inputFile, File geneFile, String save, String chrToUse, boolean exonsOnly)throws Exception{
		//AbstractAnnotationCollection<? extends MappedFragment> bam=BAMFragmentCollectionFactory.createFromBam(sampleFile, true);
		//AbstractAnnotationCollection<? extends MappedFragment> input=BAMFragmentCollectionFactory.createFromBam(inputFile, true);
		
		BAMPairedFragmentCollection bam = new BAMPairedFragmentCollection(sampleFile);
		BAMPairedFragmentCollection input = new BAMPairedFragmentCollection(inputFile);
		//AnnotationCollection<Gene> genes=BEDFileIO.loadFromFile(geneFile, bam.getReferenceCoordinateSpace());
		BEDFileIO io =  new BEDFileIO("/Users/cburghard/Downloads/sizes");
		String FeatureFile = "/Users/cburghard/Downloads/RefSeq.bed";
		AnnotationCollection<Gene> genes = io.loadFromFile(FeatureFile);
		CloseableIterator<Gene> iter = genes.sortedIterator();
		
		FileWriter writer=new FileWriter(save);
		
		System.err.println("initialized");
		
		/* TODO Add back
		int inputTotal=input.getNumAnnotations();
		int sampleTotal=bam.getNumAnnotations();
		int genomeSize=(int)bam.getReferenceCoordinateSpace().getTotalReferenceLength();
		*/
		
		System.err.println("computed bulk stats");
		
		//FeatureCollection<Gene> geneCollection=genes.get(chrToUse);
		
		int counter=0;
		//for(Gene gene: geneCollection){
		while(iter.hasNext()){
			Gene gene = iter.next();
			//Annotation gene=new SingleInterval(geneRegion.getReferenceName(), geneRegion.getReferenceStartPosition(), geneRegion.getReferenceEndPosition(), geneRegion.getOrientation());
			
			
			int readsPerGene=bam.numOverlappers(gene, true);
			int readsPerGeneInput=input.numOverlappers(gene, true);
			double lambda=(double)readsPerGene/(double)gene.size();
			CloseableIterator<? extends PopulatedWindow<? extends MappedFragment>> windows=bam.getPopulatedWindows(gene, windowSize, 1);
			
			//System.err.println(gene.getName()+" "+gene.toUCSC()+" "+counter+" "+geneCollection.size()+" "+readsPerGeneInput+" "+" "+lambda);
			
			
			/*while(windows.hasNext()){
				//TODO This is where we should compute significance
				PopulatedWindow w=windows.next();
						
				double controlCount=input.numOverlappers(w, false);
				double sampleCount=w.getNumberOfAnnotationsInWindow();
						
				double p1=ScanStat.scanPVal(w.getNumberOfAnnotationsInWindow(), lambda, w.size(), gene.size()); //Window enrichment over expected gene coverage
				double scanPLocal=ScanStat.scanPValue(controlCount, sampleCount, readsPerGeneInput, readsPerGene, w.size(), gene.size());
						
				double scanP=ScanStat.scanPValue(controlCount, sampleCount, inputTotal, sampleTotal, w.size(), genomeSize);
						
				//double enrichment=(sampleCount/controlCount)/(readsPerGene/readsPerGeneInput);
				double enrichment=(sampleCount/w.size())/((double)readsPerGene/gene.size());
						
				if(p1<0.01 && scanP<0.01 && scanPLocal<0.01){
					writer.write(w.getReferenceName()+"\t"+w.getReferenceStartPosition()+"\t"+w.getReferenceEndPosition()+"\t"+enrichment+"\n");
				}
			}*/
			
			counter++;
		}
		writer.close();
	}
	
	
	/*public static void localEnrichments(File sampleFile, File inputFile, File geneFile, String save, String chrToUse, boolean exonsOnly)throws IOException{
		AbstractAnnotationCollection<? extends MappedFragment> bam=BAMFragmentCollectionFactory.createFromBam(sampleFile, true);
		AbstractAnnotationCollection<? extends MappedFragment> input=BAMFragmentCollectionFactory.createFromBam(inputFile, true);
		
		Map<String, FeatureCollection<Gene>> genes=BEDFileIO.loadFromFileByReferenceName(geneFile.getAbsolutePath(), bam.getReferenceCoordinateSpace());
		
		FileWriter writer=new FileWriter(save);
		
		System.err.println("initialized");
		
		
		int inputTotal=input.getNumAnnotations();
		int sampleTotal=bam.getNumAnnotations();
		int genomeSize=(int)bam.getReferenceCoordinateSpace().getTotalReferenceLength();
		
		
		System.err.println("computed bulk stats");
		
		FeatureCollection<Gene> geneCollection=genes.get(chrToUse);
		
		int counter=0;
		for(Gene geneRegion: geneCollection){
			
			Collection<Annotation> scan=new ArrayList<Annotation>();
			scan.add(geneRegion);
			if(!exonsOnly){
				Collection<SingleInterval> introns=geneRegion.getIntronSet();
				scan.addAll(introns);
			}
			
			Annotation region=new SingleInterval(geneRegion.getReferenceName(), geneRegion.getReferenceStartPosition(), geneRegion.getReferenceEndPosition(), geneRegion.getOrientation());
			
			System.err.println(geneRegion.getName()+" "+counter+" "+geneCollection.size()+" "+input.computeScanPValue(geneRegion) +" "+input.computeScanPValue(region));
			
			for(Annotation gene: scan){
				double expressed=input.computeScanPValue(gene);
				
				//System.err.println(gene.toUCSC()+" "+expressed);
				if(expressed<0.01){
					int readsPerGene=bam.numOverlappers(gene, true);
					int readsPerGeneInput=input.numOverlappers(gene, true);
					double lambda=(double)readsPerGene/(double)gene.size();
					CloseableIterator<? extends PopulatedWindow<? extends MappedFragment>> windows=bam.getPopulatedWindows(gene, windowSize, 1);
					while(windows.hasNext()){
						//TODO This is where we should compute significance
						PopulatedWindow w=windows.next();
						
						double controlCount=input.numOverlappers(w, false);
						double sampleCount=w.getNumberOfAnnotationsInWindow();
						
						double p1=ScanStat.scanPVal(w.getNumberOfAnnotationsInWindow(), lambda, w.size(), gene.size()); //Window enrichment over expected gene coverage
						double scanPLocal=ScanStat.scanPValue(controlCount, sampleCount, readsPerGeneInput, readsPerGene, w.size(), gene.size());
						
						double scanP=ScanStat.scanPValue(controlCount, sampleCount, inputTotal, sampleTotal, w.size(), genomeSize);
						
						//double enrichment=(sampleCount/controlCount)/(readsPerGene/readsPerGeneInput);
						double enrichment=(sampleCount/w.size())/((double)readsPerGene/gene.size());
						
						if(p1<0.01 && scanP<0.01 && scanPLocal<0.01){
							writer.write(w.getReferenceName()+"\t"+w.getReferenceStartPosition()+"\t"+w.getReferenceEndPosition()+"\t"+enrichment+"\n");
						}
					}
				}
			}
			counter++;
		}
		writer.close();
		
	}*/
	
	public static void geneLevelEnrichments(File sampleFile, File inputFile, File geneFile, String save) throws Exception{
		//AbstractAnnotationCollection<? extends MappedFragment> bam=BAMFragmentCollectionFactory.createFromBam(sampleFile, true);
		//AbstractAnnotationCollection<? extends MappedFragment> input=BAMFragmentCollectionFactory.createFromBam(inputFile, true);
		BAMPairedFragmentCollection bam = new BAMPairedFragmentCollection(sampleFile);
		BAMPairedFragmentCollection input = new BAMPairedFragmentCollection(inputFile);
		//AnnotationCollection<Gene> genes=BEDFileIO.loadFromFile(geneFile, bam.getReferenceCoordinateSpace());
		BEDFileIO io =  new BEDFileIO("/Users/cburghard/Downloads/sizes");
		String FeatureFile = "/Users/cburghard/Downloads/RefSeq.bed";
		CloseableIterator<Gene> iter = io.loadFromFile(FeatureFile).sortedIterator();
		
		double p=(double)bam.getNumAnnotations();///(bam.getNumAnnotations()+input.getNumAnnotations());
		
		System.err.println(p);
		FileWriter writer=new FileWriter(save);
		
		//CloseableIterator<Gene> iter=genes.sortedIterator();
		
		int counter=0;
		while(iter.hasNext()){
			Gene gene=iter.next();
			int k=bam.numOverlappers(gene, true);
			int n= k + input.numOverlappers(gene, true);
			
			double p1=1;//bam.computeScanPValue(gene);
			double p2=1;//input.computeScanPValue(gene);
			double scanP=1;//ScanStat.scanPValue((double)(n-k), (double)k, (double)input.getNumAnnotations(), (double)bam.getNumAnnotations(), (double)gene.size(), (int)bam.getReferenceCoordinateSpace().getTotalReferenceLength());
			double enrichment=((double)k/(double)n)/p;
			
			//System.err.println(k+"\t"+((double)k/n)+"\t"+p1+"\t"+p2+"\t"+scanP);
			
			if(n>50 && p1<.01 && p2<0.01 && scanP<0.01 && enrichment>1.5){
				writer.write(gene.getReferenceName()+"\t"+gene.getReferenceStartPosition()+"\t"+gene.getReferenceEndPosition()+"\t"+enrichment+"\n");
			}
			counter++;
			if(counter%10000 ==0){System.err.println(counter);}
		}
		writer.close();
	}
	
	public static void main(String[] args) throws Exception{
		if(args.length>5){
			//geneLevelEnrichments(new File(args[0]), new File(args[1]), new File(args[2]), args[3]);
			localEnrichments(new File(args[0]), new File(args[1]), new File(args[2]), args[3], args[4], new Boolean(args[5]));
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=sample BAM file \n args[1]=input BAM file \n args[2]=gene BED file \n args[3]=save \n args[4]=chrToUse \n args[5]=exons only (true=exons, false=exons+introns";
	
}

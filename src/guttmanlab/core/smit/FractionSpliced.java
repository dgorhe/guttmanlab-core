package guttmanlab.core.smit;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import net.sf.samtools.util.CloseableIterator;

public class FractionSpliced {
	
	Map<Gene, Integer> totalCounts;
	Map<Gene, Integer> intronCounts;
	

	public FractionSpliced(BAMPairedFragmentCollection bamData, AnnotationCollection<Gene> genes) throws IOException{
		//Map<Gene, Double> map=new TreeMap<Gene, Double>();
		System.err.println("started");
		
		totalCounts=new TreeMap<Gene, Integer>();
		intronCounts=new TreeMap<Gene, Integer>();
		
		CloseableIterator<PairedMappedFragment<SAMFragment>> iter=bamData.sortedIterator();
		
		int counter=0;
		
		while(iter.hasNext()){
			PairedMappedFragment<SAMFragment> fragment=iter.next();
			//TODO Classify into intron/exon
			CloseableIterator<Gene> geneIter=genes.sortedIterator(fragment, true);
			
			
			while(geneIter.hasNext()){
				Gene gene=geneIter.next();
				if(gene.getOrientation().equals(fragment.getOrientation())){
					boolean intron=overlapsIntron(fragment, gene);
					//System.err.println(fragment.getRead1().toUCSC(fragment.getOrientation())+" "+fragment.getRead2().toUCSC(fragment.getOrientation())+" "+gene.getName()+" "+gene.getOrientation()+" "+intron);
					int totalCount=0;
					int intronCount=0;
					if(totalCounts.containsKey(gene)){totalCount=totalCounts.get(gene);}
					if(intronCounts.containsKey(gene)){intronCount=intronCounts.get(gene);}
					totalCount++;
					if(intron){intronCount++;}
					intronCounts.put(gene, intronCount);
					totalCounts.put(gene, totalCount);
				}
			}
			counter++;
			if(counter%100000==0){System.err.println(counter+" "+fragment.toUCSC());}
		}
		
		
		
		//For each gene, calculate percent spliced
		/*for(Gene gene: genes){
			System.err.println(gene.getName());
			int intronCount=0;
			int totalCount=0;
			int count=bamData.numOverlappers(gene.getGenomicRegion(), true);
			System.err.println(gene.getName()+" "+gene.toUCSC()+" "+count);
			CloseableIterator<SAMFragment> iter=bamData.sortedIterator(gene.getGenomicRegion(), true);
			while(iter.hasNext()){
				SAMFragment read=iter.next();
				if(overlapsIntron(read, gene)){intronCount++;}
				totalCount++;
			}
			double ratio=(double)intronCount/(double)totalCount;
			writer.write(gene.getName()+"\t"+gene.getGenomicLength()+"\t"+intronCount+"\t"+totalCount+"\t"+ratio);
			System.err.println(gene.getName()+"\t"+gene.getGenomicLength()+"\t"+intronCount+"\t"+totalCount+"\t"+ratio);
		}*/
		
	}
	
	private void write(String save, Map<Gene, Integer> totalCounts, Map<Gene, Integer> intronCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: totalCounts.keySet()){
			int total=totalCounts.get(gene);
			int intron=0;
			if(intronCounts.containsKey(gene)){
				intron=intronCounts.get(gene);
			}
			writer.write(gene.getName()+"\t"+gene.toUCSC()+"\t"+gene.getGenomicLength()+"\t"+intron+"\t"+total+"\n");
		}
		
		writer.close();
	}

	private boolean overlapsIntron(PairedMappedFragment<SAMFragment> fragment, Gene gene) {
		//If read1 or read2 is in intron
		Annotation introns=new BlockedAnnotation(gene.getIntrons(), gene.getName()+"_introns");
		if(introns.overlaps(fragment.getRead1()) || introns.overlaps(fragment.getRead2())){return true;}
		return false;
	}

	public static void main(String[] args)throws IOException{
		if(args.length>3){
			BAMPairedFragmentCollection data1=new BAMPairedFragmentCollection(new File(args[0]));
			BAMPairedFragmentCollection data2=new BAMPairedFragmentCollection(new File(args[1]));
			//BAMSingleReadCollection data=new BAMSingleReadCollection(new File(args[0]));
			AnnotationCollection<Gene> genes=BEDFileIO.loadFromFile((args[2]));
			String save=args[3];
			
			
			
			FractionSpliced WT=new FractionSpliced(data1, genes);
			FractionSpliced KO=new FractionSpliced(data2, genes);
			
			write(save, WT, KO);
			writeBED(save+".bed", WT, KO);
		}
		else{System.err.println(usage);}
		
	}
	
	static String usage=" args[0]=WT bam \n args[1]=KO bam \n args[2]=BED file (genes) \n args[3]=save";

	private static void writeBED(String save, FractionSpliced WT, FractionSpliced KO) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Gene> allGenes=new TreeSet<Gene>();
		allGenes.addAll(WT.totalCounts.keySet());
		allGenes.addAll(KO.totalCounts.keySet());
		
		for(Gene gene: allGenes){
			int wtTotal=0;
			int wtIntron=0;
			int koTotal=0;
			int koIntron=0;
			
			if(WT.totalCounts.containsKey(gene)){wtTotal=WT.totalCounts.get(gene);}
			if(WT.intronCounts.containsKey(gene)){wtIntron=WT.intronCounts.get(gene);}
			if(KO.totalCounts.containsKey(gene)){koTotal=KO.totalCounts.get(gene);}
			if(KO.intronCounts.containsKey(gene)){koIntron=KO.intronCounts.get(gene);}
			
			double ratio1=(double)wtIntron/(double)wtTotal;
			double ratio2=(double)koIntron/(double)koTotal;
			double diff=ratio2-ratio1;
			
			gene.setName("SR="+diff);
			writer.write(gene.toBED()+"\n");
		}
		
		
		writer.close();
		
	}

	private static void write(String save, FractionSpliced WT, FractionSpliced KO) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Gene> allGenes=new TreeSet<Gene>();
		allGenes.addAll(WT.totalCounts.keySet());
		allGenes.addAll(KO.totalCounts.keySet());
		
		for(Gene gene: allGenes){
			int wtTotal=0;
			int wtIntron=0;
			int koTotal=0;
			int koIntron=0;
			
			if(WT.totalCounts.containsKey(gene)){wtTotal=WT.totalCounts.get(gene);}
			if(WT.intronCounts.containsKey(gene)){wtIntron=WT.intronCounts.get(gene);}
			if(KO.totalCounts.containsKey(gene)){koTotal=KO.totalCounts.get(gene);}
			if(KO.intronCounts.containsKey(gene)){koIntron=KO.intronCounts.get(gene);}
			
			writer.write(gene.getName()+"\t"+gene.toUCSC(gene.getOrientation())+"\t"+gene.getGenomicLength()+"\t"+gene.getNumberOfBlocks()+"\t"+wtIntron+"\t"+wtTotal+"\t"+koIntron+"\t"+koTotal+"\n");
		}
		
		
		writer.close();
	}
	
}

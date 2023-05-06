package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import net.sf.samtools.util.CloseableIterator;

public class ConvertBEDGraphToTranscriptome {

	public ConvertBEDGraphToTranscriptome(File bedGraph, AnnotationCollection<Gene> genes, File save) throws IOException{
		FileWriter writer=new FileWriter(save);
		Map<SingleInterval, Double> map=BEDFileIO.loadbedgraph(bedGraph);
		
		for(SingleInterval region: map.keySet()){
			//System.err.println(region.toUCSC());
			Collection<Annotation> newRegions=convert(region, genes);
			double score=map.get(region);
			for(Annotation newRegion: newRegions){
				writer.write(newRegion.tobedgraph(score)+"\n");
			}
		}
		writer.close();
	}

	

	private Collection<Annotation> convert(SingleInterval region, AnnotationCollection<Gene> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		CloseableIterator<Gene> overlappingGenes=genes.sortedIterator(region, true);
		
		//For each overlapping gene --> convert
		while(overlappingGenes.hasNext()){
			Gene gene=overlappingGenes.next();
			Annotation newRead=convert(region, gene);
			if(newRead!=null){rtrn.add(newRead);}
		}
		return rtrn;
	}

	private Annotation convert(SingleInterval region, Gene gene) {
		Annotation newInterval1=gene.convertToFeatureSpace(region);
		return newInterval1;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File bam=new File(args[0]);
			AnnotationCollection<Gene> genes=BEDFileIO.loadFromFile(args[1]);
			String save=args[2];
			new ConvertBEDGraphToTranscriptome(bam, genes, new File(save));
		} 
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=bedgraph \n args[1]=BED file \n args[2]=save";
	
}

package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class GeneHeatmap {

	public GeneHeatmap(BarcodingData data, Collection<Gene> genes, String save) throws IOException{
		List<String> geneNames=getNames(genes);
		MatrixWithHeaders matrix=new MatrixWithHeaders(geneNames, geneNames);
		int counter=0;
		for(Gene gene1: genes){
			for(Gene gene2: genes){
				InteractionScore score=data.getInteractionScore(gene1.getGenomicRegion(), gene2.getGenomicRegion());
				double count=score.getSharedBarcodes().size();
				matrix.set(gene1.getName(), gene2.getName(), count);
			}
			counter++;
			if(counter%1000 ==0){System.err.println(counter);}
		}
		
		matrix.writeGCT(save);
	}

	private List<String> getNames(Collection<Gene> genes) {
		List<String> rtrn=new ArrayList<String>();
		for(Gene gene: genes){
			rtrn.add(gene.getName());
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		BarcodingData data=new BarcodingData(new File(args[0]));
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[1]);
		String save=args[2];
		new GeneHeatmap(data, genes, save);
	}
	
}

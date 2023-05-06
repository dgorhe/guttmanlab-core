package guttmanlab.core.sars;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

public class GetCDSUTRForGene {

	public static void main(String[] args) throws IOException{
		FileWriter writer=new FileWriter("/Users/mguttman/Downloads/proteinAnnotations.bed");
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile("/Users/mguttman/Downloads/RefSeq.hg38.proteins.bed");
		for(Gene gene: genes){
				int utr5=0;
				if(gene.get5UTR()!=null){
					utr5=gene.get5UTR().size();
				}
				int cds=gene.getCodingRegion().size();
				int utr3=0;
				if(gene.get3UTR()!=null){
					utr3=gene.get3UTR().size();
				}
				writer.write(gene.getName()+"\t"+0+"\t"+utr5+"\tUTR5\n");
				writer.write(gene.getName()+"\t"+utr5+"\t"+(utr5+cds)+"\tCDS\n");
				writer.write(gene.getName()+"\t"+(utr5+cds)+"\t"+gene.size()+"\tUTR3\n");
				
				
				
			
		}
		writer.close();
	}
	
	
}

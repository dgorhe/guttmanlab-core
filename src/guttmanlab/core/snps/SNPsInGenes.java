package guttmanlab.core.snps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;

public class SNPsInGenes {

	
	public SNPsInGenes(File snpFile, Map<String, IntervalTree<Gene>> genes, String save) throws IOException{
		
		
		Collection<SingleInterval> snps= BEDFileIO.loadSingleIntervalsFromFile(snpFile.getAbsolutePath());
		
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval snp: snps){
			if(overlapsGene(snp, genes)){
				writer.write(snp.getReferenceName()+"\t"+snp.getReferenceStartPosition()+"\t"+snp.getReferenceEndPosition()+"\t"+snp.getName()+"\n");
			}
		}
		
		writer.close();	
	}
	
	private boolean overlapsGene(SingleInterval record, Map<String, IntervalTree<Gene>> genes) {
		if(genes.containsKey(record.getReferenceName())){
			return genes.get(record.getReferenceName()).hasOverlappers(record.getReferenceStartPosition(), record.getReferenceEndPosition());
		}
		return false;
	}
	
	
	
	
	public static void main(String[] args) throws IOException{
		
		
		if(args.length>2){
			File snp=new File(args[0]);
			Map<String, IntervalTree<Gene>> genes=BEDFileIO.loadTree(args[1]);
			String save=args[2];
			
			new SNPsInGenes(snp, genes, save);
			
		}
		else{System.err.println(usage);}
		
		
	}
	
	static String usage=" args[0]=snp file \n args[1]=Gene BED file \n args[2]=save file name";	
	
}

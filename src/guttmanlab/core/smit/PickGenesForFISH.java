package guttmanlab.core.smit;

import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

public class PickGenesForFISH {
	int lengthCutoff=2000;
	int intronCutoff=1000;
	
	public PickGenesForFISH(Collection<Gene> genes){
		for(Gene gene: genes){
			if(passesCutoffs(gene)){
				System.out.println(gene.toBED());
			}
		}
		
	}
	
	private boolean passesCutoffs(Gene gene) {
		boolean cDNALength=gene.size()>lengthCutoff;
		boolean intronLength=false;
		
		for(Annotation intron: gene.getIntrons()){
			if(intron.size()>intronCutoff){intronLength=true;}
		}
		
		return cDNALength && intronLength;
	}

	public static void main(String[] args) throws IOException{
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		new PickGenesForFISH(genes);
	}
	
}

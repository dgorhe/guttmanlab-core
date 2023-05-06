package guttmanlab.core.annotationcollection;

import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

public class GeneCollection extends FeatureCollection<Gene>{

	public GeneCollection(String bedFile){
		super();
		try{
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(bedFile);
			this.addAll(genes);
		}catch(IOException ex){}	
	}	
	
	
	public static void main(String[] args){
		
	}
}
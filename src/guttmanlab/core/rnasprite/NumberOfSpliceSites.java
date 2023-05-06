package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

public class NumberOfSpliceSites {

	public NumberOfSpliceSites(String geneFile, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		Collection<Gene> genes= BEDFileIO.loadRegionsFromFile(geneFile);
		
		for(Gene g: genes){writer.write(g.getReferenceName()+"\t"+g.getReferenceStartPosition()+"\t"+g.getReferenceEndPosition()+"\ts="+g.getNumberOfBlocks()+"\n");}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		String in=args[0];
		String out=args[1];
		new NumberOfSpliceSites(in, out);
	}
	
	
}

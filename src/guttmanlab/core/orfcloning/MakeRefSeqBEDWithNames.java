package guttmanlab.core.orfcloning;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

public class MakeRefSeqBEDWithNames {

	public static void main(String[] args) throws IOException{
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		Map<String, String> idToName=BEDFileIO.loadIDToName(args[1]);
		FileWriter writer=new FileWriter(args[2]);
		
		for(Gene gene: genes){
			String name=idToName.get(gene.getName());
			gene.setName(name);
			writer.write(gene.toBED()+"\n");
		}
		writer.close();
	}
	
}

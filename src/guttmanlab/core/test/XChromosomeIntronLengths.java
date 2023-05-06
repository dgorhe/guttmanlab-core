package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class XChromosomeIntronLengths {

	public static void main(String[] args) throws IOException{
		
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		Map<String, String> names=BEDFileIO.loadIDToName(args[1]);
		FileWriter writer=new FileWriter(args[2]);
		
		Collection<SingleInterval> introns=new TreeSet<SingleInterval>();
		
		for(Gene gene: genes){
			String name=names.get(gene.getName());
			gene.setName(name);
			if(gene.getReferenceName().equalsIgnoreCase("chrX")){
				for(SingleInterval intron: gene.getIntronSet()){
					introns.add(intron);
					//writer.write(intron.toBED()+"\n");
				}
			}
		}
		
		for(SingleInterval intron: introns){
			writer.write(intron.getName()+"\t"+intron.toUCSC()+"\t"+intron.getGenomicLength()+"\n");
		}
		
		writer.close();
		
	}
	
}

package guttmanlab.core.spidr;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

public class GetGenes {

	public static void main(String[] args) throws IOException {
		GTFToJunctions gtf=new GTFToJunctions(new File(args[0]));
		List<String> list=BEDFileIO.loadLines(args[1]);
		
		Collection<Gene> genes= gtf.getGenes();
		for(Gene gene: genes) {
			if(list.contains(gene.getName())) {
				System.out.println(gene.toBED());
			}
		}
		
	}
	
}

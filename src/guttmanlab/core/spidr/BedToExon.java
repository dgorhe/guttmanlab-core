package guttmanlab.core.spidr;

import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class BedToExon {

	public static void main(String[] args) throws IOException {
		Collection<Gene> genes= BEDFileIO.loadRegionsFromFile(args[0]);
		for(Gene gene: genes) {
			
			Annotation cds=gene.getCodingRegion();
			
			if(cds!=null) {
			Iterator<SingleInterval> blocks=cds.getBlocks();
			while(blocks.hasNext()) {
				SingleInterval exon=blocks.next();
				exon.setName(gene.getName());
				System.out.println(exon.toBED());
			}
		}
		}
	}
	
}

package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class BEDToExon {

	public static void main(String[] args) throws IOException{
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		String save=args[1];
		
		Collection<SingleInterval> exons=new TreeSet<SingleInterval>();
		
		for(Gene gene: genes){
			Iterator<SingleInterval> iter=gene.getBlocks();
			while(iter.hasNext()){exons.add(iter.next());}
		}
		
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval exon: exons){writer.write(exon.toBED()+"\n");}
		
		writer.close();
		
	}
	
}

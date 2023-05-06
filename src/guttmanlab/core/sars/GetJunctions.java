package guttmanlab.core.sars;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class GetJunctions {

	public static void main(String[] args) throws IOException{
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile((args[0]));
		String save=args[1];
		
		Collection<Gene> allJunctions=new TreeSet<Gene>();
		
		for(Gene gene: genes){
			System.err.println(gene.getName());
			ArrayList<Annotation> exons=new ArrayList<Annotation>();
			exons.addAll(gene.getBlockSet());
			for(int i=1; i<exons.size(); i++){
				Annotation exon1=exons.get(i-1);
				Annotation exon2=exons.get(i);
				Collection<Annotation> list=new TreeSet<Annotation>();
				list.add(exon1);
				list.add(exon2);
				String name=gene.getName()+":junction"+i;
				Gene j=new Gene(list, name);
				allJunctions.add(j);
			}
		}
		
		FileWriter writer=new FileWriter(save);
		
		for(Gene j: allJunctions){
			writer.write(j.toBED()+"\n");
		}
		
		writer.close();
	}
	
}

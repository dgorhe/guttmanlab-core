package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class BEDToIntron {
	
	public static void main(String[] args) throws IOException{
		if(args.length>1){
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		
			String save=args[1];
			
			Collection<SingleInterval> introns=new TreeSet<SingleInterval>();
			
			int counter=0;
			for(Gene gene: genes){
				Collection<Annotation> temp=gene.getIntrons();
				for(Annotation intron: temp){introns.add(intron.getSingleInterval());}
				counter++;
				//if(counter%1000 ==0){System.err.println(counter);}
			}
			
			FileWriter writer=new FileWriter(save);
			
			for(SingleInterval intron: introns){writer.write(intron.toBED()+"\n");}
			
			writer.close();
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BED file \n args[1]=save";
}

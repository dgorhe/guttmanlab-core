package guttmanlab.core.util;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

public class IntronsToBed {

	public static void main(String[] args) throws IOException {
		if(args.length>1) {
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		FileWriter writer=new FileWriter(args[1]);
		for(Gene g: genes) {
			Collection<Annotation> introns=g.getIntrons();
			for(Annotation intron: introns) {
				writer.write(intron.toBED()+"\n");
			}
		}
		writer.close();
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=bed file \n args[1]=output";
}

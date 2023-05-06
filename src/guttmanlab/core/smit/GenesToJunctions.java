package guttmanlab.core.smit;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

public class GenesToJunctions {

	public GenesToJunctions(Collection<Gene> genes, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		for(Gene gene: genes){
			Collection<Annotation> junctions=getJunctions(gene);
			writer(writer, junctions);
		}
		writer.close();
	}

	private Collection<Annotation> getJunctions(Gene gene) {
		return gene.getExonIntronPairs();
	}

	private void writer(FileWriter writer, Collection<Annotation> junctions) throws IOException {
		for(Annotation junction: junctions){
			writer.write(junction.toBED()+"\n");
		}
		
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			new GenesToJunctions(BEDFileIO.loadRegionsFromFile(args[0]), args[1]);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Gene BED file \n args[1]=save";
}

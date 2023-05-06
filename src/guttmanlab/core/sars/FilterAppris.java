package guttmanlab.core.sars;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;

public class FilterAppris {

	public static void main(String[] args) throws IOException{
		Collection<String> tier1=parse(args[0]);
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[1]);
		String save=args[2];
		
		FileWriter writer=new FileWriter(save);
		for(Gene gene: genes){
			if(tier1.contains(gene.getName().split("\\.")[0])){
				System.err.println(gene.getName().split("\\.")[0]);
				writer.write(gene.toBED()+"\n");
			}
		}
		writer.close();
		
	}

	private static Collection<String> parse(String string) throws IOException {
		Collection<String> lines=BEDFileIO.loadLines(string);
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String line: lines){
			String name=line.split("\t")[2].split("\\.")[0];
			String princ=line.split("\t")[4];
			if(princ.equalsIgnoreCase("PRINCIPAL:1")){rtrn.add(name);}
		}
		
		return rtrn;
	}
	
}

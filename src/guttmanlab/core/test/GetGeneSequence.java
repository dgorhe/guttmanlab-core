package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class GetGeneSequence {

	public static void main(String[] args) throws IOException{
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		Map<String, Sequence> sequenceByChr=FastaFileIOImpl.readFromFileByName(args[1]);
		Collection<String> names=new TreeSet<String>();
		
		FileWriter writer=new FileWriter(args[2]+".mRNA.fa");
		
		for(Gene gene: genes){
			if(sequenceByChr.containsKey(gene.getReferenceName())){
				Sequence chrSeq=sequenceByChr.get(gene.getReferenceName());
				Sequence geneSeq=gene.getSequence(chrSeq);
				writer.write(">"+gene.getName()+"\n"+geneSeq.toString()+"\n");
				if(names.contains(gene.getName())){System.err.println("Duplicate "+gene.toBED());}
				names.add(gene.getName());
			}
			else{System.err.println(gene.toUCSC());}
		}
		
		writer.close();
		
		
		writer=new FileWriter(args[2]+".mRNA.bed");
		
		for(Gene gene: genes){
			int UTR1=0;
			if(gene.get5UTR()!=null){
				UTR1=gene.get5UTR().size();
			}
			int CDS=gene.getCodingRegion().size();
			
			int UTR2=0;
			if(gene.get3UTR()!=null){
				UTR2=gene.get3UTR().size();
			}
			
			if(UTR1>0){writer.write(gene.getName()+"\t"+0+"\t"+UTR1+"\t5'UTR"+"\n");}
			writer.write(gene.getName()+"\t"+UTR1+"\t"+(UTR1+CDS)+"\tCDS"+"\n");
			if(UTR2>0){writer.write(gene.getName()+"\t"+(UTR1+CDS)+"\t"+(UTR1+CDS+UTR2)+"\t3'UTR"+"\n");}
			
			/*Sequence chrSeq=sequenceByChr.get(gene.getReferenceName());
			Sequence geneSeq=gene.getSequence(chrSeq);
			writer.write(">"+gene.getName()+"\n"+geneSeq.toString()+"\n");*/
		}
		
		writer.close();
	
		
		writer=new FileWriter(args[2]+".intron.fa");
		
		Collection<SingleInterval> introns=new TreeSet<SingleInterval>();
		for(Gene gene: genes){
			introns.addAll(gene.getIntronSet());
		}
		
		System.err.println("Num introns "+introns.size());
		
		int counter=0;
		for(SingleInterval intron: introns){
			if(sequenceByChr.containsKey(intron.getReferenceName())){
				Sequence chrSeq=sequenceByChr.get(intron.getReferenceName());
				String name=intron.getName()+"_"+intron.toUCSC();
				//System.err.println(name);
				Sequence geneSeq=chrSeq.getSubsequence(intron);
				writer.write(">"+name+"\n"+geneSeq.toString()+"\n");
			}
			counter++;
			if(counter % 1000 ==0){System.err.println(counter);}
		}
		
		writer.close();
		
		
	}
	
}

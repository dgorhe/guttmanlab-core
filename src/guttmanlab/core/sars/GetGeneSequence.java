package guttmanlab.core.sars;

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
		
		FileWriter writer=new FileWriter(args[2]);
		for(Gene gene: genes){
			double val=Math.random();
			if(val<.1){
				SingleInterval region=gene.getGenomicRegion();
				String seq=sequenceByChr.get(region.getReferenceName()).getSubsequence(region).getSequenceBases();
				String name=gene.getName();
				writer.write(">"+name+"\n"+seq+"\n");
			}
		}
		
		writer.close();
		
		/*FileWriter writer=new FileWriter(args[2]+".intron.fa");
		
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
				if(names.contains(name)){System.err.println("Duplicate "+intron.toBED());}
				names.add(name);
			}
			counter++;
			if(counter % 1000 ==0){System.err.println(counter);}
		}
		
		writer.close();*/
		
		
	}
	
}

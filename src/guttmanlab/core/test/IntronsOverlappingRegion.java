package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;

public class IntronsOverlappingRegion {

	public static void main(String[] args) throws IOException{
		Map<String, IntervalTree<Gene>> genes=BEDFileIO.loadTree(args[0]);
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[1]);
		
		FileWriter writer=new FileWriter(args[2]+".intron.bed");
		FileWriter exon=new FileWriter(args[2]+".exon.bed");
		FileWriter geneWriter=new FileWriter(args[2]+".gene.bed");
		
		for(SingleInterval region: regions){
			Iterator<Node<Gene>> iter=genes.get(region.getReferenceName()).overlappers(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			while(iter.hasNext()){
				Gene g=iter.next().getValue();
				geneWriter.write(g.toBED()+"\n");
				for(SingleInterval intron: g.getIntronSet()){writer.write(intron.toBED()+"\n");}
				for(Annotation e: g.getBlockSet()){exon.write(e.toBED()+"\n");}
			}
		}
		
		geneWriter.close();
		writer.close();
		exon.close();
	}
	
}

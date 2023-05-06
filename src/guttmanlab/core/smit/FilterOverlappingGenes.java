package guttmanlab.core.smit;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.datastructures.IntervalTree;
import net.sf.samtools.util.CloseableIterator;

public class FilterOverlappingGenes {

	public FilterOverlappingGenes(AnnotationCollection<Gene> genes, String save) throws IOException{
		Collection<Gene> uniqueGenes=makeTree(genes); //TODO Make positive and negative tree separately
		System.err.println("Made gene tree");
		
		write(save, uniqueGenes);
	}
	
	private void write(String save, Collection<Gene> uniqueGenes) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: uniqueGenes){
			writer.write(gene.toBED()+"\n");
		}
		
		writer.close();
	}

	private Collection<Gene> makeTree(AnnotationCollection<Gene> genes) {
		Map<String, IntervalTree<Gene>> tree=new TreeMap<String, IntervalTree<Gene>>();
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		CloseableIterator<Gene>iter=genes.sortedIterator();
		while(iter.hasNext()){
			Gene gene=iter.next();
			IntervalTree<Gene> temp=new IntervalTree<Gene>();
			if(tree.containsKey(gene.getReferenceName())){
				temp=tree.get(gene.getReferenceName());
			}
			if(!temp.hasOverlappers(gene.getReferenceStartPosition(), gene.getReferenceEndPosition())){
				rtrn.add(gene);
				
			}
			else{
				//TODO Take longer of the 2 isoforms
				/*Gene gene1=temp.overlappers(gene.getReferenceStartPosition(), gene.getReferenceEndPosition()).next().getValue();
				if(gene1.getGenomicLength()<gene.getGenomicLength()){
					//Remove gene1 add gene
					rtrn.remove(gene1);
					rtrn.add(gene);
				}*/
			}
			temp.put(gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene);
			tree.put(gene.getReferenceName(), temp);
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
		AnnotationCollection<Gene> genes=BEDFileIO.loadFromFile((args[0]));
		String save=args[1];
		new FilterOverlappingGenes(genes, save);}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=BED file \n args[1]=save";
}

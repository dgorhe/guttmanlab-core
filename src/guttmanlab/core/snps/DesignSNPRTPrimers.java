package guttmanlab.core.snps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;

public class DesignSNPRTPrimers {

	public DesignSNPRTPrimers(Collection<Gene> genes, Map<String, IntervalTree<String>> snpTree, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		
		for(Gene gene: genes) {
			System.err.println(gene.getName());		
			Collection<SingleInterval> regions=getSNPs(gene, snpTree);
			write(writer, gene, regions);
		}
		
		writer.close();
	}

	private Collection<SingleInterval> getSNPs(Gene gene, Map<String, IntervalTree<String>> snpTree) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		Iterator<SingleInterval> iter=gene.getBlocks(); 
		
		while(iter.hasNext()){
			SingleInterval exon=iter.next();
			rtrn.addAll(getSNPs(exon, snpTree));
		}
		
		return rtrn;
	}

	private Collection<SingleInterval> getSNPs(SingleInterval exon, Map<String, IntervalTree<String>> snpTree) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		if(snpTree.containsKey(exon.getReferenceName())) {
			IntervalTree<String> tree=snpTree.get(exon.getReferenceName());
			Iterator<Node<String>> iter=tree.overlappers(exon.getReferenceStartPosition(), exon.getReferenceEndPosition());
			while(iter.hasNext()) {
				Node<String> next=iter.next();
				SingleInterval region=new SingleInterval(exon.getReferenceName(), next.getStart(), next.getEnd());
				region.setName(next.getValue());
				rtrn.add(region);}
		}
		
		return rtrn;
	}

	private void write(FileWriter writer, Gene gene, Collection<SingleInterval> regions) throws IOException {
		for(SingleInterval region: regions) {writer.write(gene.getName()+"\t"+region.toUCSC()+"\t"+region.getName()+"\n");}
	}
	
	public static void main(String[] args) throws IOException {
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		Collection<String> chrs=new TreeSet<String>();
		for(Gene gene: genes) {chrs.add(gene.getReferenceName());}
		
		Map<String, IntervalTree<String>> snpTree=BEDFileIO.parseSNPTree(new File(args[1]), chrs);
		String save=args[2];
		new DesignSNPRTPrimers(genes, snpTree, save);
	}
	
}

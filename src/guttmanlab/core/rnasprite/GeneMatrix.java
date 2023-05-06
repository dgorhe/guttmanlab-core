package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.barcoding.analysis.BarcodingDataStreaming;
import guttmanlab.core.barcoding.analysis.Cluster;
import guttmanlab.core.datastructures.IntervalTree;

public class GeneMatrix {

	public GeneMatrix(BarcodingDataStreaming data, Gene gene, Map<String, IntervalTree<Gene>> genes, String save) throws IOException{
		Map<Gene, Collection<String>> geneToBarcode=geneToBarcode(data, genes);
		
		write(save, gene, geneToBarcode);
		
		
	}
	
	private void write(String save, Gene gene, Map<Gene, Collection<String>> geneToBarcode) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene2: geneToBarcode.keySet()){
			Collection<String> barcodes1=geneToBarcode.get(gene);
			Collection<String> barcodes2=geneToBarcode.get(gene2);
			System.err.println(barcodes1.size()+" "+barcodes2.size());
			int count=intersect(barcodes1, barcodes2);
			if(count>1){
				writer.write(gene2.getName()+"\t"+count+"\t"+barcodes2.size()+"\n");
			}
		}
		
		
		/*
		//writer.write("Name");
		for(Gene gene: geneToBarcode.keySet()){
			if(geneToBarcode.get(gene).size()>10){
				System.err.println(gene.getName()+" "+geneToBarcode.get(gene).size());
				writer.write(gene.getName()+"\t");
			}
		}
		writer.write("\n");
		
		for(Gene gene1: geneToBarcode.keySet()){
			if(geneToBarcode.get(gene1).size()>10){
				writer.write(gene1.getName());
				System.err.println(gene1.getName());
				for(Gene gene2: geneToBarcode.keySet()){
					if(geneToBarcode.get(gene2).size()>10){
					Collection<String> barcodes1=geneToBarcode.get(gene1);
					Collection<String> barcodes2=geneToBarcode.get(gene2);
					int count=intersect(barcodes1, barcodes2);
					
					writer.write("\t"+count);
					}
				}
				writer.write("\n");
			}
		}*/
		
		writer.close();
	}

	private int intersect(Collection<String> barcodes1, Collection<String> barcodes2) {
		int counter=0;
		for(String barcode: barcodes1){
			if(barcodes2.contains(barcode)){counter++;}
		}
		return counter;
	}

	private Map<Gene, Collection<String>> geneToBarcode(BarcodingDataStreaming data, Map<String, IntervalTree<Gene>> genes) {
		Map<Gene, Collection<String>> geneToBarcode=new TreeMap<Gene, Collection<String>>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			for(SingleInterval interval: c.getAllIntervals()){
				Collection<Gene> geneList=toGene(interval, genes);
				if(!geneList.isEmpty()){
					for(Gene gene: geneList){
						Collection<String> barcodes=new TreeSet<String>();
						if(geneToBarcode.containsKey(gene)){barcodes=geneToBarcode.get(gene);}
						barcodes.add(c.getBarcode());
						geneToBarcode.put(gene, barcodes);
					}
				}
			}
		}
		return geneToBarcode;
	}

	private Collection<Gene> toGene(SingleInterval interval, Map<String, IntervalTree<Gene>> genes) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		if(genes.containsKey(interval.getReferenceName())){
			Iterator<Gene> iter=genes.get(interval.getReferenceName()).overlappingValueIterator(interval.getReferenceStartPosition(), interval.getReferenceEndPosition());
			if(iter.hasNext()){
				rtrn.add(iter.next());
			}
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Map<String, IntervalTree<Gene>> tree=BEDFileIO.loadTree(args[1]);
			String geneName=args[2];
			String save=args[3];
			Gene gene=getGene(tree, geneName);
			if(gene==null){throw new IllegalArgumentException("No gene "+geneName);}
			else{System.err.println(gene.toBED());}
			new GeneMatrix(data, gene, tree, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=barcode file \n args[1]=BED file \n args[2]=gene name \n args[3]=save";

	private static Gene getGene(Map<String, IntervalTree<Gene>> tree, String geneName) {
		for(String chr: tree.keySet()){
			Iterator<Gene> iter=tree.get(chr).valueIterator();
			while(iter.hasNext()){
				Gene gene=iter.next();
				if(gene.getName().equalsIgnoreCase(geneName)){return gene;}
			}
		}
		return null;
	}
	
}

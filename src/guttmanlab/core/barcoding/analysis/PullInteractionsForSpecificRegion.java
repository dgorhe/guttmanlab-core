package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;

public class PullInteractionsForSpecificRegion {

	//Given a region, get all interacting clusters
	public PullInteractionsForSpecificRegion(Annotation region, File data, String save, int minClusterSize, int maxClusterSize, boolean excludeQueryRegion) throws NumberFormatException, IOException{
		Map<String, IntervalTree<String>> trees=BarcodeFileParser.parseBarcodeFile(data, region, minClusterSize, maxClusterSize);
		
		FileWriter writer=new FileWriter(save);
		
		for(String chr: trees.keySet()){
			IntervalTree<String> tree=trees.get(chr);
			Iterator<Node<String>> iter=tree.iterator();
			while(iter.hasNext()){
				Node<String> next=iter.next();
				SingleInterval interval=new SingleInterval(chr, next.getStart(), next.getEnd());
				if(!excludeQueryRegion || !interval.overlaps(region)){
					Collection<String> barcodes=next.getContainedValues();
					for(String barcode: barcodes){
						writer.write(chr+"\t"+next.getStart()+"\t"+next.getEnd()+"\t"+barcode+"\n");
					}
				}
			}
		}
		writer.close();
	}
	
	public PullInteractionsForSpecificRegion(AnnotationCollection<Gene> regions, File data, String save, int minClusterSize, int maxClusterSize) throws NumberFormatException, IOException{
		
		Map<String, IntervalTree<String>> trees=BarcodeFileParser.parseBarcodeFile(data, regions, minClusterSize, maxClusterSize);
		
		FileWriter writer=new FileWriter(save);
		
		for(String chr: trees.keySet()){
			IntervalTree<String> tree=trees.get(chr);
			Iterator<Node<String>> iter=tree.iterator();
			while(iter.hasNext()){
				Node<String> next=iter.next();
				SingleInterval interval=new SingleInterval(chr, next.getStart(), next.getEnd());
				Collection<String> barcodes=next.getContainedValues();
				for(String barcode: barcodes){
					writer.write(chr+"\t"+next.getStart()+"\t"+next.getEnd()+"\t"+barcode+"\n");
				}
			}
		}
		writer.close();
	}


	public static void main(String[] args)throws IOException{
		if(args.length>4){
		File file=new File(args[0]);
		String save=args[3];
		
		String bedFile=args[4];
		AnnotationCollection<Gene> regions=BEDFileIO.loadFromFile(bedFile);
		
		int minClusterSize=new Integer(args[1]);
		int maxClusterSize=new Integer(args[2]);
		
		new PullInteractionsForSpecificRegion(regions, file, save, minClusterSize, maxClusterSize);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Barcode file \n args[1]=min cluster size \n args[2]=max cluster size \n args[3]=save \n args[4]=regions (BED file, will look for cluster with ALL regions)";
	
}

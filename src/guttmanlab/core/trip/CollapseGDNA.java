package guttmanlab.core.trip;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.sequence.Sequence;

public class CollapseGDNA {

	//static int distance=1;
	
	

	private static Map<String, Collection<String>> collapseByBarcode(List<String> lines) {
		Map<String, Collection<String>> barcodes=new TreeMap<String, Collection<String>>();
		
		int counter=0;
		for(String line: lines) {
			String[] tokens=line.split(" ");
			String barcode=tokens[0];
			String pos=tokens[1]+" "+tokens[2];
			if(!barcodes.containsKey(barcode)) {barcodes.put(barcode, new TreeSet<String>());}
			Collection<String> set=barcodes.get(barcode);
			set.add(pos);
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+lines.size());}
		}
		return barcodes;
	}
	
	private static void write(String save, Map<String, Collection<String>> map, Map<String, IntervalTree<Double>> genes) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: map.keySet()) {
			Collection<String> list=map.get(barcode);
			
			Collection<SingleInterval> positions=new TreeSet<SingleInterval>();
			Collection<Double> overlappingGenes=new TreeSet<Double>();
			
			for(String pos: list) {
				String chr=pos.split(" ")[0];
				int start=Integer.parseInt(pos.split(" ")[1]);
				SingleInterval region=new SingleInterval(chr, start, start+1);
				SingleInterval binned=region.bin(1000);
				positions.add(binned);
				
				overlappingGenes.addAll(getOverlappingGenes(chr, start, genes));
			}
			
			if(positions.size()==1) {
				double maxVal=max(overlappingGenes);
				writer.write(barcode+"\t"+list.size()+"\t"+positions.iterator().next().toUCSC()+"\t"+overlappingGenes.size()+"\t"+maxVal+"\n");
			}
		}
		writer.close();
	}
	
	private static double max(Collection<Double> overlappingGenes) {
		if(overlappingGenes.isEmpty()) {return -1;}
		return Statistics.max(overlappingGenes);
	}

	private static Collection<Double> getOverlappingGenes(String chr, int start, Map<String, IntervalTree<Double>> genes) {
		Collection<Double> rtrn=new TreeSet<Double>();
		if(genes.containsKey(chr)) {
			IntervalTree<Double> tree=genes.get(chr);
			Iterator<Double> iter=tree.overlappingValueIterator(start, start+1);
			while(iter.hasNext()) {
				double val=iter.next();
				rtrn.add(val);
			}
		}
		return rtrn;
	}

	private static Map<String, Collection<String>> findClosest(Map<String, Collection<String>> map, int distance) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		for(String barcode1: map.keySet()) {
			Collection<String> list=new TreeSet<String>();
			for(String barcode2: map.keySet()) {
				if(barcode1.length()==barcode2.length()) {
					int dist=distance(barcode1, barcode2);
					if(dist<=distance) {list.add(barcode2);}
				}
			}
			rtrn.put(barcode1, list);
		}
		return rtrn;
	}
	
	
	public static Collection<String> getCollapsedSets(Map<String, Collection<String>> barcodes, int distance){
		Map<String, Collection<String>> nearest=findClosest(barcodes, distance);
		
		Collection<String> sets=makeSets(nearest);
		return sets;
	}
	
	private static Collection<String> makeSets(Map<String, Collection<String>> nearest) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String barcode: nearest.keySet()) {
			System.err.println(barcode+" "+nearest.get(barcode).size());
			
			Collection<String> visited=new TreeSet<String>();
			Collection<String> set=new TreeSet<String>();
			set.addAll(nearest.get(barcode));
			
			//TODO while difference in set and visited is not empty, search for the remainder
			Collection<String> diff= difference(set, visited);
			
			while(diff.size()>0) {
				for(String barcode2: diff) {
					set.addAll(nearest.get(barcode2));
					visited.add(barcode2);
				}
				diff=difference(set, visited);
			}
			
			String setString=toString(set);
			rtrn.add(setString);
		}
			
			
		return rtrn;
	}
	
	private static String toString(Collection<String> set) {
		String rtrn="";
		for(String s: set) {
			rtrn+=s+" ";
		}
		return rtrn;
	}
	
	private static Collection<String> difference(Collection<String> set, Collection<String> visited) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String val: set) {
			if(!visited.contains(val)) {rtrn.add(val);}
		}
		
		return rtrn;
	}
	
	private static int distance(String barcode1, String barcode2) {
		return Sequence.distance(barcode1, barcode2);
	}
	
	
	private static Map<String, Collection<String>> getPositions(Collection<String> sets, Map<String, Collection<String>> map) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String set: sets) {
			String[] tokens=set.split(" ");
			Collection<String> list=new TreeSet<String>();
			for(int i=0; i<tokens.length; i++) {
				list.addAll(map.get(tokens[i]));
			}
			rtrn.put(set, list);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		if(args.length>3) {
			
			List<String> lines=BEDFileIO.loadLines(args[0]);
			Map<String, Collection<String>> map=collapseByBarcode(lines);
			Collection<String> sets=getCollapsedSets(map, Integer.parseInt(args[1]));
			Map<String, Collection<String>> positions=getPositions(sets, map);
			Map<String, IntervalTree<Double>> genes=BEDFileIO.loadTreePlusExpression(args[2]);
			write(args[3], positions, genes);
			
		}
		else {System.err.println(usage);}
	}

	static String usage=" args[0]=original gDNA file \n args[1]=hamming distance to collapse \n args[2]=gene expression \n args[3]=save";

	

	

	
	
}

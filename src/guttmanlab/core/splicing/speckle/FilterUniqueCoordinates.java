package guttmanlab.core.splicing.speckle;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;

public class FilterUniqueCoordinates {

	public static void filter(Map<String, SingleInterval> regions, Map<String, String> geneToFullLine, String out) throws IOException {
		Map<String, IntervalTree<String>> trees=makeTree(regions);
		
		FileWriter writerUnique=new FileWriter(out+".unique");
		FileWriter writerNonunique=new FileWriter(out+".nonunique");
		
		for(String gene: regions.keySet()) {
			if(overlaps(trees, regions.get(gene), gene)) {
				writerNonunique.write(geneToFullLine.get(gene)+"\n");
			}
			else {writerUnique.write(geneToFullLine.get(gene)+"\n");}
		}
		
		writerUnique.close();
		writerNonunique.close();
	}

	private static Map<String, IntervalTree<String>> makeTree(Map<String, SingleInterval> regions) {
		Map<String, IntervalTree<String>> rtrn=new TreeMap<String, IntervalTree<String>>();
		for(String gene: regions.keySet()) {
			SingleInterval region=regions.get(gene);
			if(!rtrn.containsKey(region.getReferenceName())) {
				rtrn.put(region.getReferenceName(), new IntervalTree<String>());
			}
			IntervalTree<String> tree=rtrn.get(region.getReferenceName());
			tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), gene);
		}
		return rtrn;
	}

	private static boolean overlaps(Map<String, IntervalTree<String>> trees, SingleInterval region, String gene) {
		if(trees.containsKey(region.getReferenceName())) {
			IntervalTree<String> tree=trees.get(region.getReferenceName());
			Iterator<String> iter=tree.overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			while(iter.hasNext()) {
				String r2=iter.next();
				if(!gene.equals(r2)) {return true;}
			}
		}
		return false;
	}

	
	public static void main(String[] args) throws IOException {
		Map<String, SingleInterval> regions=parseCoordinates(args[0]);
		Map<String, String> lines=parseFullLine(args[0]);
		String save=args[1];
		filter(regions, lines, save);
	}

	private static Map<String, SingleInterval> parseCoordinates(String file) throws IOException {
		Map<String, SingleInterval> rtrn=new TreeMap<String, SingleInterval>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			rtrn.put(tokens[0], new SingleInterval(tokens[1]));
		}
		reader.close();
		return rtrn;
	}

	private static Map<String, String> parseFullLine(String file) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			rtrn.put(tokens[0], nextLine);
		}
		reader.close();
		return rtrn;
	}
	
	
}

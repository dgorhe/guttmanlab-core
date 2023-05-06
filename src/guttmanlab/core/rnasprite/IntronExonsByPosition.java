package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.Pair;

public class IntronExonsByPosition {

	public static void getIntronExonsByPosition(BarcodingDataStreaming data, int resolution) {
		Map<String, Integer> geneCount=new TreeMap<String, Integer>();
		int counter=0;
		double total=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Cluster binned=c.bin(resolution);
			Collection<RNAInterval> rnas= c.getAllRNARegions();
			for(RNAInterval rna: rnas) {
				if(overlaps(rna, binned)) {
					Collection<String> names=rna.getGeneNames();
					for(String name: names) {
						int count=0;
						if(geneCount.containsKey(name)) {count=geneCount.get(name);}
						count++;
						geneCount.put(name, count);
					}
					
				//	if(rna.isExon() || rna.isIntron()) {
						//System.out.println(rna.getReferenceName()+"\t"+rna.getReferenceStartPosition()+"\t"+rna.getReferenceEndPosition()+"\t"+rna.getType());
					//}	
					total++;
				}
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		data.close();
		
		write(geneCount, total);
		
	}

	private static void write(Map<String, Integer> geneCount, double total) {
		for(String gene: geneCount.keySet()) {
			double count=geneCount.get(gene);
			double ratio=(count/total)*1000000;
			System.out.println(gene+"\t"+count+"\t"+total+"\t"+ratio);
		}
		
	}

	private static boolean overlaps(RNAInterval rna, Cluster binned) {
		for(SingleInterval region: binned.getAllDNAIntervals()) {
			if(rna.overlaps(region)) {
				//System.err.println(rna.toUCSC()+" "+region.toUCSC());
				return true;
			}
		}
		return false;
	}

	/*private static void write(Map<SingleInterval, Pair<Integer>> scores) {
		for(SingleInterval region: scores.keySet()) {
			Pair<Integer> score=scores.get(region);
			System.out.println(region.toBedgraph(score.getValue1())+"\t"+score.getValue2());
		}
		
	}*/
	
	public static void main(String[] args) throws IOException, NumberFormatException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		getIntronExonsByPosition(data, Integer.parseInt(args[1]));
	}
	
}

package guttmanlab.core.xist;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;

public class DistanceToLocus {

	public static void main (String[] args) throws IOException {
		System.err.println("rebuild");
		
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		
		int binResolution=Integer.parseInt(args[1]);
		
		int count=0;
		
		
		Map<SingleInterval, Double> scores=new TreeMap<SingleInterval, Double>();
		while(data.hasNext()) {
			Cluster c=data.next();
			Cluster binned=c.bin(binResolution);
			if(binned.getClusterSize()>=2 && binned.getClusterSize()<=100) {
				for(SingleInterval r: binned.getAllDNAIntervals()){
						double score=0;
						if(scores.containsKey(r)) {score=scores.get(r);}
						score+=(1.0/binned.getClusterSize());
						scores.put(r, score);
					}
					
					
					
				}
				count++;
			}
		
		System.err.println(count);
		BEDFileIO.writeBEDGraph(scores, args[2]);
		smooth(args[2]+".smoothed.bedgraph",scores, binResolution);
	}

	private static void smooth(String save, Map<SingleInterval, Double> scores, int binResolution) throws IOException {
		FileWriter writer=new FileWriter(save);
		Map<SingleInterval, Collection<SingleInterval>> bins=new TreeMap<SingleInterval, Collection<SingleInterval>>();
		
		for(SingleInterval r: scores.keySet()) {
			int start=r.getReferenceStartPosition()-(4*binResolution);
			int end=r.getReferenceEndPosition()+(4*binResolution);
			SingleInterval newInterval=new SingleInterval(r.getReferenceName(), start, end);
			Collection<SingleInterval> windows= newInterval.getWindowsCollection(binResolution, binResolution);
			bins.put(r, windows);
		}
		
		
		
		
		for(SingleInterval r: scores.keySet()) {
			Collection<SingleInterval> windows=bins.get(r);
			List<Double> list=new ArrayList<Double>();
			for(SingleInterval window: windows) {
				if(scores.containsKey(window)) {list.add(scores.get(window));}
			}
			writer.write(r.toBedgraph(Statistics.mean(list))+"\n");
		}
		writer.close();
	}

	private static Map<String, IntervalTree<SingleInterval>> make(String chr, int binResolution) {
		// TODO Auto-generated method stub
		return null;
	}

	private static Iterator<SingleInterval> getOverlappers(Map<String, IntervalTree<SingleInterval>> tree,SingleInterval region) {
		if(tree.containsKey(region.getReferenceName())) {
			return tree.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
		}
		return new ArrayList<SingleInterval>().iterator();
	}
	
}

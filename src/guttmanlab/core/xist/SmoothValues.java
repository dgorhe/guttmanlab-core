package guttmanlab.core.xist;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.math.Statistics;

public class SmoothValues {

	
	public static void main(String[] args) throws IOException {
		File bedgraph=new File(args[0]);
		String save=args[1];
		
		IntervalTree<Double> tree=getPositions(bedgraph);
		
		int start=tree.min().getStart();
		int end=tree.max().getEnd();
		
		Collection<SingleInterval> windows=new TreeSet<SingleInterval>();
		
		for(int i=start; i<end; i+=2500) {
			SingleInterval region=new SingleInterval("chrX",i, i+25000);
			windows.add(region);
		}
		
		Map<SingleInterval, Double> smoothed=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval window: windows) {
			double avg=getAverage(window, tree);
			SingleInterval center=new SingleInterval(window.getReferenceName(),window.getMidPoint().getReferenceStartPosition()-1250, window.getMidPoint().getReferenceStartPosition()+1250);
			smoothed.put(center, avg);
		}
		
		
		BEDFileIO.writeBEDGraph(smoothed, save);
	}

	private static double getAverage(SingleInterval window, IntervalTree<Double> tree) {
		double sum=0;
		double count=0;
		Iterator<Double> iter=tree.overlappingValueIterator(window.getReferenceStartPosition(), window.getReferenceEndPosition());
		while(iter.hasNext()) {
			sum+=iter.next();
			count++;
		}
		
		if(count==0) {
			double val1=tree.max(window.getReferenceStartPosition(), window.getReferenceEndPosition()).getValue();
			double val2=tree.min(window.getReferenceStartPosition(), window.getReferenceEndPosition()).getValue();
			return Math.max(val1, val2);
		}
		
		return sum/count;
	}

	private static Collection<SingleInterval> getMissingBins(IntervalTree<Double> tree) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		Iterator<Node<Double>> iter=tree.iterator();
		
		iter.hasNext();
		Node<Double> current=iter.next();
		while(iter.hasNext()) {
			Node<Double> next=iter.next();
			if(next.getStart()!=current.getEnd()) {
				SingleInterval missing=new SingleInterval("chrX", current.getEnd()+1, next.getStart()-1);
				rtrn.add(missing);
			}
			current=next;
		}
		return rtrn;
	}

	private static Map<SingleInterval, Double> interpolateMissing(Collection<SingleInterval> missing, IntervalTree<Double> tree) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval interval: missing) {
			double minScore=tree.min(interval.getReferenceStartPosition(), interval.getReferenceEndPosition()).getValue();
			double maxScore=tree.max(interval.getReferenceStartPosition(), interval.getReferenceEndPosition()).getValue();
			double avg=(minScore+maxScore)/2.0;
			rtrn.put(interval, avg);
		}
		
		return rtrn;
	}

	private static IntervalTree<Double> getPositions(File bedgraph) throws IOException {
		IntervalTree<Double> rtrn=new IntervalTree<Double>();
		TreeMap<SingleInterval, Double> map=BEDFileIO.loadbedgraph(bedgraph);
		
		for(SingleInterval interval: map.keySet()) {
			if(interval.getReferenceName().equals("chrX")) {
				rtrn.put(interval.getReferenceStartPosition(), interval.getReferenceEndPosition(), map.get(interval));
			}
		}
		return rtrn;
	}
	
}

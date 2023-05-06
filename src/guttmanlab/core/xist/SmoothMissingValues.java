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

public class SmoothMissingValues {

	static int binSize=5;
	public static void main(String[] args) throws IOException {
		File bedgraph=new File(args[0]);
		String save=args[1];
		
		IntervalTree<Double> tree=getPositions(bedgraph);
		Collection<SingleInterval> missing=getMissingBins(tree);
		Map<SingleInterval, Double> missingValues=interpolateMissing(missing, tree);
		BEDFileIO.writeBEDGraph(missingValues, save);
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
			//double minScore=tree.min(interval.getReferenceStartPosition(), interval.getReferenceEndPosition()).getValue();
			double avg=get(tree, interval, binSize);
			//double maxScore=tree.max(interval.getReferenceStartPosition(), interval.getReferenceEndPosition()).getValue();
			//double avg=(minScore+maxScore)/2.0;
			rtrn.put(interval, avg);
		}
		
		return rtrn;
	}

	private static double get(IntervalTree<Double> tree, SingleInterval interval, int binSize2) {
		Iterator<Node<Double>> previous=tree.reverseIterator(interval.getReferenceStartPosition(), interval.getReferenceEndPosition());
		
		List<Double> vals=new ArrayList<Double>();
		
		int counter=0;
		while(previous.hasNext() && counter<binSize2) {
			vals.add(previous.next().getValue());
			counter++;
		}
		
		Iterator<Double> next=tree.getNodesAfterInterval(interval.getReferenceStartPosition(), interval.getReferenceEndPosition());
		
		counter=0;
		while(next.hasNext() && counter<binSize2) {
			vals.add(next.next());
			counter++;
		}
		return Statistics.mean(vals);
		
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

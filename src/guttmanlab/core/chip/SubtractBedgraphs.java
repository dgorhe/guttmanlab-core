package guttmanlab.core.chip;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class SubtractBedgraphs {

	static int min=20;
	
	public static Map<SingleInterval, Double> subtract(TreeMap<SingleInterval, Double> map1, TreeMap<SingleInterval, Double> map2){
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		Collection<SingleInterval> allIntervals=new TreeSet<SingleInterval>();
		allIntervals.addAll(map1.keySet());
		allIntervals.addAll(map2.keySet());
		
		for(SingleInterval region: allIntervals) {
			double val1=get(map1, region);
			double val2=get(map2, region);
			if(val1>min || val2>min) {
				double diff=Math.log(val1/val2)/Math.log(2);
				rtrn.put(region, diff);
			}
		}
		
		return rtrn;
	}	
	
	
	private static double get(TreeMap<SingleInterval, Double> map1, SingleInterval region) {
		double rtrn=0;
		if(map1.containsKey(region)) {rtrn=map1.get(region);}
		return rtrn;
	}


	public static void main(String[] args) throws IOException {
		TreeMap<SingleInterval, Double> map1= BEDFileIO.loadbedgraph(new File(args[0]));
		TreeMap<SingleInterval, Double> map2= BEDFileIO.loadbedgraph(new File(args[1]));
		String save=args[2];
		Map<SingleInterval, Double> diff=subtract(map1, map2);
		BEDFileIO.writeBEDGraph(diff, save);
	}
	
}

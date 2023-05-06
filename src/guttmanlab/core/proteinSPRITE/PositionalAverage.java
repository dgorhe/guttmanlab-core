package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.math.Statistics;

public class PositionalAverage {

	
	public static void main(String[] args) throws IOException {
		Map<SingleInterval, Double> map=BEDFileIO.loadbedgraph(new File(args[0]));
		FileWriter writer=new FileWriter(args[1]);
		
		Map<Integer, List<Double>> set=new TreeMap<Integer, List<Double>>(); 
		
		for(SingleInterval r: map.keySet()) {
			List<Double> vals=new ArrayList<Double>();
			if(set.containsKey(r.getReferenceStartPosition())) {vals=set.get(r.getReferenceStartPosition());}
			vals.add(map.get(r));
			set.put(r.getReferenceStartPosition(), vals);
		}
		
		for(Integer position: set.keySet()) {
			List<Double>vals=set.get(position);
			double average=Statistics.quantile(vals, 0.9);
			writer.write(position+"\t"+average+"\n");
		}
		
		writer.close();
		
	}
	
}

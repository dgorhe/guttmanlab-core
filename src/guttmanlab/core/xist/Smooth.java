package guttmanlab.core.xist;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.math.Statistics;

public class Smooth {
	
	public static void main(String[] args) throws IOException {
		Map<SingleInterval, Double> map=BEDFileIO.loadbedgraph(new File(args[0]));
		String save=args[1];
		smooth(save, map);
	}

	
	private static void smooth(String save, Map<SingleInterval, Double> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		Map<SingleInterval, Collection<SingleInterval>> bins=new TreeMap<SingleInterval, Collection<SingleInterval>>();
		
		for(SingleInterval r: scores.keySet()) {
			int binResolution=r.size();
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
	
}

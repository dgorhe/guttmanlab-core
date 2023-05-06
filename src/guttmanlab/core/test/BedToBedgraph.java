package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.math.Statistics;

public class BedToBedgraph {

	public BedToBedgraph(Collection<SingleInterval> regions, String save, double expected) throws IOException {
		Map<Integer, Double> counts=new TreeMap<Integer, Double>();
		
		
		for(SingleInterval region: regions) {
			for(int i=region.getReferenceStartPosition(); i<region.getReferenceEndPosition(); i++) {
				double count=0;
				if(counts.containsKey(i)) {count=counts.get(i);}
				count++;
				counts.put(i, count);
			}
		}
		
		counts=norm(counts, expected);
		
		write(save, counts);
	}

	private Map<Integer, Double> norm(Map<Integer, Double> counts, double expected) {
		Map<Integer, Double> rtrn=new TreeMap<Integer, Double>();
		
		for(int pos: counts.keySet()) {
			double norm=counts.get(pos)/expected;
			rtrn.put(pos, norm);
		}
		
		return rtrn;
	}

	private void write(String save, Map<Integer, Double> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		double avg=Statistics.mean(counts.values());
		
		for(int position: counts.keySet()) {
			double count=counts.get(position);
			double norm=count/avg;
			writer.write("chrX\t"+position+"\t"+(position+1)+"\t"+norm+"\n");
		}
		
		writer.close();
	}
	
	/*public static void main(String[] args) throws IOException {
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalsFromFile(args[0]);
		String save=args[1];
		new BedToBedgraph(regions, save);
	}*/
	
}

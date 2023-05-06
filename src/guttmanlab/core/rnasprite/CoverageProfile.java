package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class CoverageProfile {

	public CoverageProfile(BarcodingDataStreaming data, String save, int resolution) throws IOException {
		Map<SingleInterval, Double> counts=new TreeMap<SingleInterval, Double>();
		
		while(data.hasNext()) {
			Cluster c=data.next();
			Cluster binned=c.bin(resolution);
			for(SingleInterval region: binned.getAllDNAIntervals()) {
				double count=0;
				if(counts.containsKey(region)) {count=counts.get(region);}
				count++;
				counts.put(region, count);
			}
		}
		
		BEDFileIO.writeBEDGraph(counts, save);
	}
	
	public static void main(String[] args) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		int resolution=Integer.parseInt(args[2]);
		new CoverageProfile(data, save, resolution);
	}
	
}

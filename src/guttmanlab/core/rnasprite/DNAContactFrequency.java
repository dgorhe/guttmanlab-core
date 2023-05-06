package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class DNAContactFrequency {

	int resolution=1000;
	
	public DNAContactFrequency(BarcodingDataStreaming data, SingleInterval region, String save) throws IOException {
		Map<SingleInterval, Integer> map=new TreeMap<SingleInterval, Integer>();
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			if(c.containsOverlappingDNA(region)) {
				Cluster binned=c.bin(resolution);
				for(SingleInterval bin: binned.getAllDNAIntervals()) {
					int count=0;
					if(map.containsKey(bin)) {count=map.get(bin);}
					count++;
					map.put(bin, count);
				}
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		data.close();
		
		Collection<SingleInterval> local=region.allBins(resolution);
		for(SingleInterval l: local) {map.remove(l);}
		
		BEDFileIO.writeBEDGraphInteger(map, save);
	}

	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			SingleInterval region=new SingleInterval(args[1]);
			String save=args[2];
			new DNAContactFrequency(data, region, save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=barcodes \n args[1]=region \n args[2]=save";
}

package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class ComputeFrequencyOfBins {
	
	int maxClusterSize=20;
	
	public ComputeFrequencyOfBins(BarcodingDataStreaming data, String save) throws IOException{
		Map<SingleInterval, Collection<String>> map=getRegionMap(data);
		write(save, map);
	}

	private void write(String save, Map<SingleInterval, Collection<String>> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval interval:map.keySet()){
			Collection<String> names=map.get(interval);
			writer.write(interval.toUCSC()+"\t"+names.size()+"\n");
		}
		
		writer.close();
	}

	private Map<SingleInterval, Collection<String>> getRegionMap(BarcodingDataStreaming data) {
		Map<SingleInterval, Collection<String>> map=new TreeMap<SingleInterval, Collection<String>>();
		//Get all clusters overlapping a region
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getSize()>1 && c.getSize()<maxClusterSize){
			Collection<SingleInterval> intervals=c.getAllIntervals();
			for(SingleInterval interval: intervals){
				Collection<String> clusters=new TreeSet<String>();
				if(map.containsKey(interval)){clusters=map.get(interval);}
					clusters.add(c.getBarcode());
					map.put(interval, clusters);
				}	
			}
			counter++;
			if(counter %1000000 ==0){System.err.println(counter);}
			}
		return map;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File file=new File(args[0]);
			String save=args[1];
			int resolution=new Integer(args[2]);
			BarcodingDataStreaming data=new BarcodingDataStreaming(file, resolution);
			new ComputeFrequencyOfBins(data, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=barcoding data \n args[1]=save \n args[2]=resolution";
	
}

package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class CoveragePlot {
	int resolution=100000;

	public CoveragePlot(BarcodingDataFileTree data, Collection<SingleInterval> regions, String save) throws IOException{
		Collection<Cluster> clusters=data.getClustersOverlappingMultipleRegions(regions);
		write(save, clusters);
	}

	private void write(String save, Collection<Cluster> clusters) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		
		for(Cluster c: clusters){
			Cluster binned=c.bin(resolution);
			for(SingleInterval region: binned.getAllIntervals()){
				int count=0;
				if(counts.containsKey(region)){count=counts.get(region);}
				count++;
				counts.put(region, count);
			}
		}
		
		for(SingleInterval region: counts.keySet()){
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+counts.get(region)+"\n");			
		}
		
		writer.close();
	}
	
	private static SingleInterval parse(String string) {
		return new SingleInterval(string);
	}
	
	public static void main(String[] args)throws IOException{
		Collection<Cluster> triples=new BarcodingData(new File(args[1])).getClusters();	
		String saveDir=args[2];
		
		for(Cluster cluster: triples){
			System.err.println(cluster);
			BarcodingDataFileTree data=new BarcodingDataFileTree(new File(args[0]).listFiles());
			String save=saveDir+"/"+cluster.toFileName()+".bedgraph";
			new CoveragePlot(data, cluster.getAllIntervals(), save);
		}
	}

	
	
}

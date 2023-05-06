package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class BEDGraphOfInteractions {
	
	int resolution=1000000;
	int maxSize=1000;

	public BEDGraphOfInteractions(BarcodingDataStreaming data, SingleInterval hub,  String save) throws IOException{
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		int countOverlap=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()<maxSize){
				Cluster cluster=c.bin(resolution);
				if(cluster.overlapsInterval(hub)){
					countOverlap++;
					System.err.println(cluster.getBarcode()+" "+countOverlap+" "+counter);
					for(SingleInterval region: cluster.getAllIntervals()){
						int count=0;
						if(counts.containsKey(region)){count=counts.get(region);}
						count++;
						counts.put(region, count);
					}
				}
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		
		write(save, counts);
		
		/*for(SingleInterval region: hub.getAllIntervals()){
			System.err.println(region.toUCSC());
			Collection<Cluster> clusters=data.getClustersOverlappingRegion(region);
			write(save+"."+region.getReferenceName()+"_"+region.getReferenceStartPosition()+".bedgraph", clusters, region);
		}*/
		
	
		
		
	}
	
	
	private void write(String save, Map<SingleInterval, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: counts.keySet()){
			String chr=r.getReferenceName();
			if(chr.contains("_mouse")){
				String newChr=chr.replaceAll("_mouse", "");
				writer.write(newChr+"\t"+r.getReferenceStartPosition()+"\t"+r.getReferenceEndPosition()+"\t"+counts.get(r)+"\n");
			}
		}
		
		writer.close();
		
	}


	private void write(String save, Collection<Cluster> clusters) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		
		for(Cluster c: clusters){
			for(SingleInterval region: c.getAllIntervals()){
				writer.write(region.toBED()+"\n");
			}
		}
		
		
		writer.close();
	}

	
	
	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		
		SingleInterval region=new SingleInterval(args[1]);
		
		
		
		String save=args[2];
		new BEDGraphOfInteractions(data, region, save);
	}
	
}

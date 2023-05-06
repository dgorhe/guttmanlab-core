package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;

public class EnhancersTranscription {
	
	int resolution=100;

	public EnhancersTranscription(BarcodingDataStreaming data, SingleInterval promoterRegion, String gene, String save) throws IOException{
		//get all clusters overlapping promoter region
		Collection<Cluster> clusters=data.getClustersOverlappingRegion(promoterRegion);
		System.err.println("Retained "+clusters.size());
		
		//split by presence and absence of RNA
		Collection<Cluster> clustersWithRNA=new ArrayList<Cluster>();
		Collection<Cluster> clustersWithoutRNA=new ArrayList<Cluster>();
		
		for(Cluster c: clusters){
			if(c.containsRNA(gene)){clustersWithRNA.add(c);}
			else{clustersWithoutRNA.add(c);}
		}
		
		
		//make bedgraph
		writeBedgraph(save+".minusRNA.bedgraph", clustersWithoutRNA, resolution);
		writeBedgraph(save+".plusRNA.bedgraph", clustersWithRNA, resolution);
		
	}

	private void writeBedgraph(String save, Collection<Cluster> clustersWithoutRNA, int resolution) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<SingleInterval, Integer> map=new TreeMap<SingleInterval, Integer>();
		
		for(Cluster c: clustersWithoutRNA){
			Cluster binned=c.bin(resolution);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				int count=0;
				if(map.containsKey(region)){count=map.get(region);}
				count++;
				map.put(region, count);
			}
		}
		
		for(SingleInterval region: map.keySet()){
			int score=map.get(region);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+score+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		SingleInterval region=new SingleInterval(args[1]);
		String gene=args[2];
		String save=args[3];
		new EnhancersTranscription(data, region, gene, save);
	}
	
}

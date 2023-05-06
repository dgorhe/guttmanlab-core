package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.Cluster;

public class ClustersToCoverage {

	public ClustersToCoverage(BarcodingDataStreaming data, String save, int binSize) throws IOException{
		//Map<String, Map<SingleInterval, Double>> allMaps=new TreeMap<String, Map<SingleInterval, Double>>();
		
		Map<SingleInterval, Double> map=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			//Map<String, Double> proteins=c.getProteinWeights();
			c=c.bin(binSize);
			
			for(SingleInterval region: c.getAllRNARegions()) {
				//System.err.println(region.toUCSC());
				double count=0;
				if(map.containsKey(region)){count=map.get(region);}
				count++;
				map.put(region, count);
			}
			
			
			
			/*for(String protein: proteins.keySet()) {
				if(!allMaps.containsKey(protein)) {
					allMaps.put(protein, new TreeMap<SingleInterval, Double>());
				}
				double weight=proteins.get(protein);
				Map<SingleInterval, Double> map=allMaps.get(protein);
				for(SingleInterval region: c.getAllRNARegions()){
					double count=0;
					if(map.containsKey(region)){count=map.get(region);}
					count+=weight;
					map.put(region, count);
				}
			}*/
			counter++;
			//if(counter%1 ==0) {System.err.println(counter);}
		}
		
		data.close();
		BEDFileIO.writeBEDGraph(map, save);
		//write(save, map);
	}

	private void write(String save, Map<String, Map<SingleInterval, Double>> allMaps) throws IOException {
		for(String protein: allMaps.keySet()) {
			FileWriter writer=new FileWriter(save+"."+protein+".bedgraph");
			Map<SingleInterval, Double> map=allMaps.get(protein); 
			for(SingleInterval region: map.keySet()){
				writer.write(region.toBedgraph(map.get(region))+"\n");
			}
			
			writer.close();
		}
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		int binSize=Integer.parseInt(args[2]);
		new ClustersToCoverage(data, save, binSize);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=cluster file \n args[1]=save \n args[2]=bin size";
}

package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.rnasprite.Cluster;

public class ClustersToUniqueness {

	public ClustersToUniqueness(BarcodingDataStreaming data, String save) throws IOException{
		//Map<String, Map<SingleInterval, Double>> allMaps=new TreeMap<String, Map<SingleInterval, Double>>();
		
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			Map<String, Double> proteins=c.getProteinWeights();
			double maxPercent=getMaxPercent(proteins);
			writer.write(c.getBarcode()+"\t"+maxPercent+"\n");
			
			/*
			c=c.bin(binSize);
			for(String protein: proteins.keySet()) {
				if(!allMaps.containsKey(protein)) {
					allMaps.put(protein, new TreeMap<SingleInterval, Double>());
				}
				double weight=proteins.get(protein);
				Map<SingleInterval, Double> map=allMaps.get(protein);
				for(SingleInterval region: c.getAllDNAIntervals()){
					double count=0;
					if(map.containsKey(region)){count=map.get(region);}
					count+=weight;
					map.put(region, count);
				}
			}*/
			counter++;
			if(counter%100000 ==0) {System.err.println(counter);}
		}
		writer.close();
		data.close();
		//write(save, allMaps);
	}

	private double getMaxPercent(Map<String, Double> proteins) {
		double max=0;
		for(String p: proteins.keySet()) {
			max=Math.max(max, proteins.get(p));
		}
		return max;
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
		if(args.length>1){
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		//int binSize=Integer.parseInt(args[2]);
		new ClustersToUniqueness(data, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=cluster file \n args[1]=save";
}

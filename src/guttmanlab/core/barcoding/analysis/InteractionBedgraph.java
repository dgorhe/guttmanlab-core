package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;

public class InteractionBedgraph {

	public InteractionBedgraph(BarcodingDataStreaming data, String save, String gene) throws IOException {
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		while(data.hasNext()){
			Cluster c=data.next();
			int size=c.getClusterSize();
			if(size>1){
				if(contains(c, gene)){
					for(SingleInterval region: c.getAllIntervals()){
						int count=0;
						if(counts.containsKey(region)){count=counts.get(region);}
						count++;
						counts.put(region, count);
					}
				}
			}
		}
		write(save, counts);
	}

	private boolean contains(Cluster c, String gene) {
		for(SingleInterval region: c.getAllIntervals()){
			if(region.getReferenceName().contains(gene)){return true;}
		}
		return false;
	}

	private void write(String save, Map<SingleInterval, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval interval: counts.keySet()){
			int count=counts.get(interval);
			writer.write(interval.getReferenceName()+"\t"+interval.getReferenceStartPosition()+"\t"+interval.getReferenceEndPosition()+"\t"+count+"\n");
		}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		if(args.length>2){
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		String gene=args[2];
		new InteractionBedgraph(data, save, gene);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=cluster file \n args[1]=save name \n args[2]=gene name";
}

package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;

public class RNARNABEDGraph {

	public RNARNABEDGraph(BarcodingDataStreaming data, String rna, int resolution, String save) throws IOException{
		
		Map<SingleInterval, Double> scores=new TreeMap<SingleInterval, Double>();
		Map<SingleInterval, Double> inputScores=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getRNANames().contains(rna)){
				for(SingleInterval region:c.getAllRNARegions()){
					SingleInterval binned=region.bin(resolution);
					double score=0;
					if(scores.containsKey(binned)){score=scores.get(binned);}
					score++;
					scores.put(binned, score);
				}
			}
			for(SingleInterval region: c.getAllRNARegions()){
				SingleInterval binned=region.bin(resolution);
				double score=0;
				if(inputScores.containsKey(binned)){score=inputScores.get(binned);}
				score++;
				inputScores.put(binned, score);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		data.close();
		writeBedgraph(save, scores);
		writeBedgraph(save+".input.bedgraph", inputScores);
	}
	
	
	private void writeBedgraph(String save, Map<SingleInterval, Double> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: scores.keySet()){
			writer.write(region.toBedgraph(scores.get(region))+"\n");
		}
		
		writer.close();
	}


	public static void main(String[] args) throws IOException{
		if(args.length>3){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String rna=args[1];
			int resolution=new Integer(args[2]);
			String save=args[3];
			new RNARNABEDGraph(data, rna, resolution, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=barcoding data \n args[1]=rna name \n args[2]=resolution \n args[3]=save";
	
}

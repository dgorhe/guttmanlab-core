package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

public class GetClusterSizes {

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		double totalReadNumber=0;
		
		Map<Integer, Double> map=new TreeMap<Integer, Double>();
		
		int counter=0;
		while(data.hasNext() && counter<new Integer(args[2])){
			Cluster c=data.next();
			int clusterSize=c.getClusterSize();
			totalReadNumber+=clusterSize;
			double count=0;
			if(map.containsKey(clusterSize)){
				count=map.get(clusterSize);
			}
			count=count+clusterSize;
			map.put(clusterSize, count);
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		data.close();
		
		write(map, totalReadNumber, args[1]);
	}

	private static void write(Map<Integer, Double> map, double totalReadNumber, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Integer clusterSize: map.keySet()){
			double fraction=map.get(clusterSize)/totalReadNumber;
			writer.write(clusterSize+"\t"+fraction+"\t"+map.get(clusterSize)+"\t"+totalReadNumber+"\n");
		}
		
		writer.close();
	}
	
}

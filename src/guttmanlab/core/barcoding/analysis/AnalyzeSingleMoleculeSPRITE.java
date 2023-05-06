package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public class AnalyzeSingleMoleculeSPRITE {

	public AnalyzeSingleMoleculeSPRITE(File file, String save) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(file);
		
		Collection<Cluster> clusters=getClustersGreaterThan1(data);
		
		//enumerate all combinations
		Map<Cluster, Integer> counts=getCounts(clusters);
		
		write(save, counts);
		
		//count
		
	}

	private void write(String save, Map<Cluster, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Cluster c: counts.keySet()){
			writer.write(c.toStringNoName()+"\t"+counts.get(c)+"\n");
		}
		
		writer.close();
	}

	private Map<Cluster, Integer> getCounts(Collection<Cluster> clusters) {
		Map<Cluster, Integer> rtrn=new TreeMap<Cluster, Integer>();
		
		for(Cluster c: clusters){
			Cluster c2=new Cluster("c", c.getAllIntervals());
			int count=0;
			if(rtrn.containsKey(c2)){
				count=rtrn.get(c2);
			}
			count++;
			rtrn.put(c2, count);
		}
		
		return rtrn;
	}

	private Collection<Cluster> getClustersGreaterThan1(BarcodingDataStreaming data) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			Cluster binned=c.bin(100);
			if(binned.getClusterSize()>1){rtrn.add(binned);}
		}
		
		data.close();
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new AnalyzeSingleMoleculeSPRITE(file, save);
	}
	
}

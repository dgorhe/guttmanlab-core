package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

public class GetGeneCounts {
	
	private static void write(Map<String, Integer> counts) {
		for(String gene: counts.keySet()) {
			System.out.println(gene+"\t"+counts.get(gene));
		}
		
	}
	
	private static Map<String, Integer> getCounts(BarcodingDataStreaming data) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Collection<String> names= c.getRNANames();
			for(String name: names) {
				int count=0;
				if(rtrn.containsKey(name)) {count=rtrn.get(name);}
				count++;
				rtrn.put(name, count);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		data.close();
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		Map<String, Integer> counts=getCounts(data);
		write(counts);
	}

	

	
	
}

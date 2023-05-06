package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;

public class PullUnannotated {

	public PullUnannotated(BarcodingDataStreaming data, String save) throws IOException{
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		while(data.hasNext()){
			Cluster c=data.next();
			//if(c.containsRNA("Unassigned_NoFeatures") || c.containsRNA("Unassigned_Ambiguity")){
				//System.err.println(c.getBarcode());
				for(SingleInterval region:c.getAllRNARegions()){
					if(region.getName().contains("Unassigned")){
						//System.err.println(region.toUCSC()+" "+region.getName());
						int count=0;
						if(rtrn.containsKey(region)){count=rtrn.get(region);}
						count++;
						rtrn.put(region, count);
					}
				}
			//}
		}
		data.close();
		write(save, rtrn);
	}

	private void write(String save, Map<SingleInterval, Integer> rtrn) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: rtrn.keySet()){
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+rtrn.get(region)+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		new PullUnannotated(data, save);
	}
	
}

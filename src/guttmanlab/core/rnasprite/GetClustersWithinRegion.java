package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import guttmanlab.core.annotation.SingleInterval;

public class GetClustersWithinRegion {

	public GetClustersWithinRegion(BarcodingDataStreaming data, SingleInterval region, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			if(c.containsOverlappingDNA(region)) {
				writer.write(c.line+"\n");
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		data.close();
			
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			SingleInterval region=new SingleInterval(args[1]);
			String save=args[2];
			new GetClustersWithinRegion(data, region, save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=barcodes \n args[1]=region \n args[2]=save";
	
	
}

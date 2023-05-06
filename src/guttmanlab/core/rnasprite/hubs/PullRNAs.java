package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;

public class PullRNAs {

	
	public PullRNAs(BarcodingDataStreaming data, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			if(c.getRNAClusterSize()>=2) {
				writer.write(c.toRNAString()+"\n");
			}
			counter++;
			if(counter%100000 ==0) {System.err.println(counter);}
		}
		
		data.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		new PullRNAs(data, save);
	}
	
}

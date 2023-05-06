package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class RNACoveragePlot {

	
	public RNACoveragePlot(BarcodingDataStreaming data, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		while(data.hasNext()) {
			Cluster c=data.next();
			for(RNAInterval r: c.getAllRNARegions()) {
				writer.write(r.toBED()+"\n");
			}
		}
		
		data.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		new RNACoveragePlot(data, save);
	}
	
}

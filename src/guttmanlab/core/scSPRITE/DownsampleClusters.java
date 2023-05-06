package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.rnasprite.Cluster;

public class DownsampleClusters {

	public DownsampleClusters(BarcodingDataStreaming data, double fraction, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			double random=Math.random();
			if(random<fraction){
				write(writer,c);
			}
			counter++;
			if(counter%100000000 ==0){System.err.println(counter);}
		}
		
		data.close();
		writer.close();
	}

	private void write(FileWriter writer, Cluster c) throws IOException {
		writer.write(c.getBarcode());
		
		for(SingleInterval region: c.getAllDNAIntervals()){
			writer.write("\t"+region.getReferenceName()+":"+region.getReferenceStartPosition());
		}
		writer.write("\n");
	}
	
	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		double fraction=new Double(args[1]);
		String save=args[2];
		new DownsampleClusters(data, fraction, save);
	}
	
}

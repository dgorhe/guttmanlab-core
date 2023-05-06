package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;

public class MultiwayVisualization {

	public MultiwayVisualization(BarcodingDataStreaming data, String save, int binResolution) throws IOException{
		
		FileWriter writer=new FileWriter(save);
		
		while(data.hasNext()){
			Cluster c= data.next();
			Annotation a=makeAnnotation(c.getAllDNAIntervals(), c.getBarcode(), binResolution);
			writer.write(a.toBED()+"\n");
		}
		
		data.close();
		writer.close();
	}

	private BlockedAnnotation makeAnnotation(Collection<SingleInterval> allDNAIntervals, String name, int binResolution) {
		Collection<SingleInterval> binned=new TreeSet<SingleInterval>();
		for(SingleInterval interval: allDNAIntervals){binned.add(interval.bin(binResolution));}
		
		BlockedAnnotation a=new BlockedAnnotation(binned, name);
		return a;
	}
	
	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		int binResolution=new Integer(args[2]);
		new MultiwayVisualization(data, save, binResolution);
	}
	
}

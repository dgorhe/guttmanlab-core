package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.SingleInterval;

public class BarcodeToCoverage {

	public BarcodeToCoverage(BarcodingDataStreaming data, SingleInterval region){
		
		Collection<Cluster> clusters=data.getClustersOverlappingRegion(region, 1, 1000);
		
		System.err.println(region.toUCSC()+" "+clusters.size());
	}
	
	public static void main(String args[]) throws IOException{
		File file=new File(args[0]);
		
		new BarcodeToCoverage(new BarcodingDataStreaming(file), new SingleInterval(args[1]));
	}
	
}

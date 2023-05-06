package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import guttmanlab.core.annotation.SingleInterval;

public class EnumerateTriples {

	
	public static void main(String[] args) throws IOException{
		
		if(args.length>2){
			File file=new File(args[0]);
			int resolution=new Integer(args[1]);
			int maxClusterSize=new Integer(args[2]);
			boolean debugMode=new Boolean(args[3]);
			
			BarcodingDataStreaming data=new BarcodingDataStreaming(file, resolution);
			Map<Cluster, Set<String>> clusters=data.enumerateAllKmers(3, maxClusterSize);
			write(clusters, debugMode);
		}
		else{System.err.println(usage);}
	}
	
	private static void write( Map<Cluster, Set<String>> clusters, boolean debug) throws IOException {
		if(debug){
			for(Cluster cluster: clusters.keySet()){
				Set<String> barcodes=clusters.get(cluster);
				for(SingleInterval interval: cluster.getAllIntervals()){
					System.out.print("\t"+interval.toUCSC());
				}
				System.out.print("\t"+barcodes.toString()+"\n");
				
			}
			System.out.flush();
		}
		
		else{
			for(Cluster cluster: clusters.keySet()){
				int count=clusters.get(cluster).size();
				System.out.print(cluster.getBarcode());
				for(SingleInterval interval: cluster.getAllIntervals()){
					System.out.print("\t"+interval.toUCSC());
				}
				System.out.print("\t"+count+"\n");
				
			}
			System.out.flush();
		}
		
		
	}

	static String usage=" args[0]=barcode file \n args[1]=resolution \n args[2]=max cluster size to consider \n args[3]=debug mode (true or false)";
}

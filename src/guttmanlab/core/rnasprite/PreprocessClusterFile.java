package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;

public class PreprocessClusterFile {

	public PreprocessClusterFile(BarcodingDataStreaming data, String savedir) throws IOException{
		//Chunk into RNA and DNA - DNA by 1Mb regions
		
		Collection<String> allRNAs=new ArrayList<String>();
		Map<SingleInterval, FileWriter> dnaWriter=new TreeMap<SingleInterval, FileWriter>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()<1000){
				Cluster binned=c.bin(1000000);
				for(SingleInterval region: binned.getAllDNAIntervals()){
					FileWriter writer;
					if(dnaWriter.containsKey(region)){
						writer=dnaWriter.get(region);
					}
					else{
						writer=new FileWriter(savedir+"/"+region.getFileName(), true);
						dnaWriter.put(region, writer);
					}
					writer.write(c.toString()+"\n");
					
				}
				counter++;
				if(counter%10000 ==0){
					System.err.println(counter);
					close(dnaWriter);
					dnaWriter=new TreeMap<SingleInterval, FileWriter>();
				}
			}
		}
		
		data.close();
		close(dnaWriter);
		
	}
	
	private void close(Map<SingleInterval, FileWriter> dnaWriter) throws IOException {
		for(SingleInterval region: dnaWriter.keySet()){
			dnaWriter.get(region).close();
		}
	}

	public static void main(String[] args)throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String saveDir=args[1];
		new PreprocessClusterFile(data, saveDir);
	}
	
}

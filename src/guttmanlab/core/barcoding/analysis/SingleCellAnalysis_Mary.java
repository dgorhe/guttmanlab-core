package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class SingleCellAnalysis_Mary {

	static int resolution=1000000;
	
	public static void main (String[] args) throws IOException{
		
		if(args.length>2){
			File[] files=new File(args[0]).listFiles();
			String saveDir=args[1];
			int minClusterSize=new Integer(args[2]);
			new File(saveDir).mkdir();
			
			for(int i=0; i< files.length; i++){
				//System.err.println(files[i].getName());
				BarcodingDataStreaming data=new BarcodingDataStreaming(files[i]);
				Collection<SingleInterval> allIntervals=new TreeSet<SingleInterval>();
				
				while(data.hasNext()){
					Cluster c=data.next();
					c.bin(resolution);
					allIntervals.addAll(c.getAllIntervals());
					
				}
				data.close();
				
				if(allIntervals.size()>minClusterSize){
					System.err.println(files[i].getName()+" "+allIntervals.size());
					write(saveDir+"/"+files[i].getName()+".bed", allIntervals);
				}
				
			}
			
			
		}
		else{
			System.err.println(line);
		}
	}
	
	static String line=" args[0]=files (clusters) \n args[1]=save directory \n args[2]=min cluster size to consider";

	private static void write(String save, Collection<SingleInterval> allIntervals) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval interval: allIntervals){
			writer.write(interval.getReferenceName()+"\t"+interval.getReferenceStartPosition()+"\t"+interval.getReferenceEndPosition()+"\n");	
		}
		writer.close();
	}	
}
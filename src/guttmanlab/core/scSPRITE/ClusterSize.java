package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math3.util.CombinatoricsUtils;

import guttmanlab.core.rnasprite.Cluster;

public class ClusterSize {
	
	public ClusterSize(File[] files) throws IOException{
		double count=0;
		double total=0;
		
		//FileWriter writer=new FileWriter(save);
		for(int i=0; i<files.length; i++){
			BarcodingDataStreaming data=new BarcodingDataStreaming(files[i]);
			
			//long numContacts=data.getNumberOfContacts(maxClusterSize, binSize);
			
			long numContacts=0;
			
			while(data.hasNext()){
				Cluster c=data.next();
				if(c.getClusterSize()>10000){
					count++;
				}
				else if(c.getClusterSize()>1){numContacts+=CombinatoricsUtils.binomialCoefficient(c.getClusterSize(), 2);}
				total++;
				//writer.write(c.getClusterSize()+"\n");
			}
			System.err.println(i+"\t"+numContacts);
		}
		
		
		System.err.println(count+" "+total+" "+(count/total));
		
		//writer.close();
	}

	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		//String save=args[1];
		new ClusterSize(files);
	}
	
}

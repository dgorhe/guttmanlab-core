package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

public class FilterBarcodeFile {

	public static void main(String[] args) throws IOException{
		if(args.length>4){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			int minClusterSize=new Integer(args[1]);
			int maxClusterSize=new Integer(args[2]);
			String save=args[3];
			
			FileWriter writer=new FileWriter(save);
			while(data.hasNext()){
				Cluster c=data.next();
				int size=c.getSize();
				if(size>minClusterSize && size<maxClusterSize){
					writer.write(c.toString()+"\n");
				}
				
			}
			
			data.close();
			writer.close();
			
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=clusters \n args[1]=min cluster \n args[2]=max cluster size \n args[3]=bin resolution \n args[4]=save";
}

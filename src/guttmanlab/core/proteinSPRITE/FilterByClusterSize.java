package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import guttmanlab.core.rnasprite.Cluster;

public class FilterByClusterSize {

	
	private static void filter(BarcodingDataStreaming data, int minSize, int maxSize, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		int counter=0; 
		int filtered=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			if(c.getProteins().size()>=minSize && c.getProteins().size()<=maxSize) {
				writer.write(c.toString()+"\n");
				filtered++;
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+filtered);}
		}
		data.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>3) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			
			int minSize=Integer.parseInt(args[1]);
			int maxSize=Integer.parseInt(args[2]);
			
			String save=args[3];
			
			filter(data, minSize, maxSize, save);
		}
		else {System.err.println(usage);}
		
	}
	
	

	static String usage=" args[0]=clusters \n args[1]=min \n args[2]=max \n args[3]=save";
}

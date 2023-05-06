package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class FilterByClusterSize {
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>3) {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		FileWriter writer=new FileWriter(args[1]);
		int minSize=Integer.parseInt(args[2]);
		int maxSize=Integer.parseInt(args[3]);
		
		int counter=0;
		int retained=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()>=minSize && c.getClusterSize()<=maxSize){
				writer.write(c.getLine()+"\n");
				//writer.write(c.toString()+"\n");
				retained++;
			}
			counter++;
			if(counter%100000==0){System.err.println(counter+" "+retained);}
		}
		writer.close();
		data.close();
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=save \n args[2]=min size \n args[3]=max size";
	
}

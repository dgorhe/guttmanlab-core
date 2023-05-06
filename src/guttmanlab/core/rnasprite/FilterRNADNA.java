package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class FilterRNADNA {

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		FileWriter writer=new FileWriter(args[1]);
		int counter=0;
		int total=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.hasDNA() && c.hasRNA()){
				writer.write(c.toString()+"\n");
				counter++;
			}
			total++;
			if(total% 1000000 ==0){System.err.println(total+" "+counter);}
		}
		System.err.println("Total "+total+" Filtered to "+counter);
		data.close();
		writer.close();
	}
	
}

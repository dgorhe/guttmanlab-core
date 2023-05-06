package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

public class GetAllRNANames {

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		Collection<String> names=new TreeSet<String>();
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			names.addAll(c.getRNANames());
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		data.close();
		FileWriter writer=new FileWriter(args[1]);
		
		for(String name: names){
			writer.write(name+"\n");
		}
		
		writer.close();
	}
	
}

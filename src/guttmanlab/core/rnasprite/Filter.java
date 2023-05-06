package guttmanlab.core.rnasprite;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import guttmanlab.core.annotation.io.BEDFileIO;

public class Filter {

	public static void main(String[] args) throws IOException{
		FileWriter writer=new FileWriter(args[1]);
		List<String> lines=BEDFileIO.loadLines(args[0]);
		
		for(String line: lines){
			if(line.split("\t").length>2){
				writer.write(line+"\n");
			}
		}
		writer.close();
	}
	
}

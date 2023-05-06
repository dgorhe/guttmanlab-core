package guttmanlab.core.clap.old;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class ConvertToBEDGraph {

	public ConvertToBEDGraph(File file, String save) throws IOException{
		//chr:start-end+ score
		//chr start end score -- split by strand
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
			String[] tokens=nextLine.split("\t");
			String position=tokens[0].substring(0, tokens[0].length()-1);
			String chr=position.split(":")[0];
			String start=position.split(":")[1].split("-")[0];
			String end=position.split(":")[1].split("-")[1];
			
			writer.write(chr+"\t"+start+"\t"+end+"\t"+tokens[1]+"\n");
		}
		
	}
	
	public static void main (String[] args) throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new ConvertToBEDGraph(file, save);
	}
	
}

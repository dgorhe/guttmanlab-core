package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;

public class RenameClusterFile {

	public RenameClusterFile(BarcodingDataStreaming data, Map<String, String> names, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Cluster renamed=c.renameRNA(names);
			writer.write(renamed.toUniqueRNAString()+"\n");
			counter++;
			if(counter%1000000 ==0) {System.err.println(counter);}
		}
		
		data.close();
		writer.close();
	}
	
	private static Map<String, String> parse(String string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		List<String> lines= BEDFileIO.loadLines(string);
		
		for(String line: lines) {
			String unqiueName=line.split("\t")[1];
			String className=line.split("\t")[2];
			rtrn.put(unqiueName, className);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		Map<String, String> names=parse(args[1]);
		String save=args[2];
		new RenameClusterFile(data, names, save);
	}

	
	
}

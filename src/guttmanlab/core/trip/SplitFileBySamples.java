package guttmanlab.core.trip;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;

public class SplitFileBySamples {

	public static void main(String[] args) throws IOException {
		String input=args[0];
		String output=args[1];
		
		List<String> lines=BEDFileIO.loadLines(input, 1);
		
		Map<String, FileWriter> writers=new TreeMap<String, FileWriter>();
		
		for(String line: lines) {
			String[] tokens=line.split(",");
			String name=tokens[3];
			String x=tokens[2];
			String y=tokens[1];
			
			if(!writers.containsKey(name)) {
				writers.put(name, new FileWriter(output+"/"+name+".txt"));
			}
			
			FileWriter writer=writers.get(name);
			writer.write(x+"\t"+y+"\t"+name+"\n");
			
		}
		
		close(writers);
		
	}

	private static void close(Map<String, FileWriter> writers) throws IOException {
		for(String name: writers.keySet()) {
			writers.get(name).close();
		}
		
	}
	
}

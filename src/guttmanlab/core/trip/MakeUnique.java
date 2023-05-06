package guttmanlab.core.trip;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;

public class MakeUnique {

	public static void main(String[] args) throws IOException {
		List<String> lines=BEDFileIO.loadLines(args[0]);
		Map<String, String> set=new TreeMap<String, String>();
		for(String line: lines) {
			String barcode=line.split("\t")[0];
			set.put(barcode, line);
		}
		
		FileWriter writer=new FileWriter(args[1]);
		
		for(String barcode: set.keySet()) {
			writer.write(set.get(barcode)+"\n");
		}
		
		writer.close();
	}
	
}

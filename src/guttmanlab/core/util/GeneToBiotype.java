package guttmanlab.core.util;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

public class GeneToBiotype {

	public static void main(String[] args) throws IOException {
		Map<String, Collection<String>> map=convert(args[0]);
		write(args[1], map);
	}

	private static void write(String save, Map<String, Collection<String>> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(String name: map.keySet()) {
			writer.write(name);
			Collection<String> types=map.get(name);
			for(String type: types) {writer.write("\t"+type);}
			writer.write("\n");
		}
		writer.close();
	}

	private static Map<String, Collection<String>> convert(String string) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(string)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String info=tokens[8];
			String geneName=get(" gene_name", info);
			String biotype=get(" gene_biotype", info);
			if(!rtrn.containsKey(geneName)) {rtrn.put(geneName, new TreeSet<String>());}
			Collection<String> set=rtrn.get(geneName);
			set.add(biotype);
		}
		reader.close();
		return rtrn;
	}

	private static String get(String key, String line) {
		String[] tokens=line.split("\\;");
		for(int i=0; i<tokens.length; i++) {
			//System.err.println(tokens[i]);
			if(tokens[i].startsWith(key)) {
				//System.err.println(tokens[i].split(" ")[2]);
				String rtrn=tokens[i].split(" ")[2];
				rtrn=rtrn.replaceAll("\"", "");
				return rtrn;
			}
		}
		System.err.println(line);
		
		return null;
	}
	
}

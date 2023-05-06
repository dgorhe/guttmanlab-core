package guttmanlab.core.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;

public class MakeTable {
	
	private static void mergeAndWrite(Map<String, Collection<String>>[] maps, File[] files, String save) throws IOException {
		Collection<String> allRNAs=new TreeSet<String>();
		for(int i=0; i<maps.length; i++) {
			allRNAs.addAll(maps[i].keySet());
		}
		
		FileWriter writer=new FileWriter(save);
		
		writer.write("RNA");
		for(int i=0; i<files.length; i++) {
			writer.write("\t"+files[i].getName());
		}
		writer.write("\n");
		
		for(String rna: allRNAs) {
			writer.write(rna);
			for(int i=0; i<maps.length; i++) {
				Collection<String> vals=get(maps[i], rna);
				String val=toString(vals);
				writer.write("\t"+val);
			}
			writer.write("\n");
		}
		
		
		
		writer.close();
	}

	private static Collection<String> get(Map<String, Collection<String>> map, String rna) {
		Collection<String> rtrn=new TreeSet<String>();
		if(map.containsKey(rna)) {
			rtrn=map.get(rna);
		}
		return rtrn;
	}
	
	private static Map<String, Collection<String>>[] parse(File[] files) throws IOException {
		Map<String, Collection<String>>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			rtrn[i]=parse(files[i]);
		}
		
		return rtrn;
	}

	private static Map<String, Collection<String>> parse(File file) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens=nextLine.split("\t");
			String rna=tokens[0];
			Collection<String> list=new TreeSet<String>();
			for(int i=4; i<tokens.length; i++) {
				list.add(tokens[i]);
			}
			rtrn.put(rna, list);
			
		}
		reader.close();
		
		return rtrn;
	}

	private static String toString(Collection<String> vals) {
		String rtrn="";
		for(String val: vals) {
			rtrn+=val+";";
		}
		if(rtrn.endsWith(";")) {
			rtrn=rtrn.substring(0, rtrn.length()-1);
		}
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		/*File[] files=new File(args[0]).listFiles();
		String save=args[1];
		Map<String, Collection<String>>[] maps=parse(files);
		mergeAndWrite(maps, files, save);*/
		
		parse(args[0], args[1]);
		
	}

	private static void parse(String in, String out) throws IOException {
		FileWriter writer=new FileWriter(out);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(new File(in))));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens=nextLine.split("\t");
			String rna=tokens[0];
			String list=tokens[1];
			
			String[] genes=list.split(";");
			writer.write(rna);
			for(int i=0; i<genes.length; i++) {
				if(!genes[i].startsWith("chr")) {
					writer.write("\t"+genes[i]);
				}
			}
			writer.write("\n");
			
			//rtrn.put(rna, list);
			
		}
		reader.close();
		
		writer.close();
	}

	

	
	
}

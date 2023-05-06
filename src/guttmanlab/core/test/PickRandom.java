package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;

public class PickRandom {

	public static void main(String[] args) throws IOException{
		List<String> lines=BEDFileIO.loadLines(args[0]);
		FileWriter writer=new FileWriter(args[1]);
		
		Map<String, List<String>> linesByName=new TreeMap<String, List<String>>();
		
		for(String line: lines){
			String name=line.split("\t")[1];
			if(!name.equals("Parent_sequence")){
				List<String> list=new ArrayList<String>();
				if(linesByName.containsKey(name)){list=linesByName.get(name);}
				list.add(line);
				linesByName.put(name, list);
			}
		}
		
		for(String name: linesByName.keySet()){
			List<String> list=linesByName.get(name);
			System.err.println(name+" "+list.size());
			Collection<String> filtered=get(list);
			for(String l: filtered){
				writer.write(l+"\n");
			}
		}
		
		writer.close();
	}

	private static Collection<String> get(List<String> list) {
		Collection<String> rtrn=new ArrayList<String>();
		if(list.size()<2){return list;}
		rtrn.add(list.get(1));
		rtrn.add(list.get(list.size()-2));
		return rtrn;
	}
	
}

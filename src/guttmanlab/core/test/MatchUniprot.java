package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;

public class MatchUniprot {

	
	public static void main(String[] args) throws IOException{
		Map<String, String> uniProt=parse(args[0]);
		Map<String, String> chromatin=parse(args[1]);
		FileWriter writer=new FileWriter(args[2]);
		
		
		for(String key: chromatin.keySet()){
			String refSeq=uniProt.get(key);
			String line=chromatin.get(key);
			writer.write(line+"\t"+refSeq+"\n");
		}
		writer.close();
	}

	private static Map<String, String> parse(String string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		Collection<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			rtrn.put(line.split("\t")[0], line);
		}
		
		return rtrn;
	}
	
}

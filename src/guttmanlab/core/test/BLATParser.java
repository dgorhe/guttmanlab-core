package guttmanlab.core.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SingleInterval;

public class BLATParser {


	private static void parsePSL(File pslFile, int minMatch, String save) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(pslFile)));
		FileWriter writer=new FileWriter(save);
		Map<String, Collection<String>> map=new TreeMap<String, Collection<String>>();
		
		String nextLine;
		while((nextLine=reader.readLine())!=null){
			try{
				String[] tokens=nextLine.split("\t");
				int match=new Integer(tokens[0]);
				if(match>minMatch){
					String name=tokens[9];
					Collection<String> temp=new TreeSet<String>();
					if(map.containsKey(name)){temp=map.get(name);}
					temp.add(nextLine);
					map.put(name, temp);
				}
			}catch(NumberFormatException e){System.err.println("Excluded "+nextLine);}
		}
		
		//Get max
		for(String name: map.keySet()){
			Collection<String> lines=map.get(name);
			String maxLine="";
			int maxScore=0;
			for(String line: lines){
				String[] tokens=line.split("\t");
				int score=new Integer(tokens[0]);
				if(score>maxScore){maxScore=score; maxLine=line;}
			}
			writer.write(parseLine(maxLine)+"\n");
		}
		
		reader.close();
		writer.close();
	}

	
	private static String parseLine(String maxLine) {
		String[] tokens=maxLine.split("\t");
		String name=tokens[9];
		String match=tokens[0];
		String target=tokens[13];
		return name+"\t"+target+"\t"+match;
	}




	public static void main(String[] args)throws IOException{
		File psl=new File(args[0]);
		int minMatch=new Integer(args[1]);
		String save=args[2];
		parsePSL(psl, minMatch, save);
	}
}

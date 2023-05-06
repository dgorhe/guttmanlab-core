package guttmanlab.core.sars.leader;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.stat.inference.ChiSquareTest;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class VariantCount {
	
	Collection<String> proteins;
	
	public VariantCount (String file, String save) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file);
		
		Map<String, Double> freq=getFrequencies(lines);
		
		write(save, freq, lines.size()-1);
	}
	
	
	
	
	private void write(String save, Map<String, Double> freq, int total) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String pos: freq.keySet()) {
			double val=freq.get(pos)/(double)total;
			writer.write(pos+"\t"+val+"\n");
		}
		
		writer.close();
	}




	

	

	
	private Map<String, Double> getFrequencies(List<String> lines) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		
		int counter=0;
		for(String line: lines) {
			Collection<String> mutants=getMutants(line);
			for(String m1: mutants) {
				double count=0;
				if(rtrn.containsKey(m1)) {count=rtrn.get(m1);}
				count++;
				rtrn.put(m1, count);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		
		return rtrn;
	}
	

	


	private Collection<String> getMutants(String line) {
		String[] tokens=line.split("\t");
		Collection<String> mutants=parseMutants(tokens[5]);
		return mutants;
	}




	
	
	
	


	private Collection<String> parseMutants(String string) {
		Collection<String> rtrn=new TreeSet<String>();
		
		if(string.equals("()") || string.isEmpty()){return rtrn;}
		else if(!string.contains("(")) {return rtrn;}
		
		String trimmed=string.split("\\(")[1].split("\\)")[0];
			
		
		
		
		String[] tokens=trimmed.split(",");
		
		for(int i=0; i<tokens.length; i++) {
			
			if(tokens[i].split("_").length<2) {System.err.println(tokens[i]);} //TODO what's going on?
			else {
				String protein=tokens[i].split("_")[0];
				String pos=getNum(tokens[i].split("_")[1]);
				rtrn.add(protein+"_"+pos);
			}
		}
		
		return rtrn;
	}
	



	private String getNum(String string) {
		String rtrn=string;
		if(string.contains("del")) {
			rtrn=string.substring(1, string.length()-3);
		}
		else if(string.contains("stop")) {
			rtrn=string.substring(1, string.length()-4);
		}
		else {
			rtrn=string.substring(1, string.length()-1);
			//System.err.println(string+" "+rtrn);
		}
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException {
		new VariantCount(args[0], args[1]);
	}
	
}

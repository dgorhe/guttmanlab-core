package guttmanlab.core.sars.leader;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

public class MutationsPerProtein {
	
	
	public MutationsPerProtein (String file, String save) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file);
		Map<String, Integer> counts=computeByProtein(lines);
		write(save, counts);
	}
	

	private void write(String save, Map<String, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String protein: counts.keySet()) {
			writer.write(protein+"\t"+counts.get(protein)+"\n");
		}
		
		writer.close();
	}


	private Map<String, Integer> computeByProtein(List<String> lines) {
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			Collection<String> mutants=parseNames(tokens[5]);
			for(String mut: mutants) {
				int count=0;
				if(counts.containsKey(mut)) {count=counts.get(mut);}
				count++;
				counts.put(mut, count);
			}
			
		}
		
		
		return counts;
	}


	private Collection<String> parseNames(String string) {
		Collection<String> rtrn=new TreeSet<String>();
		
		if(string.equals("()") || string.isEmpty()){return rtrn;}
		else if(!string.contains("(")) {return rtrn;}
		
		String trimmed=string.split("\\(")[1].split("\\)")[0];
			
		
		
		
		String[] tokens=trimmed.split(",");
		
		for(int i=0; i<tokens.length; i++) {
			rtrn.add(tokens[i].split("_")[0]);
		}
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException {
		new MutationsPerProtein(args[0], args[1]);
	}
	
}

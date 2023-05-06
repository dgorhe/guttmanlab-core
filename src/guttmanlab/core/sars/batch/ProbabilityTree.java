package guttmanlab.core.sars.batch;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Kmer;

public class ProbabilityTree {

	
	static void writeProbabilityTree(BarcodingDataStreaming data, Map<String, String> allRNA, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		IntervalTree<String> probabilityTree=data.getProbabilityTree(allRNA);
		
		Iterator<Node<String>> iter=probabilityTree.iterator();
		
		while(iter.hasNext()) {
			Node<String> node=iter.next();
			writer.write(node.getValue()+"\t"+node.getStart()+"\t"+node.getEnd()+"\n");
		}
		writer.close();
	}
	
	
	public static IntervalTree<String> parseProbabilityTree(String input) throws IOException {
		IntervalTree<String> rtrn=new IntervalTree<String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String rna=tokens[0];
			int start=Integer.parseInt(tokens[1]);
			int end=Integer.parseInt(tokens[2]);
			rtrn.put(start, end, rna);
		}
		reader.close();
		return rtrn;
	}
	
	private static Map<String, String> parseRNAs(String string) throws IOException {
		List<String> lines= BEDFileIO.loadLines(string);
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String line: lines) {
			String rna=line.split("\t")[1].replaceAll("\"", "");
			String className=line.split("\t")[2];
			rtrn.put(rna, className);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Map<String, String> allRNAs=parseRNAs(args[1]);
			String save=args[2];
			
			writeProbabilityTree(data, allRNAs, save);
			
		}
		else {System.err.println(usage);}
	}
	

	private static String usage=" args[0]=clusters \n args[1]=all rnas \n args[2]=save";
	
}

package guttmanlab.core.splicing.speckle;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;


public class MergeSpeckleDistance4PL {

	private static Map<SingleInterval, String> parse4PL(String string) throws IOException {
		Map<SingleInterval, String> rtrn=new TreeMap<SingleInterval, String>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(string)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			SingleInterval gene=new SingleInterval(nextLine.split("\t")[1]);
			rtrn.put(gene, nextLine);
		}
		reader.close();
		
		return rtrn;
	}
	
	private static Map<String, String> parseSpeckle(String string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(string)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String gene=nextLine.split("\t")[0];
			gene=gene.replaceAll("\"","");
			rtrn.put(gene, nextLine);
		}
		reader.close();
		
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException {
		Map<String, String> specklesByName=parseSpeckle(args[0]);
		//Map<String, String> fitByName=parseSpeckle(args[1]);
		//Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[1]);
		
		FileWriter writer=new FileWriter(args[2]);
		
		List<String> lines=BEDFileIO.loadLines(args[1]);
		
		for(String line: lines) {
			String name=line.split("\t")[0];
			//String line1=fitByName.get(name);
			String line2=specklesByName.get(name);
			writer.write(line2+"\t"+line+"\n");
		}
		
		writer.close();
		
		/*FileWriter writer=new FileWriter(args[2]);
		for(SingleInterval gene: fitByName.keySet()) {
			if(overlaps(gene, regions)) {
				writer.write(fitByName.get(gene)+"\n");
			}
		}
		writer.close();
		*/
		
		
		
	}

	private static boolean overlaps(SingleInterval gene, Collection<SingleInterval> regions) {
		for(SingleInterval bin: regions) {
			if(gene.overlaps(bin)) {return true;}
		}
		return false;
	}

	
	
	
}

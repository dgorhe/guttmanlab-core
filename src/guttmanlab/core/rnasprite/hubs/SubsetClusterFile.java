package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;

public class SubsetClusterFile {

	public SubsetClusterFile(BarcodingDataStreaming data, Map<String, String> allRNA, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			if(containsAny(c, allRNA)) {
				write(writer, c);
			}
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		data.close();
		writer.close();
	}

	private boolean containsAny(Cluster c, Map<String, String> allRNA) {
		for(String rna: c.getRNANames()) {
			if(allRNA.containsKey(rna)) {return true;}
		}
		return false;
	}

	private void write(FileWriter writer, Cluster c) throws IOException {
		writer.write(c.toString()+"\n");
		
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
			Map<String, String> rnas=parseRNAs(args[1]);
			String save=args[2];
			new SubsetClusterFile(data, rnas, save);
		}
		
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=barcodes \n args[1]=all rnas \n args[2]=save";
}

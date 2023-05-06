package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class GetClusters {

	private static Map<String, Kmer> parseSets(String string) throws IOException {
		Map<String, Kmer> rtrn=new TreeMap<String, Kmer>();
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			String group=tokens[1];
			String name=tokens[0].replaceAll("\"", "");
			Kmer kmer=new Kmer();
			if(rtrn.containsKey(group)){
				kmer=rtrn.get(group);
			}
			kmer.addRegion(name);
			rtrn.put(group, kmer);
		}
		
		return rtrn;
	}
	
	private static Kmer parse(String string) {
		String[] tokens=string.split("_");
		Kmer rtrn=new Kmer();
		
		for(int i=0; i<tokens.length; i++){
			rtrn.addRegion(tokens[i]);
		}
		
		return rtrn;
	}
	
	private static void write(String string, Collection<Cluster> clusters) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Cluster c: clusters){
			writer.write(c.getLine()+"\n");
		}
		
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>4){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String save=args[1];
			Kmer kmer=parse(args[2]);
			Map<String, Kmer> collapseSet=parseSets(args[3]);
			int cutoff=new Integer(args[4]);
			
			Collection<Cluster> clusters=data.getClusters(kmer, collapseSet, cutoff);
			
			System.err.println("Retained clusters "+clusters.size());
			
			write(save+".clusters", clusters);
			
			data.close();
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=save \n args[2]=kmer (_) \n args[3]=collapse sets \n args[4]=number of RNAs";
	
}

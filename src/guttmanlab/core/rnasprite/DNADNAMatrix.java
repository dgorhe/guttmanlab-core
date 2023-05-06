package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class DNADNAMatrix {

	private static Kmer parse(String string) {
		String[] tokens=string.split("_");
		Kmer rtrn=new Kmer();
		
		for(int i=0; i<tokens.length; i++){
			rtrn.addRegion(tokens[i]);
		}
		
		return rtrn;
	}
	
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
	
	
	private static void write(String string, Collection<Cluster> clusters) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Cluster c: clusters){
			writer.write(c.getLine()+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>5){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			int resolution=new Integer(args[1]);
			String save=args[2];
			Kmer kmer=parse(args[3]);
			boolean weight=new Boolean(args[4]);
			SingleInterval region=new SingleInterval(args[5]);
			
			
			Collection<Cluster> clusters=data.getClusters(kmer);
			
			System.err.println("Retained clusters "+clusters.size());
			
			write(save+".clusters", clusters);
			
			
			
			MatrixWithHeaders mwh=data.getDNADNAContactMatrix(resolution, weight, clusters, region);
			mwh.write(save);
			data.close();
		}
		else{System.err.println(usage);}
	}
	
	

	



	static String usage=" args[0]=barcodes \n args[1]=resolution \n args[2]=save \n args[3]=gene (multiple with _) \n args[4]=weight \n args[5]=region";
	
}

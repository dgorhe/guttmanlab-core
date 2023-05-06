package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.rnasprite.Kmer;

public class ClusterHeatmap {
	
	int minCount=0;

	public ClusterHeatmap(BarcodingDataStreaming data, Map<String, String> allRNA, Collection<Kmer> hubs, String save) throws IOException {
		List<String> columns=getColumns(allRNA);
		List<String> rows=getRows(data);
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		Map<String, Collection<String>> clusterToGroups=new TreeMap<String, Collection<String>>();
		
		List<String> subRows=new ArrayList<String>();
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Kmer hits=getHits(c, allRNA);
			Collection<String> hubHits=count(hits, hubs);
			if(!hubHits.isEmpty()) {
				for(String rna: hits.getRegions()) {
					mwh.incrementCount(c.getBarcode(), rna);
				}
				subRows.add(c.getBarcode());
				clusterToGroups.put(c.getBarcode(), hubHits);
			}
			
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		mwh=mwh.submatrixByRowNames(subRows);
		data.close();
		mwh.write(save);
		write(save+".annotation", clusterToGroups);
	}
	
	private void write(String save, Map<String, Collection<String>> clusterToGroups) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String name: clusterToGroups.keySet()) {
			Collection<String> set=clusterToGroups.get(name);
			String val="";
			int counter=0;
			for(String s: set) {
				if(counter>0) {val+="_";}
				val+=s;
				counter++;
			}
			writer.write(name+"\t"+val+"\n");
		}
		
		writer.close();
	}

	private int count(Kmer hits, Kmer hub) {
		int counter=0;
		for(String rna: hits.getRegions()) {
			if(hub.containsRegion(rna)) {counter++;}
		}
		return counter;
	}
	
	private Collection<String> count(Kmer hits, Collection<Kmer> hubs) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(Kmer hub: hubs) {
			int count=count(hits, hub);
			if(count>minCount) {rtrn.add(hub.getName());}
		}
		
		return rtrn;
	}

	private List<String> getColumns(Map<String, String> allRNA) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String name: allRNA.keySet()) {
			String className=allRNA.get(name);
			if(!rtrn.contains(className)) {
				rtrn.add(className);
			}
		}
		
		return rtrn;
	}

	private List<String> getRows(BarcodingDataStreaming data) {
		List<String> rtrn=new ArrayList<String>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			rtrn.add(c.getBarcode());
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		data.close();
		
		
		return rtrn;
	}

	private Kmer getHits(Cluster c, Map<String, String> allRNA) {
		Kmer rtrn=new Kmer();
		
		for(String name: c.getRNANames()) {
			if(allRNA.containsKey(name)) {
				String className=allRNA.get(name);
				rtrn.addRegion(className);
			}
		}
		
		return rtrn;
	}

	private boolean containsAny(Cluster c, Map<String, String> allRNA, Kmer hub) {
		for(String rna: c.getRNANames()) {
			if(allRNA.containsKey(rna)) {
				String className=allRNA.get(rna);
				if(hub.containsRegion(className)) {return true;}
			}
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
	
	private static Collection<Kmer> parseHubs(String string) throws IOException {
		List<String> lines= BEDFileIO.loadLines(string);
		
		Map<String, Collection<String>> temp=new TreeMap<String, Collection<String>>();
		for(String line: lines) {
			String hub=line.split("\t")[0];
			String className=line.split("\t")[2];
			Collection<String> set=new TreeSet<String>();
			if(temp.containsKey(hub)) {set=temp.get(hub);}
			set.add(className);
			temp.put(hub, set);
		}
		
		
		Collection<Kmer> rtrn=new ArrayList<Kmer>();
		for(String hub: temp.keySet()) {
			Kmer kmer=new Kmer();
			kmer.setName(hub);
			kmer.addRegions(temp.get(hub));
			rtrn.add(kmer);
		}
		
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Map<String, String> rnas=parseRNAs(args[1]);
			Collection<Kmer> hubs=parseHubs(args[1]);
			String save=args[2];
			new ClusterHeatmap(data, rnas, hubs, save);
		}
		
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=barcodes \n args[1]=all rnas \n args[2]=save";
}

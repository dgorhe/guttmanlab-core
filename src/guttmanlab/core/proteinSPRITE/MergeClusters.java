package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.rnasprite.Cluster;

public class MergeClusters {

	public MergeClusters(BarcodingDataStreaming data, String save) throws IOException {
		Map<String, Integer> counts=getCounts(data);
		write(save+".counts", counts);
		counts=filter(counts, 1);
		
		Map<String, Cluster> merged=mergeClusters(data, counts, save);
		
		//write(save, counts);
		//write(save, merged, counts);
	}

	private Map<String, Integer> filter(Map<String, Integer> counts, int i) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(String name: counts.keySet()) {
			int count=counts.get(name);
			if(count>i) {rtrn.put(name, count);}
		}
		
		return rtrn;
	}

	private Map<String, Integer> getCounts(BarcodingDataStreaming data) {
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			int count=0;
			if(counts.containsKey(c.getBarcode())) {
				count=counts.get(c.getBarcode());
			}
			count++;
			counts.put(c.getBarcode(), count);
			counter++;
			if(counter%100000 ==0) {System.err.println(counter);}
		}
		data.close();
		return counts;
	}

	private Map<String, Cluster> mergeClusters(BarcodingDataStreaming data, Map<String, Integer> counts, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		Map<String, Cluster> merged=new TreeMap<String, Cluster>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			String barcode=c.getBarcode();
			if(counts.containsKey(barcode)) {
				if(merged.containsKey(barcode)) {
					Cluster c1=merged.get(barcode);
					c=merge(c, c1);
				}
				merged.put(barcode, c);
			}
			else {
				writer.write(c.toString()+"\n");
			}
			
			counter++;
			if(counter%100000 ==0) {System.err.println(counter);}
		}
		data.close();
		
		write(writer, merged);
		
		writer.close();
		
		return merged;
	}

	private Cluster merge(Cluster c, Cluster c1) {
		Cluster rtrn=new Cluster(c.getBarcode());
		
		rtrn.addDNAReads(c.getAllDNAIntervals());
		rtrn.addDNAReads(c1.getAllDNAIntervals());
		
		rtrn.addProteins(c.getProteins());
		rtrn.addProteins(c1.getProteins());
		
		return rtrn;
	}

	private void write(String save, Map<String, Cluster> merged, Map<String, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String name: merged.keySet()) {
			Cluster c=merged.get(name);
			int count=counts.get(name);
			if(count>1) {
				writer.write(c.toString()+"\n");
			}
		}
		
		writer.close();
	}
	
	private void write(FileWriter writer, Map<String, Cluster> merged) throws IOException {
		for(String name: merged.keySet()) {
			Cluster c=merged.get(name);
			writer.write(c.toString()+"\n");
		}
	}
	
	private void write(String save, Map<String, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String name: counts.keySet()) {
			//Cluster c=merged.get(name);
			int count=counts.get(name);
			writer.write(name+"\t"+count+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		new MergeClusters(data, save);
	}
	
}

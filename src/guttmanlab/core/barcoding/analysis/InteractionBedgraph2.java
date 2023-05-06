package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class InteractionBedgraph2 {
	int binResolution=100;

	public InteractionBedgraph2(BarcodingDataStreaming data, String rna, String save) throws IOException {
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		Map<String, Integer> allCounts=new TreeMap<String, Integer>();
		while(data.hasNext()){
			Cluster c=data.next();
			int size=c.getClusterSize();
			if(size>1 && c.containsChr(rna)){
				for(SingleInterval region: c.getAllIntervals()){
					int count=0;
					String chr=region.getReferenceName();
					if(counts.containsKey(chr)){count=counts.get(chr);}
					count++;
					counts.put(chr, count);
				}
			}
			
			for(SingleInterval region: c.getAllIntervals()){
				int count=0;
				String chr=region.getReferenceName();
				if(allCounts.containsKey(chr)){count=allCounts.get(chr);}
				count++;
				allCounts.put(chr, count);
			}
			
		}
		write(save, counts, allCounts);
	}

	
	
	public InteractionBedgraph2(BarcodingDataStreaming data, String rna1, String rna2, String save) throws IOException {
		
		//int length1=data.getReferenceLength(rna1);
		//int length2=data.getReferenceLength(rna2);
		
		List<String> rows=makeRows(data, rna1);
		List<String> columns=makeRows(data, rna2);
		
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		while(data.hasNext()){
			Cluster c=data.next();
			c=c.bin(binResolution);
			for(SingleInterval region1: c.getAllIntervals()) {
				for(SingleInterval region2: c.getAllIntervals()) {
					if(!region1.equals(region2) && mwh.containsRow(region1.toUCSC()) && mwh.containsColumn(region2.toUCSC())) {
						mwh.incrementCount(region1.toUCSC(), region2.toUCSC());
					}
				}
			}
			
			
			
		}
		
		data.close();
		
		mwh.write(save);
	}

	

	private List<String> makeRows(BarcodingDataStreaming data, String rna1) {
		Collection<SingleInterval> list= new TreeSet<SingleInterval>();
		while(data.hasNext()) {
			Cluster c=data.next();
			c=c.bin(binResolution);
			for(SingleInterval region: c.getAllIntervals()) {
				if(region.getReferenceName().equalsIgnoreCase(rna1)) {
					list.add(region);
				}
			}
		}
		data.close();
		
		List<String> rtrn=new ArrayList<String>();
		
		for(SingleInterval region: list) {rtrn.add(region.toUCSC());}
		
		return rtrn;
	}



	private List<String> makeRows(String rna1, int length1) {
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<length1; i++) {
			String name=rna1+":"+i+"-"+(i+1);
			rtrn.add(name);
		}
		
		return rtrn;
	}



	private void write(String save, Map<String, Integer> counts, Map<String, Integer> allCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String rna: counts.keySet()){
			int count=counts.get(rna);
			int allCount=allCounts.get(rna);
			writer.write(rna+"\t"+count+"\t"+allCount+"\n");
		}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		if(args.length>2){
			/*BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String rna1=args[1];
			String rna2=args[2];
			String save=args[3];
			new InteractionBedgraph2(data, rna1, rna2, save);*/
			
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String rna1=args[1];
			String save=args[2];
			new InteractionBedgraph2(data, rna1, save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=cluster file \n args[1]=RNA 1 \n args[2]=save name";
}

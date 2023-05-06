package guttmanlab.core.chipdip;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class SingleCellAnalysis {
	
	int binResolution=10;

	public SingleCellAnalysis(File file, String save) throws IOException {
		
		Map<String, Collection<String>> cellToRegions=singleCell(file, save);
		
		//write(save, cellToRegions, minCount);
		Map<SingleInterval, Integer> counts=ensemble(cellToRegions);
		BEDFileIO.writeBEDGraphInteger(counts, save);
	}
	
	

	
	private Map<SingleInterval, Integer> ensemble(Map<String, Collection<String>> cellToRegions) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
	
		for(String cell: cellToRegions.keySet()) {
			for(String r: cellToRegions.get(cell)) {
				SingleInterval region=new SingleInterval(r);
				Collection<SingleInterval> bins=region.allBins(binResolution);
				for(SingleInterval bin: bins) {
					int count=0;
					if(rtrn.containsKey(bin)) {
						count=rtrn.get(bin);
					}
					count++;
					rtrn.put(bin, count);
				}
			}
		}
		
		return rtrn;
	}




	private Map<String, Collection<String>> singleCell(File file, String save) throws IOException {
		Map<String, Collection<String>> cellToRegions=new TreeMap<String, Collection<String>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		int counter=0;
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			parse(nextLine, cellToRegions);
			counter++;
			if(counter%10000 ==0) {System.err.println(counter);}
			
		}
		reader.close();
		System.err.println("done reading file");	
		
		return cellToRegions;
		
		
		
	}


	private void write(String save, Map<String, Collection<String>> cellToRegions, int minCount) throws IOException {
		cellToRegions=filter(cellToRegions, minCount);
		
		
		int counter=0;
		for(String barcode: cellToRegions.keySet()) {
			FileWriter writer=new FileWriter(save+"/"+barcode+".bedgraph");
			Collection<String> list=cellToRegions.get(barcode);
			Map<SingleInterval, Integer> counts=count(list);
			for(SingleInterval region: counts.keySet()) {
				//SingleInterval r=new SingleInterval(region).bin(binResolution);
				int count=counts.get(region);
				writer.write(region.toBedgraph(count)+"\n");
			}
			writer.close();
			counter++;
			if(counter%1000==0) {System.err.println(barcode+" "+counter+" "+cellToRegions.size());}
		}
		
	
	}
	
	private Map<SingleInterval, Integer> count(Collection<String> list) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(String region: list) {
			SingleInterval r=new SingleInterval(region).bin(binResolution);
			int count=0;
			if(rtrn.containsKey(r)) {count=rtrn.get(r);}
			count++;
			rtrn.put(r, count);
		}
		
		return rtrn;
	}


	/*private void write(String save, Map<String, Collection<String>> cellToRegions, int minCount) throws IOException {
		cellToRegions=filter(cellToRegions, minCount);
		Collection<SingleInterval> allBins=getAllBins(cellToRegions);
		List<String> columns=toBins(allBins);
		List<String> rows=new ArrayList<String>();
		rows.addAll(cellToRegions.keySet());
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		int counter=0;
		for(String barcode: cellToRegions.keySet()) {
			Collection<String> list=cellToRegions.get(barcode);
			for(String region: list) {
				String column=new SingleInterval(region).bin(binResolution).toUCSC();
				mwh.set(barcode, column, 1.0);
			}
			counter++;
			if(counter%1000==0) {System.err.println(barcode+" "+counter+" "+cellToRegions.size());}
		}
		
		mwh.write(save);
	}*/

	private Map<String, Collection<String>> filter(Map<String, Collection<String>> cellToRegions, int minCount) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String cell: cellToRegions.keySet()) {
			Collection<String> list=cellToRegions.get(cell);
			if(list.size()>minCount) {rtrn.put(cell, list);}
		}
		
		return rtrn;
	}

	private List<String> toBins(Collection<SingleInterval> allBins) {
		List<String> rtrn=new ArrayList<String>();
		
		for(SingleInterval r: allBins) {
			rtrn.add(r.toUCSC());
		}
		
		return rtrn;
	}

	private Collection<SingleInterval> getAllBins(Map<String, Collection<String>> cellToRegions) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(String s: cellToRegions.keySet()) {
			Collection<String> list=cellToRegions.get(s);
			for(String val: list) {
				SingleInterval r=new SingleInterval(val).bin(binResolution);
				rtrn.add(r);
			}
		}
		
		return rtrn;
	}

	private String toList(Collection<SingleInterval> list) {
		String rtrn="";
		
		for(SingleInterval i:list) {
			rtrn+=i.toUCSC()+",";
		}
		
		if(rtrn.endsWith(",")) {
			rtrn=rtrn.substring(0, rtrn.length()-1);
		}
		
		return rtrn;
	}

	private void parse(String nextLine, Map<String, Collection<String>> cellToRegions) {
		String[] tokens=nextLine.split("\t");
		
		for(int i=0; i<tokens.length; i++) {
			if(tokens[i].startsWith("DPM")) {
				update(tokens[i], cellToRegions);
			}
		}
		
	}

	private void update(String string, Map<String, Collection<String>> cellToRegions) {
		String[] tokens=string.split("]_");
		String region=tokens[1];
		
		String scBarcode=tokens[0].split("\\[")[1];
		
		if(!scBarcode.contains("FOUND")) {
			Collection<String> list=new TreeSet<String>();
			if(cellToRegions.containsKey(scBarcode)) {
				list=cellToRegions.get(scBarcode);
			}
			list.add(region);
			cellToRegions.put(scBarcode, list);
		}
	}
	
	public static void main(String[] args) throws IOException {
		File file=new File(args[0]);
		String saveDir=args[1];
		new SingleCellAnalysis(file, saveDir);
	}
	
}

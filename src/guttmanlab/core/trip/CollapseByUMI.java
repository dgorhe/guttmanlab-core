package guttmanlab.core.trip;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.nanobody.Cluster;
import guttmanlab.core.sequence.Sequence;

public class CollapseByUMI {
	
	public CollapseByUMI(String input, Collection<String> gDNABarcodes, int minDistance, String save) throws IOException {
		Map<String, Collection<String>> rawCounts=parse(input);
		
		Map<String, Collection<String>> barcodes=getClose(gDNABarcodes, rawCounts, minDistance);
		
		Map<String, Collection<Cluster>> collapsedUMIs=collapseUMIs(barcodes, minDistance);
			
		write(save, collapsedUMIs);
	}
	
	public CollapseByUMI(String input, int minDistance, String save) throws IOException {
		Map<String, Collection<String>> rawCounts=parse(input);
		
		Map<String, Collection<String>> barcodes=getClose(rawCounts, minDistance);
		
		Map<String, Collection<Cluster>> collapsedUMIs=collapseUMIs(barcodes, minDistance);
			
		write(save, collapsedUMIs);
	}
	
	
	private Map<String, Collection<Cluster>> collapseUMIs(Map<String, Collection<String>> barcodes, int minDistance) {
		Map<String, Collection<Cluster>> rtrn=new TreeMap<String, Collection<Cluster>>();
		for(String barcode: barcodes.keySet()) {
			Collection<Cluster> clusters=Sequence.cluster(barcodes.get(barcode), minDistance);
			rtrn.put(barcode, clusters);
		}
		
		
		return rtrn;
	}

	
	private Map<String, Collection<String>> getClose(Map<String, Collection<String>> rawCounts, int minDistance) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		
		for(String gDNABarcode: rawCounts.keySet()) {
			Collection<String> close=getClose(gDNABarcode, rawCounts, minDistance);
			Collection<String> list=new TreeSet<String>();
			for(String c: close) {
				list.addAll(rawCounts.get(c));
			}
			System.out.println(gDNABarcode+" "+list.size()+" "+close.size());
			rtrn.put(gDNABarcode, list);
			
		}
		
		
		return rtrn;
	}
	

	private Map<String, Collection<String>> getClose(Collection<String> gDNABarcodes, Map<String, Collection<String>> rawCounts, int minDistance) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		
		for(String gDNABarcode: gDNABarcodes) {
			Collection<String> close=getClose(gDNABarcode, rawCounts, minDistance);
			Collection<String> list=new TreeSet<String>();
			for(String c: close) {
				list.addAll(rawCounts.get(c));
			}
			System.out.println(gDNABarcode+" "+list.size()+" "+close.size());
			rtrn.put(gDNABarcode, list);
			
		}
		
		
		return rtrn;
	}


	private Collection<String> getClose(String gDNABarcode, Map<String, Collection<String>> rawCounts, int minDistance) {
		Collection<String> rtrn=new TreeSet<String>();
		for(String s2: rawCounts.keySet()) {
			if(distance(gDNABarcode, s2, minDistance)) {rtrn.add(s2);}
		}
		return rtrn;
	}


	private void write(String save, Map<String, Collection<Cluster>> collapsedUMIs) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: collapsedUMIs.keySet()) {
			Collection<Cluster> set=collapsedUMIs.get(barcode);
			writer.write(barcode+"\t"+set.size());
			for(Cluster s: set) {writer.write("\t"+s.toString());}
			writer.write("\n");
		}
		
		writer.close();
	}


	private static Map<String, Collection<String>> unique(String input, int minDistance) throws IOException {
		Map<String, Collection<String>> counts=parse(input);
		
		Map<String, Collection<String>> barcodes=collapseList(counts.keySet(), minDistance, counts);
		//for(String b: barcodes) {System.out.println(b);}
		
		Map<String, Collection<String>> collapsed=collapse(barcodes, minDistance);
		
		return collapsed;
	}



	
	private static Map<String, Collection<String>> collapse(Map<String, Collection<String>> counts, int minDistance) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		int count=0;
		for(String barcode: counts.keySet()) {
			System.err.println(count+" "+barcode);
			Collection<String> list=counts.get(barcode);
			Map<String, Collection<String>> collapsed=collapse(list, minDistance);
			Collection<String> collapsedList=list(collapsed, minDistance);
			
			//System.out.println(barcode+"\t"+list.size()+"\t"+collapsedList.size());
			
			
			rtrn.put(barcode, collapsedList);
			count++;
		}
		
		return rtrn;
	}

	
	private static Map<String, Collection<String>> collapseList(Collection<String> list, int minDistance, Map<String, Collection<String>> input) {
		Map<String, Collection<String>> collapsed=collapse(list, minDistance);
		Collection<String> collapsedList=list2(collapsed);
		
		
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String c: collapsedList) {
			Collection<String> set=getOverlaps(c, input);
			rtrn.put(c, set);
		}
		
		return rtrn;
	}
	



	private static Collection<String> getOverlaps(String c, Map<String, Collection<String>> input) {
		Collection<String> list=parseCollection(c);
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String s: list) {
			rtrn.addAll(input.get(s));
		}
		
		return rtrn;
	}


	
	private static Collection<String> list2(Map<String, Collection<String>> collapsed) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String barcode: collapsed.keySet()) {
			/*Collection<String> full=new TreeSet<String>();
			for(String other: collapsed.get(barcode)) {
				full.addAll(collapsed.get(other));
			}*/
			
			String string=toString(collapsed.get(barcode));
			rtrn.add(string);
		}
		
		//rtrn=removeSubset(rtrn);
		
		return rtrn;
	}


	private static Collection<String> list(Map<String, Collection<String>> collapsed, int minDistance) {
		
		//TODO if close is 1 then add
		Collection<String> rtrn=new TreeSet<String>();
	
		
		
		
		
		
		for(String barcode: collapsed.keySet()) {
			Map<String, Collection<String>> temp=new TreeMap<String, Collection<String>>();
			
			Collection<String> full=new TreeSet<String>();
			for(String other: collapsed.get(barcode)) {
				full.addAll(collapsed.get(other));
			}
			
			for(String f: full) {
				if(!temp.containsKey(f)) {temp.put(f, new TreeSet<String>());}
				Collection<String> set=temp.get(f);
				set.addAll(full);
			}
			
			for(String f: temp.keySet()) {
				String string=toString(temp.get(f));
				rtrn.add(string);
			}
		}
		
		rtrn=removeSubset(rtrn);
		
		return rtrn;
	}




	private static boolean isClose(Collection<String> collection, Collection<String> collection2, int minDistance) {
		for(String c1: collection) {
			for(String c2: collection2) {
				if(distance(c1, c2, minDistance)) {return true;}
			}
		}
		return false;
	}




	private static Collection<String> removeSubset(Collection<String> list) {
		Collection<String> remove=new TreeSet<String>();
		
		for(String str1: list) {
			for(String str2: list) {
				if(subset(str2, str1)) {remove.add(str2);}
			}
		}
		
		Collection<String> rtrn=new TreeSet<String>();
		for(String r: list) {
			if(!remove.contains(r)) {rtrn.add(r);}
		}
		return rtrn;
	}




	private static boolean subset(String str2, String str1) {
		Collection<String> tokens2=parseCollection(str2);
		Collection<String> tokens1=parseCollection(str1);
	
		if(tokens2.size()>=tokens1.size()) {return false;}
		
		for(String val2: tokens2) {
			if(!tokens1.contains(val2)) {return false;}
		}
		
		return true;
	}




	private static Collection<String> parseCollection(String str2) {
		Collection<String> rtrn=new TreeSet<String>();
		
		String[] tokens=str2.split("_");
		for(int i=0; i<tokens.length; i++) {
			rtrn.add(tokens[i]);
		}
		
		return rtrn;
	}




	private static String toString(Collection<String> collection) {
		String rtrn="";
		for(String r: collection) {
			rtrn+=r+"_";
		}
		return rtrn.substring(0, rtrn.length()-1);
	}




	private static Map<String, Collection<String>> collapse(Collection<String> list, int minDistance) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		for(String string1: list) {
			Collection<String> close=new TreeSet<String>();
			close.add(string1);
			for(String string2: list) {
				if(distance(string1, string2, minDistance)) {
					close.add(string2);
				}
			}
			rtrn.put(string1, close);
		}
		
		return rtrn;
	}




	private static boolean distance(String seq1, String seq2, int minDistance) {
		char[] char1=seq1.toCharArray();
		char[] char2=seq2.toCharArray();
		
		int distance=0;
		for(int i=0; i<char1.length; i++) {
			if(char1[i]!=char2[i]) {
				distance++;
				if(distance>minDistance) {return false;}
			}
		}
		return true;
	}




	private static Map<String, Collection<String>> parse(String input) throws IOException {
		List<String> lines=BEDFileIO.loadLines(input);
		Map<String, Collection<String>> map=new TreeMap<String, Collection<String>>();
		for(String line: lines) {
			String[] tokens=line.split(" ");
			//System.err.println(line+" "+tokens.length+" "+tokens[tokens.length-1]);
			String cellBarcode=tokens[tokens.length-1].split("_")[0];
			String umi=tokens[tokens.length-1].split("_")[1];
			if(!map.containsKey(cellBarcode)) {map.put(cellBarcode, new TreeSet<String>());}
			Collection<String> list=map.get(cellBarcode);
			list.add(umi);
		}
		return map;
	}




	private static void write(Map<String, Collection<String>> counts, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: counts.keySet()) {
			Collection<String> list=counts.get(barcode);
			writer.write(barcode);
			for(String val: list) {
				writer.write("\t"+val);
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			String input=args[0];
			//Collection<String> gDNA=BEDFileIO.loadLines(args[1]);
			String save=args[1];
			int minDistance=Integer.parseInt(args[2]);
			new CollapseByUMI(input, minDistance, save);
		}
		
		
		
		
		else {System.err.println(usage);}
	}
	
	

	
	static String usage=" args[0]=input \n args[1]=save \n args[2]=min distance";
	
}

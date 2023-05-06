package guttmanlab.core.nanobody;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.sequence.Sequence;

public class Utils {

	
	static Collection<String> collapseBarcodes(Collection<String> barcodes, int distance){
		Map<String, Collection<String>> nearest=findClosest(barcodes, distance);
		Collection<String> sets=makeSets(nearest);
		return sets;
	}
	
	
	private static Map<String, Collection<String>> findClosest(Collection<String> barcodes, int n) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		for(String barcode1: barcodes) {
			Collection<String> list=new TreeSet<String>();
			for(String barcode2: barcodes) {
				int distance=Sequence.distance(barcode1, barcode2);
				if(distance<=n) {
					list.add(barcode2);
				}
			}
			rtrn.put(barcode1, list);
		}
		return rtrn;
	}
	
	private static Collection<String> makeSets(Map<String, Collection<String>> nearest) {
		Collection<String> rtrn=new TreeSet<String>();
		
		Collection<String> visited2=new TreeSet<String>();
		
		for(String barcode: nearest.keySet()) {
			if(!visited2.contains(barcode)) {
				Collection<String> visited=new TreeSet<String>();
				Collection<String> set=new TreeSet<String>();
				set.addAll(nearest.get(barcode));
				
				Collection<String> diff= difference(set, visited);
				
				while(diff.size()>0) {
					for(String barcode2: diff) {
						visited2.add(barcode2);
						set.addAll(nearest.get(barcode2));
						visited.add(barcode2);
					}
					diff=difference(set, visited);
				}
				
				String setString=toString(set);
				rtrn.add(setString);
			}
		}
			
			
		return rtrn;
	}
	
	private static String toString(Collection<String> set) {
		String rtrn="";
		for(String s: set) {
			rtrn+=s+" ";
		}
		return rtrn;
	}
	
	private static Collection<String> difference(Collection<String> set, Collection<String> visited) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String val: set) {
			if(!visited.contains(val)) {rtrn.add(val);}
		}
		
		return rtrn;
	}
	
}

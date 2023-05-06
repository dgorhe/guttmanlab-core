package guttmanlab.core.rnasprite;

import java.io.BufferedReader;
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

import org.apache.commons.math3.util.CombinatoricsUtils;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.math.Statistics;

public class BarcodeString implements Comparable<BarcodeString>{

	//int numberOfTags;
	private ArrayList<String> tags;
	private String barcodeString;
	
	public BarcodeString(String barcode) {
		this.barcodeString=barcode;
		String[] tokens=barcode.split("\\.");
		this.tags=new ArrayList<String>();
		for(int i=0; i<tokens.length; i++) {tags.add(tokens[i]);}
	}
	
	public int numberOfTags() {return tags.size();}
	
	public List<String> getTags(){return tags;}
	
	public boolean equals(BarcodeString other) {
		if(this.numberOfTags()!=other.numberOfTags()) {return false;}
			for(int i=0; i<getTags().size(); i++) {
				if(!getTags().get(i).equals(other.getTags().get(i))) {return false;}
			}
		return true;
	}
	
	public int distance(BarcodeString other) {
		int score=0;
		
		for(int i=0; i<getTags().size(); i++) {
			if(!getTags().get(i).equals(other.getTags().get(i))) {score++;}
		}
		
		return score;
	}
	
	@Override
	public int compareTo(BarcodeString o) {
		return this.barcodeString.compareTo(o.barcodeString);
	}

	
	
	private static Map<BarcodeString, Collection<String>> parse(String file) throws IOException {
		Map<BarcodeString, Collection<String>> rtrn=new TreeMap<BarcodeString, Collection<String>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			BarcodeString b=new BarcodeString(tokens[0]);
			Collection<String> set=new TreeSet<String>();
			for(int i=1; i<tokens.length; i++) {set.add(tokens[i]);}
			rtrn.put(b, set);
		}
		reader.close();
		return rtrn;
	}
	
	public String toString() {return this.barcodeString;}
	
	
	private static List<Integer> getMin(BarcodeString b1, Map<BarcodeString, Collection<String>> map) {
		List<Integer> rtrn=new ArrayList<Integer>();
		for(BarcodeString b2: map.keySet()) {
			if(!b1.equals(b2)) {
				if(b1.distance(b2)==1) {
					System.out.println(b1+" "+map.get(b1).size());
					System.out.println(b2+" "+map.get(b2).size());
					System.out.println();
				}
				rtrn.add(b1.distance(b2));
			}
		}
		return rtrn;
	}
	
	private static Map<BarcodeString, Collection<String>> filter(Map<BarcodeString, Collection<String>> map) {
		Map<BarcodeString, Collection<String>> rtrn=new TreeMap<BarcodeString, Collection<String>>();
		
		for(BarcodeString b1: map.keySet()) {
			if(map.get(b1).size()>1) {rtrn.put(b1, map.get(b1));}
		}
		
		return rtrn;
	}

	private static Map<BarcodeString, Collection<String>> merge(Map<BarcodeString, Collection<String>> map) {
		Map<BarcodeString, Collection<String>> rtrn=new TreeMap<BarcodeString, Collection<String>>();
		
		for(BarcodeString bc: map.keySet()) {
			Collection<BarcodeString> rs=bc.remove();
			//System.err.println(bc+" "+r1);
			for(BarcodeString r1: rs) {
				Collection<String> list=new ArrayList<String>();
				if(rtrn.containsKey(r1)) {
					list=rtrn.get(r1);
				}
				list.addAll(map.get(bc));
				rtrn.put(r1, list);
			}
		}
		
		return rtrn;
	}
	
	
	private Collection<BarcodeString> remove() {
		Collection<BarcodeString> rtrn=new ArrayList<BarcodeString>();
		
		for(int i=0; i<numberOfTags(); i++) {
			rtrn.add(remove(i));
		}
		
		return rtrn;
	}

	private BarcodeString remove(int index) {
		String rtrn="";
		
		for(int i=0; i<this.getTags().size(); i++) {
			if(i!=index) {
				rtrn=rtrn+this.getTags().get(i);
				if(i!=this.getTags().size()-1) {rtrn+=".";}
			}
		}
		
		return new BarcodeString(rtrn);
	}

	private static void write(String save, Map<BarcodeString, Collection<String>> merged, Map<BarcodeString, Collection<String>> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(BarcodeString b1: map.keySet()) {
			Collection<BarcodeString> sub=b1.remove();
			
			Collection<String> set=new TreeSet<String>();
			for(BarcodeString s1: sub) {
				set.addAll(merged.get(s1));
			}
			
			writer.write(b1.barcodeString);
			
			for(String l: set) {writer.write("\t"+l);}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		Map<BarcodeString, Collection<String>> map=parse(args[0]);
		
		
		Map<BarcodeString, Collection<String>> merged=merge(map);
		
		write(args[1], merged, map);
		
		/*Map<BarcodeString, Collection<String>> filteredMap=filter(map);
		
		System.err.println(map.size()+" "+filteredMap.size());
		
		int counter=1000;
		int num=0;
		for(BarcodeString b1: filteredMap.keySet()) {
			if(num<counter) {	
				List<Integer> distance=getMin(b1, filteredMap);
				System.out.println(num+" "+b1.barcodeString+" "+map.get(b1).size()+" "+Statistics.min(distance));
			}
			num++;
		}*/
		
		
	}

	

	

	

	

	
	
	
	
}

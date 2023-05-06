package guttmanlab.core.trip;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.sequence.Sequence;

public class CollapseCDNA {
	
	/*private Map<String, Double> getFractionByBarcode(List<String> lines){
		Map<String, Pair<Double>> counts=new TreeMap<String, Pair<Double>>();
		
		for(String line: lines) {
			String[] tokens=line.split(" ");
			String barcode=tokens[0];
			String UMI=tokens[1];
			String spliced=tokens[3];
			
			
			
		}
		
		
	}*/
	
	private static Map<String, Pair<Double>> fractionSpliced(Map<String, Map<String, String>> map) {
		Map<String, Pair<Double>> rtrn=new TreeMap<String, Pair<Double>>();
		
		
		for(String barcode: map.keySet()) {
			Map<String, String> list=map.get(barcode);
			Pair<Double> fraction=fractionSpliced2(list);
			rtrn.put(barcode, fraction);
		}
		
		return rtrn;
	}
	
	private static Pair<Double> fractionSpliced2(Map<String, String> list) {
		double total=0;
		double spliced=0;
		for(String umi: list.keySet()) {
			String type=list.get(umi);
			if(type.equals("spliced")) {spliced++;}
			total++;
		}
		
		Pair<Double> rtrn=new Pair<Double>();
		rtrn.setValue1(spliced);
		rtrn.setValue2(total);
		
		return rtrn;
		
		
	}

	private static Map<String, Pair<Double>> fractionSpliced(Collection<String> sets, Map<String, Collection<String>> map) {
		Map<String, Pair<Double>> rtrn=new TreeMap<String, Pair<Double>>();
		
		
		for(String set: sets) {
			String[] tokens=set.split(" ");
			//Pair<Double> pair=new Pair<Double>(0.0,0.0);
			Collection<String> allUmis=new TreeSet<String>();
			for(int i=0; i<tokens.length; i++) {
				String barcode=tokens[i];
				Collection<String> list=map.get(barcode);
				/*Pair<Double> fraction=fractionSpliced(list);
				pair.setValue1(pair.getValue1()+fraction.getValue1());
				pair.setValue2(pair.getValue2()+fraction.getValue2());*/
				allUmis.addAll(list);
			}
			
			for(String umi: allUmis) {
				System.out.println(set+"\t"+umi);
			}
			
			Pair<Double> pair=fractionSpliced(allUmis);
			rtrn.put(set, pair);
		}
		
		return rtrn;
	}
	
	
	private static Pair<Double> fractionSpliced(Collection<String> list) {
		double total=0;
		double spliced=0;
		for(String umi: list) {
			String[] tokens=umi.split(" ");
			if(tokens[1].equals("spliced")) {spliced++;}
			total++;
		}
		
		Pair<Double> rtrn=new Pair<Double>();
		rtrn.setValue1(spliced);
		rtrn.setValue2(total);
		
		return rtrn;
	}

	
	private static Map<String, Map<String, Integer>> countByUMI(List<String> lines) {
		Map<String, Map<String, Integer>> barcodes=new TreeMap<String, Map<String, Integer>>();
		
		int counter=0;
		for(String line: lines) {
			String[] tokens=line.split(" ");
			String barcode=tokens[0];
			String UMI=tokens[2];
			String spliced=tokens[3];
			if(!barcodes.containsKey(barcode)) {barcodes.put(barcode, new TreeMap<String, Integer>());}
			Map<String, Integer> set=barcodes.get(barcode);
			int count=0;
			if(set.containsKey(UMI)) {count=set.get(UMI);}
			count++;
			set.put(UMI, count);
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+lines.size());}
		}
		return barcodes;
	}

	static Map<String, Collection<String>> collapseByUMI(List<String> lines) {
		Map<String, Collection<String>> barcodes=new TreeMap<String, Collection<String>>();
		
		int counter=0;
		for(String line: lines) {
			String[] tokens=line.split(" ");
			String barcode=tokens[0];
			String UMI=tokens[2];
			String spliced=tokens[3];
			if(!barcodes.containsKey(barcode)) {barcodes.put(barcode, new TreeSet<String>());}
			Collection<String> set=barcodes.get(barcode);
			set.add(UMI+" "+spliced);
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+lines.size());}
		}
		return barcodes;
	}
	
	private static void write(Map<String, Pair<Double>> map, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: map.keySet()) {
			Pair<Double> pair=map.get(barcode);
			double fraction=pair.getValue1()/pair.getValue2();
			writer.write(barcode+"\t"+pair.getValue1()+"\t"+pair.getValue2()+"\t"+fraction+"\n");
		}
		writer.close();
	}
	
	
	private static void writeUMI(Map<String, Collection<String>> map) {
		for(String barcode: map.keySet()) {
			Collection<String> list=map.get(barcode);
			for(String umi: list) {
				System.out.println(barcode+" "+umi);
			}
		}
		
	}
	
	

	private static Collection<String> makeSets(Map<String, Collection<String>> nearest) {
		Collection<String> rtrn=new TreeSet<String>();
		
		Collection<String> visited2=new TreeSet<String>();
		
		int counter=0;
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
					//System.err.println(barcode+" "+barcode2+" "+nearest.get(barcode2).size()+" "+set.size()+" "+visited.size());
				}
				diff=difference(set, visited);
			}
			
			String setString=toString(set);
			rtrn.add(setString);
			counter++;
			
			//if(counter%100 ==0) {System.err.println(counter+" "+nearest.size()+" "+barcode+" "+nearest.get(barcode).size()+" "+set.size());}
			//if(counter==1) {return null;}
			}
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

	private static void add(Map<String, Collection<String>> rtrn, String key, Collection<String> list) {
		if(!rtrn.containsKey(key)) {rtrn.put(key, new TreeSet<String>());}
		Collection<String> temp=rtrn.get(key);
		for(String val: list) {temp.add(val);}
	}

	private static String toString(Collection<String> set) {
		String rtrn="";
		for(String s: set) {
			rtrn+=s+" ";
		}
		return rtrn;
	}


	private static Map<String, Collection<String>> findClosest(Collection<String> barcodes, int n) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		int counter=0;
		for(String barcode1: barcodes) {
			Collection<String> list=new TreeSet<String>();
			for(String barcode2: barcodes) {
				int distance=Sequence.distance(barcode1, barcode2);
				if(distance<=n) {
					//if(distance>0) {System.err.println(barcode1+" "+barcode2+" "+distance);}
					list.add(barcode2);
				}
			}
			rtrn.put(barcode1, list);
			counter++;
			//if(counter%1000==0) {System.err.println(counter+" "+barcodes.size());}
		}
		return rtrn;
	}
	
	private static MatrixWithHeaders makeDistanceMatrix(Map<String, Collection<String>> map) {
		List<String> rows=new ArrayList<String>();
		for(String b: map.keySet()) {
			rows.add(b);
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, rows);
		
		int counter=0;
		for(String barcode1: map.keySet()) {
			for(String barcode2: map.keySet()) {
				int distance=Sequence.distance(barcode2, barcode1);
				rtrn.set(barcode1, barcode2, distance);
			}
			counter++;
			if(counter%10000 ==0) {System.err.println(counter+" "+map.size());}
		}
		
		return rtrn;
	}
	
	
	public static Collection<String> getCollapsedSets(Set<String> barcodes, int distance){
		Map<String, Collection<String>> nearest=findClosest(barcodes, distance);
		
		Collection<String> sets=makeSets(nearest);
		return sets;
	}
	
	


	/*private static Map<String, Collection<String>> collapseUMIs(Collection<String> sets, Map<String, Collection<String>> map) {
		Map<String, Map<String, String>> rtrn=new TreeMap<String, Map<String, String>>();
		
		for(String set: sets) {
			String[] tokens=set.split(" ");
			//Pair<Double> pair=new Pair<Double>(0.0,0.0);
			Map<String, String> allUmis=new TreeMap<String, String>();
			for(int i=0; i<tokens.length; i++) {
				String barcode=tokens[i];
				Collection<String> list=map.get(barcode);
				
				//allUmis.addAll(list);
				for(String val: list) {
					allUmis.put(val.split(" ")[0], val.split(" ")[1]);
				}
				
			}
			
			System.out.println(set+"\t"+allUmis.size());
			
			
			//Pair<Double> pair=fractionSpliced(allUmis);
			rtrn.put(set, allUmis);
		}
		
		Map<String, Collection<String>> collapsedMap=new TreeMap<String, Collection<String>>();
		int counter=0;
		for(String barcode: rtrn.keySet()) {
			int startingSize=rtrn.get(barcode).size();
			Map<String, Collection<String>> nearest=findClosest(rtrn.get(barcode).keySet(), 2);
			Collection<String> collapsed=makeSets(nearest);
			System.out.println(barcode+"\t"+collapsed.size());
			collapsedMap.put(barcode, collapsed);
			counter++;
			System.err.println(counter+" "+rtrn.size()+" "+barcode +" "+startingSize+" "+collapsed.size());
		}
		
		return collapsedMap;
	}*/


	static Map<String, Collection<String>> collapseUMIs(Map<String, Map<String, String>> map, int distance) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String barcode: map.keySet()) {
			Collection<String> collapsed=collapseUMI(map.get(barcode), distance);
			rtrn.put(barcode, collapsed);
			System.out.println(barcode+"\t"+map.get(barcode).size()+"\t"+collapsed.size());
		}
		
		return rtrn;
	}
	
	static Collection<String> collapseUMI(Map<String, String> umis, int distance) {
		Map<String, Collection<String>> nearest=findClosest(umis.keySet(), distance);
		Collection<String> collapsed=makeSets(nearest);
		return collapsed;
		
	}
	
	static Collection<String> collapseUMI(Collection<String> umis, int distance) {
		Map<String, Collection<String>> nearest=findClosest(umis, distance);
		Collection<String> collapsed=makeSets(nearest);
		return collapsed;
		
	}

	static Map<String, Map<String, String>> collapseBarcodes(List<String> lines, int distance) {
		Map<String, Map<String, String>> rtrn=new TreeMap<String, Map<String, String>>();
		Map<String, Collection<String>> map=collapseByUMI(lines);
		
		//System.err.println("num "+map.size());
		
		Map<String, Collection<String>> nearest=findClosest(map.keySet(), distance);
		
		Collection<String> sets=makeSets(nearest);
		
		//TODO: Keep UMI count
		for(String set: sets) {
			String[] tokens=set.split(" ");
			Map<String, String> allUmis=new TreeMap<String, String>();
			for(int i=0; i<tokens.length; i++) {
				String barcode=tokens[i];
				Collection<String> list=map.get(barcode);
				for(String val: list) {
					allUmis.put(val.split(" ")[0], val.split(" ")[1]);
				}
				
			}
			
			//System.out.println(set+"\t"+allUmis.size());
			
			rtrn.put(set, allUmis);
		}
		
		return rtrn;
	}

	
	static Collection<String> collapseBarcodes(Collection<String> barcodes, int distance){
		Map<String, Collection<String>> nearest=findClosest(barcodes, distance);
		Collection<String> sets=makeSets(nearest);
		return sets;
	}
	

	
	private static Map<String, Pair<Double>> fractionSpliced(Map<String, Collection<String>> collapsedUMIs, Map<String, Map<String, String>> sets) {
		Map<String, Pair<Double>> rtrn=new TreeMap<String, Pair<Double>>();
		
		
		for(String barcode: collapsedUMIs.keySet()) {
			Collection<String> umis=collapsedUMIs.get(barcode);
			double spliced=0;
			double total=0;
			
			//for each UMI get spliced/unspliced tag
			
			for(String umi: umis) {
				Collection<String> type=new TreeSet<String>();
				String[] tokens=umi.split(" ");
				for(int i=0; i<tokens.length; i++) {
					type.add(sets.get(barcode).get(tokens[i]));
				}
				if(type.size()==1) {
					String t=type.iterator().next();
					if(t.equals("spliced")) {spliced++;}
					total++;
				}
			}
			
			Pair<Double> pair=new Pair<Double>(spliced, total);
			rtrn.put(barcode, pair);
		}
		
		return rtrn;
	}
	
	private static void write(Map<String, Pair<Double>>[] fractionSpliced, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: fractionSpliced[0].keySet()) {
			writer.write(barcode);
			for(int i=0; i<fractionSpliced.length; i++) {
				Pair<Double> pair=fractionSpliced[i].get(barcode);
				double ratio=pair.getValue1()/pair.getValue2();
				writer.write("\t"+pair.getValue2()+"\t"+ratio);
			}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			List<String> lines=BEDFileIO.loadLines(args[0]);
			Map<String, Map<String, String>> sets=collapseBarcodes(lines, Integer.parseInt(args[2]));
			
			int[] umiDistances= {0,1,2};
			Map<String, Pair<Double>>[] fractionSpliced=new Map[umiDistances.length];
			
			for(int i=0; i<umiDistances.length; i++) {
				Map<String, Collection<String>> umis=collapseUMIs(sets, umiDistances[i]);
				fractionSpliced[i]=fractionSpliced(umis, sets);
			}
			write(fractionSpliced, args[1]);
		}
		else {System.err.println(usage);}
	}

	

	static String usage=" args[0]=cDNA file \n args[1]=save \n args[2]=distance";

	
	
}

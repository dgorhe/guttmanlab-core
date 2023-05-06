package guttmanlab.core.proteinSPRITE;

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
import guttmanlab.core.rnasprite.Cluster;

public class Random {

	public Random() throws IOException {
		
		Collection<String[]> randomBarcodes=random(3358417);
		match1Missing(randomBarcodes);
		
	}
	
	private Collection<String[]> random(int size) {
		Collection<String[]> rtrn=new ArrayList<String[]>();
		for(int i=0; i<size; i++) {
			String[] barcode=random();
			rtrn.add(barcode);
			if(i%100000==0) {System.err.println(i);}
		}
		return rtrn;
	}

	private String[] random() {
		String[] rtrn=new String[6];
		
		for(int i=0; i<rtrn.length; i++) {
			String round="Round"+i;
			rtrn[i]=round+"_"+new Double(Math.random()*24).intValue();
		}
		
		return rtrn;
	}

	private void write(String save, Map<String, Collection<String>> updated) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String sub: updated.keySet()) {
			Collection<String> list=updated.get(sub);
			writer.write(sub+"\t"+list.size());
			for(String full: list) {writer.write("\t"+full);}
			writer.write("\n");
		}
		
		writer.close();
	}

	private void match1Missing(Collection<String[]> barcodes) {
		int total=6;
		
		//Map<String, Collection<String>> truncated=new TreeMap<String, Collection<String>>();
		
		
		
		for(int skip=0; skip<total; skip++) {
			System.err.println(skip);
			Collection<String> set=new TreeSet<String>();
			int counter=0;
			for(String[] barcode: barcodes) {
				//String full=makeString(barcode);
				String sub=makeString(barcode, skip);
				
				/*if(!truncated.containsKey(sub)) {
					truncated.put(sub, new TreeSet<String>());
				}
				Collection<String> list=truncated.get(sub);
				list.add(full);
				truncated.put(sub, list);*/
				set.add(sub);
				counter++;
				if(counter%100000==0) {System.err.println(counter+" "+barcodes.size());}
			}
			System.out.println(skip+" "+barcodes.size()+" "+set.size());
		}

		//System.err.println(barcodes.size()+" "+set.size());
		
		//return truncated;
		
	}

	private void print(Map<String, Integer> dist) {
		for(String b: dist.keySet()) {
			System.out.println(b+"\t"+dist.get(b));
		}
		
	}

	private String makeString(String[] b) {
		String rtrn="";
		
		for(int i=0; i<b.length; i++) {
			if(i>0) {rtrn+=".";}
			rtrn+=b[i];
		}
		return rtrn;
	}

	private String makeString(String[] b, int skip) {
		String rtrn="";
		
		int start=0;
		for(int i=0; i<b.length; i++) {
			if(i!=skip) {
				if(start>0) {rtrn+=".";}
				rtrn+=b[i];
				start++;
			}
		}
		return rtrn;
	}
	
	
	private Map<String, Integer> minDistance(Collection<String[]> barcodes) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		Map<String, Collection<String>> indexed0=indexBy(barcodes, 0);
		Map<String, Collection<String>> indexed1=indexBy(barcodes, 1);
		Map<String, Collection<String>> indexed2=indexBy(barcodes, 2);
		Map<String, Collection<String>> indexed3=indexBy(barcodes, 3);
		Map<String, Collection<String>> indexed4=indexBy(barcodes, 4);
		Map<String, Collection<String>> indexed5=indexBy(barcodes, 5);
		
		
		for(String[] b1: barcodes) {
			Collection<String> list0=indexed0.get(b1[0]);
			Collection<String> list1=indexed1.get(b1[1]);
			Collection<String> list2=indexed2.get(b1[2]);
			Collection<String> list3=indexed3.get(b1[3]);
			Collection<String> list4=indexed4.get(b1[4]);
			Collection<String> list5=indexed5.get(b1[5]);
			
			String barcode=makeString(b1);
			int minDistance=score(list0, list1, list2, list3, list4, list5, barcode);
			System.out.println(barcode+"\t"+minDistance);
			/*int minDist=Integer.MAX_VALUE;
			
			for(String[] b2: barcodes) {
				minDist=Math.min(minDist, dist(b1, b2));
			}*/
			rtrn.put(barcode, minDistance);
		}
		
		return rtrn;
	}

	private int score(Collection<String> list0, Collection<String> list1, Collection<String> list2,
			Collection<String> list3, Collection<String> list4, Collection<String> list5, String barcode) {
		
		
		Collection<String> all=new TreeSet<String>();
		all.addAll(list0);
		all.addAll(list1);
		all.addAll(list2);
		all.addAll(list3);
		all.addAll(list4);
		all.addAll(list5);
		
		int minDist=Integer.MAX_VALUE;
		Map<String, Integer> scores=new TreeMap<String, Integer>();
		for(String b: all) {
			int score=0;
			if(list0.contains(b)) {score++;}
			if(list1.contains(b)) {score++;}
			if(list2.contains(b)) {score++;}
			if(list3.contains(b)) {score++;}
			if(list4.contains(b)) {score++;}
			if(list5.contains(b)) {score++;}
			if(!b.equals(barcode)) {
				//System.err.println(barcode+"\t"+b+"\t"+score);
				scores.put(b, 6-score);
				minDist=Math.min(minDist, 6-score);
			}
		}
		return minDist;
		
	}

	private Map<String, Collection<String>> indexBy(Collection<String[]> barcodes, int i) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String[] b: barcodes) {
			Collection<String> list=new ArrayList<String>();
			if(rtrn.containsKey(b[i])){
				list=rtrn.get(b[i]);
			}
			list.add(makeString(b));
			rtrn.put(b[i], list);
		}
		
		return rtrn;
	}

	private int dist(String[] b1, String[] b2) {
		int dist=0;
		for(int i=0; i<b1.length; i++) {
			if(b1[i]!=b2[i]) {dist++;}
		}
		return dist;
	}

	private void print(Map<String, Double> fraction1, Map<String, String> names) {
		for(String name: fraction1.keySet()) {
			System.err.println(name+"\t"+names.get(name)+"\t"+fraction1.get(name));
		}
		
	}

	private Map<String, Double> getFraction(Collection<String[]> barcodes, int i) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		for(String[] b: barcodes) {
			String name=b[i];
			double count=0;
			if(rtrn.containsKey(name)) {count=rtrn.get(name);}
			count++;
			rtrn.put(name, count);
		}
		return rtrn;
	}

	private boolean isCorrectOrder(String[] barcode, Map<String, String> barcodeToRound) {
		for(int i=1; i<barcode.length; i++) {
			if(barcodeToRound.containsKey(barcode[i])) {
				String name=barcodeToRound.get(barcode[i]);
				if(i==1 && !name.contains("ROUND5")) {return false;}
				if(i==2 && !name.contains("ROUND4")) {return false;}
				if(i==3 && !name.contains("ROUND3")) {return false;}
				if(i==4 && !name.contains("ROUND2")) {return false;}
				if(i==5 && !name.contains("ROUND1")) {return false;}
			}
		}
		return true;
	}

	private String getRounds(String[] barcode, Map<String, String> barcodeToRound) {
		String rtrn="";
		
		for(String b: barcode) {
			if(barcodeToRound.containsKey(b)) {
				if(!rtrn.isEmpty()) {rtrn+=".";}
				rtrn+=barcodeToRound.get(b);
			}
		}
		return rtrn;
	}

	
	private static Map<String, String> parse(String string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		List<String> lines=BEDFileIO.loadLines(string);
		
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			rtrn.put(tokens[0], tokens[1]);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		new Random();
	}

	
}

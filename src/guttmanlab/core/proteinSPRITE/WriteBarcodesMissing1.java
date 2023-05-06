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

public class WriteBarcodesMissing1 {

	public WriteBarcodesMissing1(BarcodingDataStreaming data, String save) throws IOException {
		FileWriter writer1=new FileWriter(save+".M1.clusters");
		FileWriter writer2=new FileWriter(save+".M2.clusters");
		FileWriter writer3=new FileWriter(save+".M3.clusters");
		FileWriter writer4=new FileWriter(save+".M4.clusters");
		FileWriter writer5=new FileWriter(save+".M5.clusters");
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			String name1=sub(c,1);
			String name2=sub(c,2);
			String name3=sub(c,3);
			String name4=sub(c,4);
			String name5=sub(c,5);
			
			writer1.write(c.toString(name1)+"\n");
			writer2.write(c.toString(name2)+"\n");
			writer3.write(c.toString(name3)+"\n");
			writer4.write(c.toString(name4)+"\n");
			writer5.write(c.toString(name5)+"\n");
			
			/*String[] barcode=c.getBarcode().split("\\.");
			String rounds=getRounds(barcode, barcodeToRound);
			boolean isCorrectOrder=isCorrectOrder(barcode, barcodeToRound);
			if(isCorrectOrder) {
				//barcodes.add(barcode);
				if(c.getClusterSize()>4) {writer4.write(c.getBarcode()+"\t"+rounds+"\t"+c.getClusterSize()+"\t"+c.getAllDNAIntervals().size()+"\t"+c.getProteins().size()+"\n");}
				else{writer3.write(c.getBarcode()+"\t"+rounds+"\t"+c.getClusterSize()+"\t"+c.getAllDNAIntervals().size()+"\t"+c.getProteins().size()+"\n");}
			}
			else if(rounds.contains("ROUND1") && rounds.contains("ROUND2") && rounds.contains("ROUND3")&& rounds.contains("ROUND4")&& rounds.contains("ROUND5")) {
				writer2.write(c.getBarcode()+"\t"+rounds+"\t"+c.getClusterSize()+"\t"+c.getAllDNAIntervals().size()+"\t"+c.getProteins().size()+"\n");
			}
			else {
				String[] corrected=correct(barcode, barcodeToRound);
				missingBarcodes.put(c.getBarcode(), corrected);
				writer1.write(c.getBarcode()+"\t"+rounds+"\t"+makeString(corrected)+"\t"+c.getClusterSize()+"\t"+c.getAllDNAIntervals().size()+"\t"+c.getProteins().size()+"\n");
			}*/
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		writer1.close();
		writer2.close();
		writer3.close();
		writer4.close();
		writer5.close();
		
	}
	
	private String sub(Cluster c, int skip) {
		String[] barcode=c.getBarcode().split("\\.");
		return makeString(barcode, skip);
	}

	private Map<String, Collection<String>> match(Map<String, String[]> missingBarcodes, BarcodingDataStreaming data) {
		Map<String, Collection<String>> subToFull=new TreeMap<String, Collection<String>>();
		
		Collection<String> errorSet=new TreeSet<String>();
		
		int count=0;
		for(String barcode: missingBarcodes.keySet()) {
			String[] missing=missingBarcodes.get(barcode);
			String sub=find(missing);
			errorSet.add(sub);
			count++;
			if(count%10000==0) {System.err.println(count+" "+missingBarcodes.size());}
		}
		
		count=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Collection<String> subMatches=match(c, errorSet);
			for(String sub: subMatches) {
				System.err.println(sub);
				add(sub, c.getBarcode(), subToFull);
			}
			count++;
			if(count%10000==0) {System.err.println(count);}
		}
		
		return subToFull;
	}
	
	private Collection<String> match(Cluster c, Collection<String> errorSet) {
		Collection<String> rtrn=new TreeSet<String>();
		Collection<String> subs=getSub(c.getBarcode().split("\\."));
		
		for(String sub: subs) {
			if(errorSet.contains(sub)) {rtrn.add(sub);}
		}
		
		return rtrn;
	}

	private Collection<String> getSub(String[] split) {
		Collection<String> rtrn=new TreeSet<String>();
		for(int i=0; i<split.length; i++) {
			rtrn.add(makeString(split, i));
		}
		return rtrn;
	}

	private void add(String sub, String barcode, Map<String, Collection<String>> subToFull) {
		if(!subToFull.containsKey(sub)) {subToFull.put(sub, new TreeSet<String>());}
		Collection<String> list=subToFull.get(sub);
		list.add(barcode);
	}

	private String find(String[] missing) {
		int count=0;
		String rtrn="";
		for(int i=0; i<missing.length; i++) {
			if(missing[i].contains("NF")) {
				count++;
			}
			else {
				if(count!=0) {rtrn+=".";}
				rtrn+=missing[i];
				count++;
			}
		}
		
		return rtrn;
		
	}

	private void find(BarcodingDataStreaming data, String[] missing, String originalBarcode) {
		int[] skip=new int[missing.length];
		int count=0;
		for(int i=0; i<missing.length; i++) {
			if(missing[i].contains("NF")) {
				skip[i]=1;
				count++;
			}
		}
		
		if(count==1) {
			
			while(data.hasNext()) {
				Cluster c=data.next();
				String[] barcode=c.getBarcode().split("\\.");
				boolean match=match(missing, barcode, skip);
				if(match) {
					System.out.println(originalBarcode+" "+c.getBarcode() +" "+c.getClusterSize()+" "+c.getAllDNAIntervals().size()+" "+c.getProteins().size());
				}
			}
			data.close();
		}
		
	}

	private boolean match(String[] missing, String[] barcode, int[] skip) {
		for(int i=0; i<skip.length; i++) {
			if(skip[i]!=1) {
				if(!missing[i].equals(barcode[i])) {return false;}
			}
		}
		return true;
	}

	private String[] correct(String[] barcode, Map<String, String> barcodeToRound) {
		int countNF=0;
		String[] rtrn=new String[barcode.length-1];
		for(int i=0; i<barcode.length-1; i++) {
			rtrn[i]=barcode[i];
			if(barcodeToRound.containsKey(barcode[i])) {
				String name=barcodeToRound.get(barcode[i]);
				if(i==1 && !name.contains("ROUND5")) {rtrn[i]="ROUND5_NF"; countNF++;}
				if(i==2 && !name.contains("ROUND4")) {rtrn[i]="ROUND4_NF"; countNF++;}
				if(i==3 && !name.contains("ROUND3")) {rtrn[i]="ROUND3_NF"; countNF++;}
				if(i==4 && !name.contains("ROUND2")) {rtrn[i]="ROUND2_NF"; countNF++;}
				if(i==5 && !name.contains("ROUND1")) {rtrn[i]="ROUND1_NF"; countNF++;}
			}
		}
		
		//if(countNF==1) {System.out.println(makeString(rtrn));}
		//System.out.println(countNF);
		
		return rtrn;
	}

	private Collection<String[]> random(int size) {
		Collection<String[]> rtrn=new ArrayList<String[]>();
		for(int i=0; i<size; i++) {
			String[] barcode=random();
			rtrn.add(barcode);
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
			System.err.println(skip+" "+barcodes.size()+" "+set.size());
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
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		new WriteBarcodesMissing1(data, save);
	}

	
}

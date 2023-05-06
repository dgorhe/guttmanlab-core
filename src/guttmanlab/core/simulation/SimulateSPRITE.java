package guttmanlab.core.simulation;

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

public class SimulateSPRITE {

	private static int numMolecules=1000000;
	private static int clusterSize=5;
	double ligationEfficiency=0.01;
	int numberRounds=8;
	int numSplits=12;
	double errorRate=0.05;
	int mismatches=1;
	
	public SimulateSPRITE(Map<String, ArrayList<String>> clusters, String save) throws IOException{
		//barcode each molecule
		FileWriter writer=new FileWriter(save);
		Map<String, String> moleculeToBarcode=new TreeMap<String, String>();
		
		Map<Integer, Integer> distanceMap=new TreeMap<Integer, Integer>();
		
		Map<String, Collection<String>> errorMap=new TreeMap<String, Collection<String>>();
		
		for(String cluster: clusters.keySet()){
			Collection<String> errors=new TreeSet<String>();
			String barcode=getBarcodeString();
			//for each molecule in cluster add barcode with efficiency
			for(String molecule: clusters.get(cluster)){
				//String barcode2=addDropout(barcode);
				String barcode2=addError(barcode);
				
				int distance=distance(barcode, barcode2);
				int count=0;
				if(distanceMap.containsKey(distance)){count=distanceMap.get(distance);}
				count++;
				distanceMap.put(distance, count);
				
				if(distance<=mismatches){
					errors.add(barcode2);
				}
				
				errorMap.put(barcode, errors);
				writer.write(molecule+"\t"+cluster+"\t"+barcode2+"\n");
				moleculeToBarcode.put(molecule, barcode2);
			}
		}
		
		for(Integer distance: distanceMap.keySet()){
			System.err.println(distance+" "+distanceMap.get(distance));
		}
		
		writer.close();
		
		Map<String, Collection<String>> barcodes=new TreeMap<String, Collection<String>>();
		
		for(String molecule: moleculeToBarcode.keySet()){
			String cluster=moleculeToBarcode.get(molecule);
			Collection<String> list=new ArrayList<String>();
			if(barcodes.containsKey(cluster)){list=barcodes.get(cluster);}
			list.add(molecule);
			barcodes.put(cluster, list);
		}
		
		writeBarcodes(save+".clusters", barcodes, errorMap);
	
		//findCloseBarcodes(barcodes, save+".collapsed");
	}
	
	private void findCloseBarcodes(Map<String, Collection<String>> barcodes, String save) throws IOException {
		//Collapse close
		Map<String, String> collapsed=new TreeMap<String, String>();
		
		int count=0;
		for(String barcode: barcodes.keySet()){
			String collapsedCluster="cluster"+count;
			Collection<String> close=findClose(barcode, barcodes.keySet(), 1);
			count++;
			
			for(String other: close){
				if(!collapsed.containsKey(other)){
					collapsed.put(other, collapsedCluster);
				}
			}
		}
		
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: collapsed.keySet()){
			String collapsedName=collapsed.get(barcode);
			writer.write(barcode+"\t"+collapsedName+"\t"+barcodes.get(barcode)+"\n");
		}
		
		writer.close();
	}

	private int size(Map<String, Collection<String>> barcodes, String barcode, Collection<String> close) {
		Collection<String> list=new TreeSet<String>();
		list.addAll(barcodes.get(barcode));
		for(String other: close){
			list.addAll(barcodes.get(other));
		}
		return list.size();
	}

	private Collection<String> findClose(String barcode, Set<String> keySet, int i) {
		Collection<String> rtrn=new TreeSet<String>(); 
		for(String other: keySet){
			if(distance(barcode, other)<=i){rtrn.add(other);}
		}
		return rtrn;
	}

	private int distance(String barcode, String other) {
		int count=0;
		for(int i=0; i<barcode.toCharArray().length; i++){
			if(barcode.toCharArray()[i]!=other.toCharArray()[i]){count++;}
		}
		return count;
	}

	private void writeBarcodes(String save, Map<String, Collection<String>> barcodes, Map<String, Collection<String>> errorMap) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		
		
		for(String cluster: barcodes.keySet()){
			writer.write(cluster);
			for(String molecule: barcodes.get(cluster)){
				writer.write("\t"+molecule);
			}
			writer.write("\n");
		}
		
		writer.close();
		
		
		
		writer=new FileWriter(save+".corrected");
		
		for(String barcode: errorMap.keySet()){
			Collection<String> allMolecules=new TreeSet<String>();
			for(String corrected: errorMap.get(barcode)){
				Collection<String> molecules=barcodes.get(corrected);
				allMolecules.addAll(molecules);
			}
			writer.write(barcode);
			for(String molecule: allMolecules){
				writer.write("\t"+molecule);
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	private String addError(String barcode) {
		String rtrn="";
		for(int i=0; i<barcode.length(); i++){
			String pos=barcode.substring(i, i+1);
			if(Math.random()<errorRate){
				pos=toChar(new Double(Math.random()*numSplits).intValue());
			}
			rtrn+=pos;
		}
		return rtrn;
	}
	
	private String addDropout(String barcode) {
		String rtrn="";
		for(int i=0; i<barcode.length(); i++){
			String pos=barcode.substring(i, i+1);
			if(Math.random()<ligationEfficiency){
				pos="-";
				rtrn+=pos;
				return rtrn;
			}
			rtrn+=pos;
		}
		return rtrn;
	}

	private String getBarcodeString() {
		String barcode="";
		for(int i=0; i<numberRounds; i++){
			int well=new Double(Math.random()*numSplits).intValue();
			String round=toChar(well); 
			barcode+=round;
		}
		return barcode;
	}

	

	private String toChar(int well) {
		if(well==0){return "A";}
		if(well==1){return "B";}
		if(well==2){return "C";}
		if(well==3){return "D";}
		if(well==4){return "E";}
		if(well==5){return "F";}
		if(well==6){return "G";}
		if(well==7){return "H";}
		if(well==8){return "I";}
		if(well==9){return "J";}
		if(well==10){return "K";}
		if(well==11){return "L";}
		return "M";
	}

	private void write(String save, Map<String, String> moleculeToBarcode) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String molecule: moleculeToBarcode.keySet()){
			String barcode=moleculeToBarcode.get(molecule);
			writer.write(molecule+"\t"+barcode+"\n");
		}
		
		writer.close();
	}
	
	private static Map<String, ArrayList<String>> parseClusters(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		
		Map<String, ArrayList<String>> temp=new TreeMap<String, ArrayList<String>>();
		for(String line: lines){
			String molecule=line.split("\t")[0];
			String cluster=line.split("\t")[1];
			ArrayList<String> list=new ArrayList<String>();
			if(temp.containsKey(cluster)){list=temp.get(cluster);}
			list.add(molecule);
			temp.put(cluster, list);
		}
		
		return temp;
	}
	
	private static Map<String, ArrayList<String>> generateClusters() {
		Map<String, ArrayList<String>> rtrn=new TreeMap<String, ArrayList<String>>();
		for(int i=0; i<numMolecules; i++){
			String molecule="M"+i;
			String cluster="C"+(i/clusterSize);
			ArrayList<String> list=new ArrayList<String>();
			if(rtrn.containsKey(cluster)){
				list=rtrn.get(cluster);
			}
			list.add(molecule);
			rtrn.put(cluster, list);
		}
		return rtrn;
	}
	

	public static void main(String[] args) throws IOException{
		
		
		Map<String, ArrayList<String>> clusters=generateClusters();
		
		System.err.println("Generated clusters "+clusters.size());
		
		
		String save="/Users/mguttman/Desktop/test";
		new SimulateSPRITE(clusters, save);
	}

	
}

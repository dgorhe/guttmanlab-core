package guttmanlab.core.util;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;

public class BarcodeGenerator {

	public BarcodeGenerator(int length, int number, int distance, List<String> listToStartWith, String save) throws IOException{
		Collection<String> allBarcodes=generateBarcodes(length, 1000*number);
		
		List<String> barcodes=exclude(allBarcodes);
		
		barcodes=filterHomopolymer(barcodes);
		
		Collection<String> list=pick(barcodes, number, distance, listToStartWith);
		
		FileWriter writer=new FileWriter(save);
		for(String b: list){writer.write(b+"\n");}
		writer.close();
	}

	private List<String> filterHomopolymer(List<String> barcodes) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String barcode: barcodes){
			if(!homopolymer(barcode)){rtrn.add(barcode);}
		}
		
		return rtrn;
	}

	private boolean homopolymer(String barcode) {
		if(barcode.contains("AAA")|| barcode.contains("TTT") || barcode.contains("CCC")|| barcode.contains("GGG")){return true;}
		return false;
	}

	private Collection<String> pick(List<String> barcodes, int number, int distance, List<String> listToStartWith) {
		Collection<String> rtrn=new TreeSet<String>();
		
		//rtrn.addAll(listToStartWith);
		barcodes=exclude(barcodes, listToStartWith, distance);
		
		while(rtrn.size()<number && barcodes.size()>0){
			//randomly pick
			String barcode=randomPick(barcodes);
			rtrn.add(barcode);
			//exclude all within distance
			barcodes=exclude(barcodes, barcode, distance);
			System.err.println(rtrn.size()+" "+barcode+" "+barcodes.size());
		}
		
		return rtrn;
	}

	private String randomPick(List<String> barcodes) {
		int pos=new Double(Math.random()*barcodes.size()).intValue();
		return barcodes.get(pos);
	}

	
	private List<String> exclude(List<String> barcodes, List<String> listToStartWith, int distance) {
		for(String barcode: listToStartWith){
			barcodes=exclude(barcodes, barcode, distance);
		}
		
		return barcodes;
	}
	
	private List<String> exclude(List<String> barcodes, String barcode, int distance) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String barcode1: barcodes){
			int d=distance(barcode1, barcode);
			if(d>distance){rtrn.add(barcode1);}
		}
		
		return rtrn;
	}

	private int distance(String barcode1, String barcode) {
		int distance=0;
		//System.err.println(barcode1+" "+barcode);
		for(int i=0; i<barcode.length(); i++){
			if(barcode.toCharArray()[i]!=barcode1.toCharArray()[i]){distance++;}
		}
		return distance;
	}

	
	private List<String> exclude(Collection<String> barcodes) {
		List<String> rtrn=new ArrayList<String>();
		
		
		for(String barcode: barcodes){
			rtrn.add(barcode);
		}
		
		System.err.println(barcodes.size());
		
		return rtrn;
	}
	
	private List<String> exclude(Collection<String> barcodes, Collection<String> exclude) {
		List<String> rtrn=new ArrayList<String>();
		
		
		for(String barcode: barcodes){
			
			if(!overlaps(exclude, barcode)){rtrn.add(barcode);}
			
			//if(!exclude.contains(barcode)){rtrn.add(barcode);}
		}
		
		System.err.println(barcodes.size()+" "+exclude.size() +" "+rtrn.size());
		
		return rtrn;
	}

	private boolean overlaps(Collection<String> exclude, String barcode) {
		for(String barcode1: exclude){
			if(overlap(barcode1, barcode)){
				//System.err.println(barcode+" b1 "+barcode1);
				return true;
			}
		}
		return false;
	}

	private boolean overlap(String barcode1, String barcode) {
		Pair<String> largerSmaller=getLargerSmaller(barcode, barcode1);
		//String larger=getLarger(barcode, barcode1);
		
		int index=largerSmaller.getValue1().indexOf(largerSmaller.getValue2());
		//System.err.println(barcode1 +" "+barcode+" "+index +" "+largerSmaller.getValue1()+" "+largerSmaller.getValue2());
		//System.err.println(i);
		if(index>=0){return true;}
		return false;
	}

	private Pair<String> getLargerSmaller(String barcode, String barcode1) {
		Pair<String> rtrn=new Pair<String>();
		
		rtrn.setValue1(barcode1);
		rtrn.setValue2(barcode);
		
		if(barcode.length()>barcode1.length()){
			rtrn.setValue1(barcode);
			rtrn.setValue2(barcode1);
		}
		return rtrn;
	}
	
	private String getLarger(String barcode, String barcode1) {
		if(barcode.length()>barcode1.length()){return barcode;}
		return barcode1;
	}

	private Collection<String> generateBarcodes(int length, int num) {
		Collection<String> rtrn=new TreeSet<String>();
		for(int i=0; i<num; i++){
			String barcode=random(length);
			rtrn.add(barcode);
		}
		return rtrn;
	}

	private String random(int length) {
		String rtrn="";
		for(int i=0; i<length; i++){
			rtrn+=random();
		}
		return rtrn;
	}

	private String random() {
		double ran=Math.random();
		if(ran<0.25){return "A";}
		if(ran>.25 && ran<0.5){return "C";}
		if(ran>.5 && ran<0.75){return "G";}
		return "T";
	}
	
	private static List<String> parse(String file) throws IOException {
		List<String> rtrn=new ArrayList<String>();
		List<String> lines=BEDFileIO.loadLines(file,1);
		
		for(String line: lines){
			//System.err.println("line: "+line);
			rtrn.add(line.split("\t")[3]);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		int length=16;
		int number=336;
		int distance=5;
		
		List<String> listToStartWith=parse(args[0]);
		String save=args[1];
		
		new BarcodeGenerator(length, number, distance, listToStartWith, save);
	}

	


}

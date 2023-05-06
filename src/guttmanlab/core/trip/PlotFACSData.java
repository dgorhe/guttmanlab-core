package guttmanlab.core.trip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;

public class PlotFACSData {
	
	public PlotFACSData(File[] files, double intervalSize, String save) throws IOException {
		TreeSet<Integer> bins=new TreeSet<Integer>();
		Map<Integer, List<Double>>[] smoothed=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			List<Pair<Double>> list=parse(files[i]);
			smoothed[i]=smooth(list, intervalSize);
			bins.addAll(smoothed[i].keySet());
		}
		write(save, bins, smoothed, files);
	}

	private void write(String save, TreeSet<Integer> bins, Map<Integer, List<Double>>[] smoothed, File[] files) throws IOException {
		FileWriter writer=new FileWriter(save);

		writer.write("bin");
		for(int i=0; i<files.length; i++) {
			writer.write("\t"+files[i].getName());
		}
		writer.write("\n");
		
		for(Integer bin: bins) {
			writer.write(bin.toString());
			for(int i=0; i<smoothed.length; i++) {
				double val=get(smoothed[i], bin);
				if(val>0) {writer.write("\t"+val);}
				else {writer.write("\t ");}
			}
			writer.write("\n");
		}
		
		writer.close();
		
	}

	private double get(Map<Integer, List<Double>> map, Integer bin) {
		if(map.containsKey(bin)) {
			List<Double> vals=map.get(bin);
			double avg=Statistics.quantile(vals, 0.5);
			return avg;
		}
		return 0;
	}

	

	private Map<Integer, List<Double>> smooth(List<Pair<Double>> list, double intervalSize){
		Map<Integer, List<Double>> rtrn=new TreeMap<Integer, List<Double>>();
		
		for(Pair<Double> pair: list){
			double ratio=pair.getValue1()/intervalSize;
			int bin=Double.valueOf(ratio).intValue();
			//System.err.println(pair.getValue1()+" "+ratio+" "+bin);
			if(!rtrn.containsKey(bin)) {rtrn.put(bin, new ArrayList<Double>());}
			
			List<Double> temp=rtrn.get(bin);
			temp.add(pair.getValue2());
			rtrn.put(bin, temp);
		}
		
		return rtrn;
	}
	
	

	private List<Pair<Double>> parse(File input) throws IOException {
		List<Pair<Double>> rtrn=new ArrayList<Pair<Double>>();
		List<String> lines=BEDFileIO.loadLines(input.getAbsolutePath(), 1);
		
		for(String line: lines) {
			String[] tokens=line.split(",");
			Pair<Double> pair=new Pair<Double>();
			pair.setValue1(Double.parseDouble(tokens[1]));
			pair.setValue2(Double.parseDouble(tokens[0]));
			rtrn.add(pair);
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		File[] input=new File(args[0]).listFiles();
		double intervalSize=Double.valueOf(args[1]);
		String save=args[2];
		new PlotFACSData(input, intervalSize, save);
	}
	
}

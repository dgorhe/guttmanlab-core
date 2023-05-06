package guttmanlab.core.trip;

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
import guttmanlab.core.math.Statistics;

public class MS2FACSData {

	public static void smooth(File[] files, String save) throws IOException {
		
		Map<Integer, Collection<Double>>[] values=new Map[files.length];
		for(int i=0; i< files.length; i++) {
			List<String> lines=BEDFileIO.loadLines(files[i].getAbsolutePath(),1);
			values[i]=parse(lines);
		}
		
		write(save, values, files);
		
	}

	private static void write(String save, Map<Integer, Collection<Double>>[] values, File[] files) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		
		Collection<Integer> positions=getPositions(values);
		
		writer.write("Index");
		
		for(int i=0; i<files.length; i++) {
			writer.write("\t"+files[i].getName());
		}
		writer.write("\n");
		
		for(Integer pos: positions) {
			writer.write(pos.toString());
			for(int i=0; i<values.length; i++) {
				Collection<Double> vals=get(values[i],pos);
				if(vals.isEmpty()) {writer.write("\t");}
				else {
					double avg=Statistics.mean(vals);
					writer.write("\t"+avg);
				}
			}
			writer.write("\n");
		}
		
		writer.close();
		
	}

	private static Collection<Integer> getPositions(Map<Integer, Collection<Double>>[] values) {
		Collection<Integer> rtrn=new TreeSet<Integer>();
		
		for(int i=0; i<values.length; i++) {rtrn.addAll(values[i].keySet());}
		
		return rtrn;
	}

	private static Collection<Double> get(Map<Integer, Collection<Double>> map, Integer pos) {
		if(map.containsKey(pos)) {return map.get(pos);}
		return new ArrayList<Double>();
	}

	private static Map<Integer, Collection<Double>> parse(List<String> lines) {
		Map<Integer, Collection<Double>> map=new TreeMap<Integer, Collection<Double>>();
		
		for(String line: lines) {
			String[] tokens=line.split(",");
			Double y=Double.parseDouble(tokens[2]);
			Double x=Double.parseDouble(tokens[1])*10;
			int binned=x.intValue();
			//System.out.println(x+" "+binned);
			
			
			
			if(!map.containsKey(binned)) {map.put(binned, new ArrayList<Double>());}
			
			Collection<Double> list=map.get(binned);
			list.add(y);
		}
		
		return map;
	}

	private static Collection<Integer> getPositions(Map<String, Map<Integer, Collection<Double>>> valuesBySample) {
		Collection<Integer> rtrn=new TreeSet<Integer>();
		
		for(String sample: valuesBySample.keySet()) {
			rtrn.addAll(valuesBySample.get(sample).keySet());
		}
		
		return rtrn;
	}

	private static Collection<Double> get(Map<String, Map<Integer, Collection<Double>>> valuesBySample, String sample, Integer pos) {
		if(valuesBySample.get(sample).containsKey(pos)) {return valuesBySample.get(sample).get(pos);}
		return new ArrayList<Double>();
	}
	
	
	public static void main(String[] args) throws IOException {
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		smooth(files, save);
	}
	
}

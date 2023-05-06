package guttmanlab.core.spidr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class EnrichmentPerWindow {
	double c=100000;

	public EnrichmentPerWindow(File sample, File[] controls) throws NumberFormatException, IOException {
		Map<String, Double> sampleScores=getScores(sample);
		Map<String, Double> controlScores=getScores(controls);
		
		writeEnrichment(sampleScores, controlScores);
		
	}
	
	public EnrichmentPerWindow(File[] files, File[] controls, String save) throws NumberFormatException, IOException {
		Map<String, Double>[] enrichments=new Map[files.length];
		Map<String, Double> controlScores=getScores(controls);
		for(int i=0; i<files.length; i++) {
			System.err.println(i+" "+files.length+" "+files[i].getName());
			Map<String, Double> sampleScores=getScores(files[i]);
			enrichments[i]=getEnrichments(sampleScores, controlScores);
		}
		
		MatrixWithHeaders matrix=makeMatrix(enrichments, files);
		matrix.write(save);
	}
	
	
	public EnrichmentPerWindow(File[] files, File[] controls, List<String> regions, String save) throws NumberFormatException, IOException {
		Map<String, Double> controlScores=getScores(controls);
		double controlMin=Statistics.min(controlScores.values());
		
		MatrixWithHeaders matrix=new MatrixWithHeaders(regions, getColumns(files));
		
		for(int i=0; i<files.length; i++) {
			System.err.println(i+" "+files.length+" "+files[i].getName());
			Map<String, Double> sampleScores=getScores(files[i]);
			double sampleMin=Statistics.min(sampleScores.values());
			double min=Math.max(sampleMin, controlMin);
			
			for(String region: regions) {
				double enrich=getEnrichment(sampleScores, controlScores, region, min);
				matrix.set(region, files[i].getName(), enrich);
			}
		}
		
		matrix.write(save);
	}

	private double getEnrichment(Map<String, Double> sampleScores, Map<String, Double> controlScores, String region, double min) {
		double o=Math.max(min, get(sampleScores,region));
		double e=Math.max(min, get(controlScores,region));
		double enrich=o/e;
		double diff=Math.max(0, o-e);
		double log=Math.max(0, Math.log(enrich)/Math.log(2));
		return diff;
	}

	private List<String> getColumns(File[] files) {
		List<String> columns=new ArrayList<String>();
		for(int i=0; i<files.length; i++) {columns.add(files[i].getName());}
		return columns;
	}

	private MatrixWithHeaders makeMatrix(Map<String, Double>[] enrichments, File[] files) {
		List<String> columns=new ArrayList<String>();
		for(int i=0; i<files.length; i++) {columns.add(files[i].getName());}
		
		Collection<String> set=new TreeSet<String>();
		for(int i=0; i<enrichments.length; i++) {set.addAll(enrichments[i].keySet());}
		List<String> rows=new ArrayList<String>();
		rows.addAll(set);
		
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(int i=0; i<enrichments.length; i++) {
			String col=files[i].getName();
			for(String row: enrichments[i].keySet()) {
				double score=enrichments[i].get(row);
				rtrn.set(row, col, score);
			}
		}
		
		return rtrn;
	}

	private File[] getFiles(File[] files, int pos) {
		File[] rtrn=new File[files.length-1];
		
		int counter=0;
		for(int i=0; i<files.length; i++) {
			if(i!=pos) {
				rtrn[counter]=files[i];
				counter++;
			}
		}
		
		return rtrn;
	}

	private Map<String, Double> getScores(File[] controls) throws NumberFormatException, IOException {
		Map<String, List<Double>> map=new TreeMap<String, List<Double>>();
		
		for(int i=0; i<controls.length; i++) {
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(controls[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				double score=Double.parseDouble(tokens[1]);
				if(!map.containsKey(tokens[0])) {map.put(tokens[0], new ArrayList<Double>());}
				List<Double> list=map.get(tokens[0]);
				list.add(score);
			}
			reader.close();
		}
		
		
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		double total=sum(map.get("total"));
		for(String region: map.keySet()) {
			double score=sum(map.get(region));
			double norm=c*(score/total);
			rtrn.put(region, norm);
		}
		
		return rtrn;
	}

	private double sum(List<Double> list) {
		return Statistics.sum(list);
	}

	private Map<String, Double> getScores(File sample) throws NumberFormatException, IOException {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		double total=0;
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(sample)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			if(tokens[0].equals("total")) {total=Double.parseDouble(tokens[1]);}
			else {
				double score=Double.parseDouble(tokens[1]);
				double norm=c*(score/total);
				rtrn.put(tokens[0], norm);
			}
		}
		reader.close();
		return rtrn;
	}
	
	private Map<String, Double> getEnrichments(Map<String, Double> sampleScores, Map<String, Double> controlScores) {
		double sampleMin=Statistics.min(sampleScores.values());
		double controlMin=Statistics.min(controlScores.values());
		double min=Math.max(sampleMin, controlMin);
		
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		for(String region: sampleScores.keySet()) {
			double o=Math.max(min, get(sampleScores,region));
			double e=Math.max(min, get(controlScores,region));
			
			double enrich=o/e;
			double diff=Math.max(0, o-e);
			double log=Math.max(0, Math.log(enrich)/Math.log(2));
			rtrn.put(region, diff);
		}
		return rtrn;
	}

	private void writeEnrichment(Map<String, Double> sampleScores, Map<String, Double> controlScores) {
		double sampleMin=Statistics.min(sampleScores.values());
		double controlMin=Statistics.min(controlScores.values());
		double min=Math.max(sampleMin, controlMin);
		for(String region: sampleScores.keySet()) {
			double o=Math.max(min, get(sampleScores,region));
			double e=Math.max(min, get(controlScores,region));
			
			double enrich=o/e;
			double diff=Math.max(0, o-e);
			double log=Math.max(0, Math.log(enrich)/Math.log(2));
			System.out.println(region+"\t"+o+"\t"+e+"\t"+enrich+"\t"+log+"\t"+diff);
		}
		
	}
	
	private double get(Map<String, Double> sampleScores, String region) {
		if(sampleScores.containsKey(region)) {return sampleScores.get(region);}
		return 0;
	}

	public static void main(String[] args) throws NumberFormatException, IOException {
		File[] samples=new File(args[0]).listFiles();
		File[] controls=new File(args[1]).listFiles();
		List<String> regions=BEDFileIO.loadLines(args[2]);
		String save=args[3];
		new EnrichmentPerWindow(samples, controls, regions, save);
	}
	
}

package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.math.Statistics;

public class ReadsByRNA {
	
	private static List<String> parseList(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		
		List<String> rtrn=new ArrayList<String>();
		for(String line: lines) {
			rtrn.add(line.split("\t")[0]);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		int binResolution=Integer.parseInt(args[1]);
		String save=args[2];
		
		int numberOfRegions=Integer.parseInt(args[3]);
		
		
		List<SingleInterval> allRegions=new ArrayList<SingleInterval>();
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			//if(c.containsRNA(RNA)) {
			Collection<SingleInterval> regions=c.getAllDNAIntervals();
			allRegions.addAll(regions);
				//for(SingleInterval r: regions) {System.out.println(r.toShortBED());}
			//}
			counter++;
			if(counter%10000==0) {System.err.println(counter);}
		}
		data.close();
	
		int numPerm=100;
		
		double[] dist=new double[numPerm];
		for(int i=0; i<numPerm; i++) {
			Collection<SingleInterval> random=sample(allRegions, numberOfRegions);
			Map<SingleInterval, Integer> map=count(random, binResolution);
			dist[i]=max(map);
		}
		
		
		write(save, dist);
		
		//BEDFileIO.writeBEDGraphInteger(map, save);
	
	}

	private static double max(Map<SingleInterval, Integer> map) {
		return Statistics.max(map.values());
	}

	private static void write(String save, double[] dist) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(int i=0; i<dist.length; i++) {
			writer.write(dist[i]+"\n");
		}
		
		writer.close();
	}

	private static Collection<SingleInterval> sample(List<SingleInterval> allRegions, int numberOfRegions) {
		Collection<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		for(int i=0; i<numberOfRegions; i++) {
			int index=new Double(Math.random()*(allRegions.size()-1)).intValue();
			rtrn.add(allRegions.get(index));
		}
		return rtrn;
	}

	private static Map<SingleInterval, Integer> count(Collection<SingleInterval> allRegions, int binResolution) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(SingleInterval r: allRegions) {
			SingleInterval binned=r.bin(binResolution);
			int count=0;
			if(rtrn.containsKey(binned)) {count=rtrn.get(binned);}
			count++;
			rtrn.put(binned, count);
		}
		
		return rtrn;
	}

	private static void getClusters(BarcodingDataStreaming data, List<String> rnas, String save) throws IOException {
		//iterate over all clusters
		FileWriter writer=new FileWriter(save);
		int counter=0;
		while(data.hasNext()) {
			boolean use=false;
			Cluster c=data.next();
			Collection<String> rnaNames=c.getRNANames();
			for(String rna: rnaNames) {
				if(rnas.contains(rna)) {
					use=true;
				}
			}
			if(use) {writer.write(c.line+"\n");}	
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		data.close();
		writer.close();
	}
	
}

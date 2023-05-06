package guttmanlab.core.trip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.rnasprite.Kmer;

public class BinPositions {

	public BinPositions(File input, String save, int binSize, BarcodingDataStreaming data, Kmer nuclearBody) throws IOException {
		nuclearBody=bin(nuclearBody, binSize);
		Map<SingleInterval, Pair<Integer>> vals=parse(input);
		vals=binAndMerge(vals, binSize);
		
		
		Map<SingleInterval, Double> speckleCounts=getCounts(data, nuclearBody, binSize, vals.keySet());
		
		write(save, vals, speckleCounts);
	}
	
	private Kmer bin(Kmer nuclearBody, int binResolution2) {
		Kmer rtrn=new Kmer();
		
		for(SingleInterval region: nuclearBody.getIntervals()) {
			Collection<SingleInterval> set=region.allBins(binResolution2);
			rtrn.addIntervals(set);
		}
		
		return rtrn;
	}
	
	private Map<SingleInterval, Double> getCounts(BarcodingDataStreaming data, Kmer nuclearBody, int binSize, Collection<SingleInterval> tripRegions) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			c=c.bin(binSize);
			Cluster body=hitsInBody(c, nuclearBody);
			
			for(SingleInterval region: c.getAllDNAIntervals()) {
				if(tripRegions.contains(region)) {
					//TODO ensure that this cluster contains nuclear body region besides ones on same chromosome
					boolean hasNonChr=hasNonChr(body, region.getReferenceName());
					if(hasNonChr) {
						double count=0;
						if(rtrn.containsKey(region)) {count=rtrn.get(region);}
						//count++;
						count+=(2.0/c.getClusterSize());
						rtrn.put(region, count);
					}
				}
			}
			counter++;
			if(counter%100000 ==0) {System.err.println(counter);}
		}
		
		data.close();
		
		
		return rtrn;
	}
	
	private boolean hasNonChr(Cluster body, String chr) {
		for(SingleInterval region: body.getAllDNAIntervals()) {
			if(!region.getReferenceName().equals(chr)) {return true;}
		}
		return false;
	}
	
	private Cluster hitsInBody(Cluster c, Kmer nuclearBody) {
		Cluster rtrn=new Cluster("body");
		
		for(SingleInterval region: c.getAllDNAIntervals()) {
			if(nuclearBody.getIntervals().contains(region)) {rtrn.addDNARead(region);}
		}
		
		return rtrn;
	}
	

	private Map<SingleInterval, Pair<Integer>> parse(File input) throws IOException {
		Map<SingleInterval, Pair<Integer>> rtrn= new TreeMap<SingleInterval, Pair<Integer>>();
		
		List<String> lines=BEDFileIO.loadLines(input.getAbsolutePath());
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			SingleInterval r=new SingleInterval(tokens[3], Integer.parseInt(tokens[4]), Integer.parseInt(tokens[4])+1);
			Pair<Integer> pair=new Pair<Integer>();
			pair.setValue1(Integer.parseInt(tokens[1]));
			pair.setValue2(Integer.parseInt(tokens[2]));
			rtrn.put(r, pair);
		}
		
		
		return rtrn;
	}

	private Map<SingleInterval, Pair<Integer>> binAndMerge(Map<SingleInterval, Pair<Integer>> vals, int binSize) {
		Map<SingleInterval, Pair<Integer>> rtrn= new TreeMap<SingleInterval, Pair<Integer>>();
		
		for(SingleInterval r: vals.keySet()) {
			SingleInterval bin=r.bin(binSize);
			Pair<Integer> pair=vals.get(r);
			Pair<Integer> newPair=new Pair<Integer>();
			newPair.setValue1(pair.getValue1());
			newPair.setValue2(pair.getValue2());
			if(rtrn.containsKey(bin)) {
				Pair<Integer> updatedPair=rtrn.get(bin);
				newPair.setValue1(pair.getValue1()+updatedPair.getValue1());
				newPair.setValue2(pair.getValue2()+updatedPair.getValue2());
			}
			rtrn.put(bin, newPair);
		}
		
		return rtrn;
	}

	private void write(String save, Map<SingleInterval, Pair<Integer>> vals, Map<SingleInterval, Double> speckleCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: vals.keySet()) {
			Pair<Integer> pair=vals.get(r);
			double speckleScore=-1;
			if(speckleCounts.containsKey(r)) {speckleScore=speckleCounts.get(r);}
			writer.write(r.getReferenceName()+"\t"+r.getReferenceStartPosition()+"\t"+r.getReferenceEndPosition()+"\t"+pair.getValue1()+"\t"+pair.getValue2()+"\t"+speckleScore+"\n");
		}
		
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		if(args.length>4) {
			File input=new File(args[0]);
			String save=args[1];
			int size=Integer.parseInt(args[2]);
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[3]));
			Kmer kmer=new Kmer();
			kmer.addIntervals(BEDFileIO.loadSingleIntervalFromFile(args[4]));
			new BinPositions(input, save, size, data, kmer);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=input \n args[1]=save \n args[2]=binsize \n args[3]=barcoding data \n args[4]=active hub";
	
}

package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.math.ScanStat;
import guttmanlab.core.math.Statistics;

public class PullRNAsInGenomicIntervals {

	public PullRNAsInGenomicIntervals(BarcodingDataStreaming data, SingleInterval region, String save) throws IOException {
		
		int binResolution=region.getGenomicLength();
		
		Map<SingleInterval, Map<String, Integer>> rnaCountsByBin=new TreeMap<SingleInterval, Map<String, Integer>>();
		
		Map<String, Integer> totalCounts=new TreeMap<String, Integer>();
		int counter=0;
		Collection<Cluster> clusters=new HashSet<Cluster>();
		while(data.hasNext()) {
			Cluster c=data.next();
			increment(c, totalCounts);
			
			//TODO for each cluster compute RNA distribution
			Map<String, Integer> rnaCounts=getRNADistribution(c);
			//get bins, 
			Collection<SingleInterval> bins=getBins(c, binResolution, region);
			//add distribution to bins
			incrementCounts(rnaCountsByBin, bins, rnaCounts);
			
			
			
			if(c.containsOverlappingDNA(region)) {
				clusters.add(c);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		for(Cluster c: clusters) {
			Collection<String> names=c.getRNANames();
			for(String name: names) {
				int count=0;
				if(counts.containsKey(name)) {
					count=counts.get(name);
				}
				count++;
				counts.put(name, count);
			}
		}
		
		//BEDFileIO.write(counts, save);
		write(save, counts, totalCounts, rnaCountsByBin);
	}
	
	
	private Collection<SingleInterval> getBins(Cluster c, int binResolution, SingleInterval region) {
		String chr=region.getReferenceName();
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		Cluster binned=c.bin(binResolution);
		for(SingleInterval r: binned.getAllDNAIntervals()) {
			if(r.getReferenceName().equals(chr) && !r.overlaps(region)) {
				rtrn.add(r);
			}
		}
		
		return rtrn;
	}


	private void incrementCounts(Map<SingleInterval, Map<String, Integer>> rnaCountsByBin, Collection<SingleInterval> bins, Map<String, Integer> rnaCounts) {
		for(SingleInterval bin: bins) {
			Map<String, Integer> rtrn=new TreeMap<String, Integer>();
			if(rnaCountsByBin.containsKey(bin)) {
				rtrn=rnaCountsByBin.get(bin);
			}
			for(String rna: rnaCounts.keySet()) {
				int count=0;
				if(rtrn.containsKey(rna)) {
					count=rtrn.get(rna);
				}
				count+=rnaCounts.get(rna);
				rtrn.put(rna, count);
			}
			rnaCountsByBin.put(bin, rtrn);
		}
		
	}


	private Map<String, Integer> getRNADistribution(Cluster c) {
		Map<String, Integer> totalCounts=new TreeMap<String, Integer>();
		Collection<String> names=c.getRNANames();
		for(String name: names) {
			int count=0;
			if(totalCounts.containsKey(name)) {
				count=totalCounts.get(name);
			}
			count++;
			totalCounts.put(name, count);
		}
		return totalCounts;
	}


	private void write(String save, Map<String, Integer> counts, Map<String, Integer> totalCounts, Map<SingleInterval, Map<String, Integer>> rnaCountsByBin) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<String, List<Double>> map=makeMap(rnaCountsByBin, totalCounts);
		
		for(String rna: counts.keySet()) {
			int count=counts.get(rna);
			int totalCount=totalCounts.get(rna);
			double average=getAverage(map, rna, rnaCountsByBin.keySet().size());
			//double stdev=getStdev(map, rna, average);
			//double p95=getPercentile(map, rna, 0.95);
			//double avg2=(double)totalCount/(double)rnaCountsByBin.keySet().size();
			
			//average=Math.max(average, avg2);
			
			double p=ScanStat.poissonPValue(count, average);
			double bonf=Math.min(1.0, p*counts.size());
			
			List<Double> scores=map.get(rna);
			double percentile=getPercentile(count, scores, rnaCountsByBin.keySet().size());
			
			writer.write(rna+"\t"+count+"\t"+totalCount+"\t"+average+"\t"+p+"\t"+bonf+"\t"+percentile+"\n");
		}
		
		writer.close();
	}

	
	
	private double getPercentile(int count, List<Double> scores, int size) {
		double lessThan=size-scores.size();
		
		for(double val: scores) {
			if(count<val) {lessThan++;}
		}
		
		return lessThan/(double)size;
	}


	private Map<String, List<Double>> makeMap(Map<SingleInterval, Map<String, Integer>> rnaCountsByBin, Map<String, Integer> totalCounts) {
		Map<String, List<Double>> rtrn=new TreeMap<String, List<Double>>();
		
		for(String rna: totalCounts.keySet()) {
			rtrn.put(rna, new ArrayList<Double>());
		}
		
		
		for(SingleInterval bin: rnaCountsByBin.keySet()) {
			//System.out.println(bin.toShortBED());
			for(String rna: rnaCountsByBin.get(bin).keySet()) {
				double score=(double)rnaCountsByBin.get(bin).get(rna);
				List<Double> list=rtrn.get(rna);
				list.add(score);
			}
		}
		
		
		return rtrn;
	}


	private double getStdev(Map<String, List<Double>> map, String rna, double avg) {
		List<Double> list=new ArrayList<Double>();
		if(map.containsKey(rna)) {
			list=map.get(rna);
		}
		
		return Statistics.sem(list, avg);
	}

	private double getAverage(Map<String, List<Double>> map, String rna, int size) {
		List<Double> list=new ArrayList<Double>();
		
		if(map.containsKey(rna)) {
			list=map.get(rna);
			double fractionZero=(double)list.size()/(double)size;
			//System.err.println(rna+" "+list.size()+" "+size+" "+fractionZero);
		}
		
		return Statistics.sum(list)/(double)size;
	}
	
	
	private double getPercentile(Map<String, List<Double>> map, String rna, double pct) {
		List<Double> list=new ArrayList<Double>();
		if(map.containsKey(rna)) {
			list=map.get(rna);
		}
		
		return Statistics.quantile(list, pct);
	}


	private Double get(Map<String, Integer> counts, String rna) {
		if(counts.containsKey(rna)) {
			return (double)counts.get(rna);
		}
		return 0.0;
	}


	private void increment(Cluster c, Map<String, Integer> totalCounts) {
		Collection<String> names=c.getRNANames();
		for(String name: names) {
			int count=0;
			if(totalCounts.containsKey(name)) {
				count=totalCounts.get(name);
			}
			count++;
			totalCounts.put(name, count);
		}
		
	}


	public static void main(String[] args) throws IOException{
		if(args.length>2){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			SingleInterval region=new SingleInterval(args[1]);
			String save=args[2];
			new PullRNAsInGenomicIntervals(data, region, save);
			
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=clusters \n args[1]=gene name \n args[2]=bin resolution \n args[3]=save \n args[4]=max cluster Size";
	
	
}

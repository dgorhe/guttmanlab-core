package guttmanlab.core.splicing.speckle;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.rnasprite.Kmer;

public class DistanceToNuclearBody {

	
	/*public DistanceToNuclearBody(BarcodingDataStreaming data, Kmer nuclearBody, String save, int binResolution) throws IOException {
		Map<SingleInterval, Double> counts=getCounts(data, nuclearBody, binResolution);
		
		//Map<SingleInterval, Double> norm=normalizeToHub(counts, nuclearBody, binResolution);
		
		BEDFileIO.writeBEDGraph(counts, save);
		//BEDFileIO.writeBEDGraph(norm, save+".norm.bedgraph");
	}*/

	
	/*private void write(String save, Map<Kmer, Map<SingleInterval, Double>> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<SingleInterval> allIntevals=new TreeSet<SingleInterval>();
		
		for(Kmer k: counts.keySet()) {
			allIntevals.addAll(counts.get(k).keySet());
		}
		
		
		for(SingleInterval region: allIntevals) {
			writer.write(region.toShortBED());
			for(Kmer k: counts.keySet()) {
				Map<SingleInterval, Double> vals=counts.get(k);
				double val=get(vals, region);
				writer.write("\t"+val);
			}
			writer.write("\n");
		}
			
		
		writer.close();
	}*/


	private static double get(Map<SingleInterval, Double> vals, SingleInterval region) {
		if(vals.containsKey(region)) {return vals.get(region);}
		return 0;
	}


	/*public DistanceToNuclearBody(BarcodingDataStreaming data, Kmer nuclearBody, String save, int binResolution, Collection<Gene> genes) throws IOException {
		Map<SingleInterval, Double> counts=getCounts(data, nuclearBody, binResolution);
		
		//Map<SingleInterval, Double> norm=normalizeToHub(counts, nuclearBody, binResolution);
		
		BEDFileIO.writeBEDGraph(counts, save+".bedgraph");
		
		Map<Annotation, Double> geneScores=scoreGenes(counts, genes, binResolution);
		write(geneScores, save+".genes.scores");
		
		//BEDFileIO.writeBEDGraph(norm, save+".norm.bedgraph");
	}*/
	
	
	/*public static Map<SingleInterval, Double> distance(BarcodingDataStreaming data, Kmer nuclearBody, int binResolution, Collection<SingleInterval> genes) throws IOException {
		Map<SingleInterval, Double> counts=getCounts(data, nuclearBody, binResolution);
		
		Map<SingleInterval, Double> geneScores=scoreRegions(counts, genes, binResolution);
		return geneScores;
	}*/
	
	
	
	private void write(Map<Gene, Double> geneScores, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene r: geneScores.keySet()) {
			writer.write(r.getName()+"\t"+r.toUCSC()+"\t"+geneScores.get(r)+"\n");
		}
		
		writer.close();
	}


	private static Map<SingleInterval, Double> scoreRegions(Map<SingleInterval, Double> counts, Collection<SingleInterval> genes, int binResolution) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval gene: genes) {
			Collection<SingleInterval> bins=gene.allBins(binResolution);
			List<Double> vals=new ArrayList<Double>();
			for(SingleInterval bin: bins) {
				if(counts.containsKey(bin)) {
					double val=counts.get(bin);
					vals.add(val);
				}
			}
			double score=Statistics.quantile(vals, 0.5);
			rtrn.put(gene, score);
		}
		
		return rtrn;
	}
	
	private static Map<Annotation, Double> scoreGenes(Map<SingleInterval, Double> counts, Collection<? extends Annotation> genes, int binResolution) {
		Map<Annotation, Double> rtrn=new TreeMap<Annotation, Double>();
		
		for(Annotation gene: genes) {
			Collection<SingleInterval> bins=gene.getSingleInterval().allBins(binResolution);
			List<Double> vals=new ArrayList<Double>();
			for(SingleInterval bin: bins) {
				if(counts.containsKey(bin)) {
					double val=counts.get(bin);
					vals.add(val);
				}
			}
			double score=Statistics.quantile(vals, 0.5);
			rtrn.put(gene, score);
		}
		
		return rtrn;
	}


	private Map<SingleInterval, Double> normalizeToHub(Map<SingleInterval, Double> counts, Kmer nuclearBody, int binSize) {
		nuclearBody=bin(nuclearBody, binSize);
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();

		List<Double> values=new ArrayList<Double>();
		for(SingleInterval region: nuclearBody.getIntervals()) {
			if(counts.containsKey(region)) {
				double val=counts.get(region);
				values.add(val);
			}
		}
		
		double denom=Statistics.quantile(values, 0.5);
		
		for(SingleInterval r: counts.keySet()) {
			double val=counts.get(r);
			double norm=val/denom;
			rtrn.put(r, norm);
		}
				
				
		return rtrn;
	}


	public DistanceToNuclearBody(BarcodingDataStreaming data, Kmer nuclearBody, String save, int binResolution, String chr, String gene) throws IOException {
		nuclearBody=bin(nuclearBody, binResolution);
		
		Collection<Cluster> clusters=data.getDNAClusters(nuclearBody, binResolution);
		clusters=filterClusters(clusters, gene);
		
		Map<SingleInterval, Double> counts=getCounts(nuclearBody, clusters);
		
		BEDFileIO.writeBEDGraph(counts, save+".bedgraph", chr);
	}
	
	public DistanceToNuclearBody(BarcodingDataStreaming data, Kmer nuclearBody, String save, int binResolution, String chr) throws IOException {
		nuclearBody=bin(nuclearBody, binResolution);
		
		Collection<Cluster> clusters=data.getDNAClusters(nuclearBody, binResolution);
		
		Map<SingleInterval, Double> counts=getCounts(nuclearBody, clusters);
		
		BEDFileIO.writeBEDGraph(counts, save+".bedgraph", chr);
	}
	
	
	public static Map<SingleInterval, Double> getDistanceToNuclearBody(BarcodingDataStreaming data, Kmer nuclearBody, int binResolution) throws IOException {
		nuclearBody=bin(nuclearBody, binResolution);
		
		Collection<Cluster> clusters=data.getDNAClusters(nuclearBody, binResolution);
		
		Map<SingleInterval, Double> counts=getCounts(nuclearBody, clusters);
		
		return counts;
	}
	
	
	private Collection<Cluster> filterClusters(Collection<Cluster> clusters, String gene) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Cluster c: clusters) {
			if(c.containsRNA(gene)) {rtrn.add(c);}
		}
		
		return rtrn;
	}


	

	private static Kmer bin(Kmer nuclearBody, int binResolution2) {
		Kmer rtrn=new Kmer();
		
		for(SingleInterval region: nuclearBody.getIntervals()) {
			Collection<SingleInterval> set=getRegions(region, binResolution2);
			rtrn.addIntervals(set);
		}
		
		return rtrn;
	}
	
	
	
	
	private static Collection<SingleInterval> getRegions(SingleInterval region, int binResolution) {
		return region.allBins(binResolution);
	}
	
	
	private static Map<SingleInterval, Double>[] getCounts(BarcodingDataStreaming data, Kmer nuclearBody, int binSize) {
		nuclearBody=bin(nuclearBody, binSize);
		
		//for(SingleInterval region: nuclearBody.getIntervals()) {System.out/groups/guttman/mguttman/human.clusters.println(region.toBED());}
		
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Map<SingleInterval, Double> input=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		double total=0;
		double totalInput=0;
		
		while(data.hasNext()) {
			Cluster c=data.next();
			c=c.bin(binSize);
			
			Cluster body=hitsInBody(c, nuclearBody);
			for(SingleInterval region: c.getAllDNAIntervals()) {
				//TODO ensure that this cluster contains nuclear body region besides ones on same chromosome
				boolean hasNonChr=hasNonChr(body, region.getReferenceName());
				if(hasNonChr) {
					double count=0;
					if(rtrn.containsKey(region)) {count=rtrn.get(region);}
					//count++;
					count+=(2.0/c.getClusterSize());
					total+=(2.0/c.getClusterSize());
					rtrn.put(region, count);
				}
				else {
					double count=0;
					if(input.containsKey(region)) {count=input.get(region);}
					count+=(2.0/c.getClusterSize());
					totalInput+=(2.0/c.getClusterSize());
					input.put(region, count);
				}
			}
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		data.close();
		
		rtrn=norm(rtrn, total);
		input=norm(input, totalInput);
		
		Map<SingleInterval, Double>[] array= new Map[2];
		array[0]=rtrn;
		array[1]=input;
		
		return array;
	}
	
	
	private static Map<SingleInterval, Double> getRawCounts(BarcodingDataStreaming data, Kmer nuclearBody, int binSize) {
		nuclearBody=bin(nuclearBody, binSize);
		
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		//Map<SingleInterval, Double> input=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		
		while(data.hasNext()) {
			Cluster c=data.next();
			c=c.bin(binSize);
			
			Cluster body=hitsInBody(c, nuclearBody);
			for(SingleInterval region: c.getAllDNAIntervals()) {
				boolean hasNonChr=hasNonChr(body, region.getReferenceName());
				if(hasNonChr) {
					double count=0;
					if(rtrn.containsKey(region)) {count=rtrn.get(region);}
					//count++;
					count+=(2.0/c.getClusterSize());
					rtrn.put(region, count);
				}
			}
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		data.close();
		
		
		return rtrn;
	}
	
	
	
	private static Map<SingleInterval, Double>[] getRawCounts(BarcodingDataStreaming data, Kmer nuclearBody, int binSize, double fraction, int numPerm) {
		nuclearBody=bin(nuclearBody, binSize);
		
		Map<SingleInterval, Double>[] rtrn=new Map[numPerm];
		
		for(int i=0; i<rtrn.length; i++) {rtrn[i]=new TreeMap<SingleInterval, Double>();}
		
		//Map<SingleInterval, Double> input=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		
		while(data.hasNext()) {
			Cluster c=data.next();
			double[] rands=getRandom(numPerm);
			c=c.bin(binSize);
			Cluster body=hitsInBody(c, nuclearBody);
			for(SingleInterval region: c.getAllDNAIntervals()) {
				boolean hasNonChr=hasNonChr(body, region.getReferenceName());
				if(hasNonChr) {
					for(int i=0; i<rands.length; i++) {
						if(rands[i]<fraction) {
							double count=0;
							if(rtrn[i].containsKey(region)) {count=rtrn[i].get(region);}
							count+=(2.0/c.getClusterSize());
							rtrn[i].put(region, count);
						}
					}
				}
			}
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		data.close();
		
		
		return rtrn;
	}
	
	
	private static double[] getRandom(int numPerm) {
		double[] rtrn=new double[numPerm];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=Math.random();
		}
		
		return rtrn;
	}


	private static Map<SingleInterval, Integer>[] sampleClusterCounts(BarcodingDataStreaming data, Kmer nuclearBody, int binSize, int numPerm) {
		nuclearBody=bin(nuclearBody, binSize);
		
		//for(SingleInterval region: nuclearBody.getIntervals()) {System.out/groups/guttman/mguttman/human.clusters.println(region.toBED());}
		
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		Map<SingleInterval, Integer> input=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		double total=0;
		double totalInput=0;
		
		while(data.hasNext()) {
			Cluster c=data.next();
			c=c.bin(binSize);
			
			Cluster body=hitsInBody(c, nuclearBody);
			for(SingleInterval region: c.getAllDNAIntervals()) {
				//TODO ensure that this cluster contains nuclear body region besides ones on same chromosome
				boolean hasNonChr=hasNonChr(body, region.getReferenceName());
				if(hasNonChr) {
					total++;
					int count=0;
					if(rtrn.containsKey(region)) {count=rtrn.get(region);}
					count++;
					//count+=(2.0/c.getClusterSize());
					//total+=(2.0/c.getClusterSize());
					rtrn.put(region, count);
				}
				else {
					totalInput++;
					int count=0;
					if(input.containsKey(region)) {count=input.get(region);}
					count++;
					//count+=(2.0/c.getClusterSize());
					//totalInput+=(2.0/c.getClusterSize());
					input.put(region, count);
				}
			}
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		data.close();
		
		double p=total/totalInput;
		System.err.println(total+" "+totalInput+" "+p);
		
		Map<SingleInterval, Integer>[] array= new Map[numPerm+1];
		array[0]=rtrn;
		
		for(int i=1; i<array.length; i++) {
			array[i]=sample(input, p);
		}
		
		return array;
	}
	
	
	private static Map<SingleInterval, Double>[] sampleClusterCounts(File sample, Collection<File> controls, double fraction, Kmer nuclearBody, int binSize, int numPerm) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(sample);
		Map<SingleInterval, Double> observedCounts=getRawCounts(data, nuclearBody, binSize);
		
		/*for(SingleInterval r: observedCounts.keySet()) {
			System.out.println(r.toBedgraph(observedCounts.get(r)));
		}*/
		
		Map<SingleInterval, Double>[] expectedCounts=getRandomCounts(controls, nuclearBody, binSize, fraction, numPerm);
		
		
		
		for(SingleInterval region: observedCounts.keySet()) {
			double observed=observedCounts.get(region);
			double[] expected=getExpected(expectedCounts, region);
			double p=Statistics.percentLessThan(observed, expected);
			double eMean=Statistics.mean(expected);
			if(p>0.95 && eMean>0) {
				double p95=Statistics.quantile(expected, 0.95);
				double diff=Math.max(0, observed-p95);
				if(diff>0) {System.out.println(region.toBedgraph(diff));}
			}
		}
		
		return expectedCounts;
	}
	
	
	/*private static Map<SingleInterval, Integer>[] getRandomCounts(Collection<File> controls, Kmer nuclearBody, int binSize, double fraction, int numPerm) throws IOException {
		Map<SingleInterval, Integer>[] rtrn=new Map[numPerm];
		
		for(int i=0; i<numPerm; i++) {
			rtrn[i]=getRandomCounts(controls, nuclearBody, binSize, fraction);
		}
		
		return rtrn;
	}*/


	private static Map<SingleInterval, Double>[] getRandomCounts(Collection<File> controls, Kmer nuclearBody, int binSize, double fraction, int numPerm) throws IOException {
		Map<SingleInterval, Double>[] scores=new Map[numPerm];
		
		for(int i=0; i<scores.length; i++) {scores[i]=new TreeMap<SingleInterval, Double>();}
		
		
		for(File control: controls) {
			System.err.println(control);
			Map<SingleInterval, Double>[] map=getRawCounts(new BarcodingDataStreaming(control), nuclearBody, binSize, fraction, numPerm);
			scores=sum(scores, map);
		}
		return scores;
	}


	private static Map<SingleInterval, Double>[] sum(Map<SingleInterval, Double>[] maps1, Map<SingleInterval, Double>[] maps2) {
		Map<SingleInterval, Double>[] rtrn=new Map[maps1.length];
		for(int i=0; i<maps1.length; i++) {
			rtrn[i]=sum(maps1[i], maps2[i]);
		}
		return rtrn;
	}
	


	private static Map<SingleInterval, Double> sum(Map<SingleInterval, Double> scores, Map<SingleInterval, Double> map) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Collection<SingleInterval> keys=new TreeSet<SingleInterval>();
		keys.addAll(scores.keySet());
		keys.addAll(map.keySet());
		
		for(SingleInterval key: keys) {
			double score1=get(scores, key);
			double score2=get(map, key);
			double sum=score1+score2;
			rtrn.put(key, sum);
		}
		return rtrn;
	}


	private static int getI(Map<SingleInterval, Integer> scores, SingleInterval key) {
		if(scores.containsKey(key)) {return scores.get(key);}
		return 0;
	}


	private static double[] getExpected(Map<SingleInterval, Double>[] expectedCounts, SingleInterval region) {
		double[] rtrn=new double[expectedCounts.length];
		
		for(int i=0; i<expectedCounts.length; i++) {
			double score=0;
			if(expectedCounts[i].containsKey(region)) {
				score=expectedCounts[i].get(region);
			}
			rtrn[i]=score;
		}
		
		return rtrn;
	}


	private static Map<SingleInterval, Integer> sample(Map<SingleInterval, Integer> input, double p) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(SingleInterval region: input.keySet()) {
			int count=input.get(region);
			int sampledCount=sample(count, p);
			rtrn.put(region, sampledCount);
		}
		
		return rtrn;
	}


	private static int sample(int count, double p) {
		int rtrn=0;
		for(int i=0; i<count; i++) {
			double rand=Math.random();
			if(rand<p) {rtrn++;}
		}
		return rtrn;
	}


	private static Map<SingleInterval, Double> getCounts(BarcodingDataStreaming data, int binSize) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		int counter=0;
		double total=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			c=c.bin(binSize);
			
			for(SingleInterval region: c.getAllDNAIntervals()) {
				double count=0;
				if(rtrn.containsKey(region)) {count=rtrn.get(region);}
				count+=(2.0/c.getClusterSize());
				total+=(2.0/c.getClusterSize());
				rtrn.put(region, count);
			}
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		
		rtrn=norm(rtrn, total);
		
		data.close();
		return rtrn;
	}
	
	
	private static Map<SingleInterval, Double> norm(Map<SingleInterval, Double> rtrn, double total) {
		Map<SingleInterval, Double> newMap=new TreeMap<SingleInterval, Double>();
		for(SingleInterval r: rtrn.keySet()) {
			double val=rtrn.get(r);
			double norm=(val/total)*10000000;
			newMap.put(r, norm);
		}
		
		return newMap;
	}


	private static Map<Kmer, Map<SingleInterval, Double>> getCounts(BarcodingDataStreaming data, Collection<Kmer> nuclearBodies, int binSize) {
		nuclearBodies=bin(nuclearBodies, binSize);
		System.err.println(nuclearBodies.size());
		
		
		Map<Kmer, Map<SingleInterval, Double>> map=new TreeMap<Kmer, Map<SingleInterval, Double>>();
		
		int counter=0;
		
		while(data.hasNext()) {
			Cluster c=data.next();
			c=c.bin(binSize);
			for(Kmer nuclearBody: nuclearBodies) {
				Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
				if(map.containsKey(nuclearBody)) {rtrn=map.get(nuclearBody);}
				
				Cluster body=hitsInBody(c, nuclearBody);
				for(SingleInterval region: c.getAllDNAIntervals()) {
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
				map.put(nuclearBody, rtrn);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		data.close();
		return map;
	}
	
	private static Collection<Kmer> bin(Collection<Kmer> nuclearBodies, int binSize) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		int count=0;
		for(Kmer n: nuclearBodies) {
			Kmer k=bin(n, binSize);
			k.setName("k"+count);
			count++;
			rtrn.add(k);
		}
		
		return rtrn;
	}


	private static Map<SingleInterval, Double> getCounts(Kmer nuclearBody, Collection<Cluster> clusters) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(Cluster c: clusters) {
			Cluster body=hitsInBody(c, nuclearBody);
			for(SingleInterval region: c.getAllDNAIntervals()) {
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
		return rtrn;
	}
	
	private static Cluster hitsInBody(Cluster c, Kmer nuclearBody) {
		Cluster rtrn=new Cluster("body");
		
		for(SingleInterval region: c.getAllDNAIntervals()) {
			if(nuclearBody.getIntervals().contains(region)) {rtrn.addDNARead(region);}
		}
		
		return rtrn;
	}
	
	private static boolean hasNonChr(Cluster body, String chr) {
		for(SingleInterval region: body.getAllDNAIntervals()) {
			if(!region.getReferenceName().equals(chr)) {return true;}
		}
		return false;
	}
	
	
	
	/*private static Map<Annotation, Double> scoreGenes(BarcodingDataStreaming data, Kmer kmer, int binResolution, Collection<Gene> genes) {
		Map<SingleInterval, Double> counts=getCounts(data, kmer, binResolution);
		
		//Map<SingleInterval, Double> norm=normalizeToHub(counts, nuclearBody, binResolution);
		
		Collection<SingleInterval> hubRegions=getHubRegions(kmer);
		
		Map<Annotation, Double> kmerScores=scoreGenes(counts, hubRegions, binResolution);
		
		Map<Annotation, Double> geneScores=scoreGenes(counts, genes, binResolution);
		//write(geneScores, save+".genes.scores");
		
		//BEDFileIO.writeBEDGraph(norm, save+".norm.bedgraph");
		
		geneScores.putAll(kmerScores);
		
		return geneScores;
	}*/
	
	
	
	public static Map<SingleInterval, Double>[] scoreBins(BarcodingDataStreaming data, Kmer kmer, int binResolution) {
		return getCounts(data, kmer, binResolution);
	}
	
	public static Map<SingleInterval, Double> scoreBins(BarcodingDataStreaming data, int binResolution) {
		return getCounts(data, binResolution);
	}
	
	
	private static double median(Map<Annotation, Double> kmerScores) {return quantile(kmerScores, 0.5);}
	
	private static double quantile(Map<Annotation, Double> kmerScores, double percentile) {
		List<Double> vals=new ArrayList<Double>();
		
		for(Annotation r: kmerScores.keySet()) {
			vals.add(kmerScores.get(r));
		}
		
		return Statistics.quantile(vals, percentile);
	}


	private static Collection<SingleInterval> getHubRegions(Kmer kmer) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		int counter=0;
		for(SingleInterval interval: kmer.getIntervals()) {
			interval.setName("hub"+counter);
			rtrn.add(interval);
			counter++;
		}
		
		return rtrn;
	}


	/*private static void write(String save, Map<Annotation, Double>[] scores, int[] binResolutions) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Annotation> genes=new TreeSet<Annotation>();
		for(int i=0; i<scores.length; i++) {
			genes.addAll(scores[i].keySet());
		}
		
		for(Annotation g: genes) {
			double[] vals=getVals(scores, g);
			writer.write(g.getName()+"\t"+g.toUCSC());
			for(int i=0; i<vals.length; i++) {writer.write("\t"+vals[i]);}
			writer.write("\n");
		}
		
		writer.close();
	}*/
	
	private static double[] getVals(Map<Annotation, Double>[] scores, Annotation g) {
		double[] rtrn=new double[scores.length];
		
		for(int i=0; i<scores.length; i++) {
			rtrn[i]=get(scores[i], g);
		}
		
		return rtrn;
	}


	private static double get(Map<Annotation, Double> map, Annotation g) {
		if(map.containsKey(g)) {return map.get(g);}
		return 0;
	}

	
	private static void write(String save, Map<SingleInterval, Integer> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: scores.keySet()) {
			writer.write(region.toBedgraph(scores.get(region))+"\n");
		}
		
		writer.close();
	}
	
	private static void write(String save, Map<SingleInterval, Integer> sample, Map<SingleInterval, Integer> input) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: sample.keySet()) {
			if(input.containsKey(region)) {
				double sampleScore=sample.get(region);
				double inputScore=input.get(region);
				double ratio=sampleScore-inputScore;
				writer.write(region.toBedgraph(ratio)+"\n");
			}
			else {System.err.println("skipped "+region.toUCSC());}
		}
		
		writer.close();
	}
	
	
	private static void write(String save, Map<SingleInterval, Integer>[] scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<SingleInterval, Integer> sample=scores[0];
		
		for(SingleInterval region: sample.keySet()) {
			double sampleScore=sample.get(region);
			double inputScore=getPercentile(scores, region, 0.95);
			double ratio=Math.max(0,sampleScore-inputScore);
			writer.write(region.toBedgraph(ratio)+"\n");
		}
		
		writer.close();
	}

	
	private static void writeInput(String save, Map<SingleInterval, Integer>[] scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<SingleInterval, Integer> sample=scores[0];
		
		for(SingleInterval region: sample.keySet()) {
			double sampleScore=sample.get(region);
			double inputScore=getPercentile(scores, region, 0.95);
			//double ratio=sampleScore-inputScore;
			writer.write(region.toBedgraph(inputScore)+"\n");
		}
		
		writer.close();
	}





	private static double getPercentile(Map<SingleInterval, Integer>[] scores, SingleInterval region, double d) {
		List<Double> vals=new ArrayList<Double>();
		for(int i=1; i<scores.length; i++) {
			double score=0;
			if(scores[i].containsKey(region)) {score=scores[i].get(region);}
			vals.add(score);
		}
		
		return Statistics.quantile(vals, d);
	}

	
	private static Collection<File> getControls(File sample, String dir) {
		Collection<File> rtrn=new ArrayList<File>();
		File[] files=new File(dir).listFiles();
		for(int i=0; i<files.length; i++) {
			if(!files[i].equals(sample)) {
				System.err.println(files[i]);
				rtrn.add(files[i]);
			}
		}
		return rtrn;
	}
	
	private static double getScore(Map<SingleInterval, Double> distances, SingleInterval bin) {
		if(distances.containsKey(bin)) {return distances.get(bin);}
		return 0;
	}


	public static void main(String[] args) throws IOException {
		if(args.length>4) {
			/*File sample=new File(args[0]);
			Collection<File> controls=getControls(sample, args[1]);
			double fraction=Double.parseDouble(args[2]);
			Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[3]);
			Kmer kmer=new Kmer();
			kmer.addIntervals(regions);
			int binResolution= Integer.parseInt(args[4]);
			int numPerm=100;
			
			Map<SingleInterval, Double>[] scores=DistanceToNuclearBody.sampleClusterCounts(sample, controls, fraction, kmer, binResolution, numPerm);*/
			
			
			BarcodingDataStreaming data1=new BarcodingDataStreaming(new File(args[0]));
			BarcodingDataStreaming data2=new BarcodingDataStreaming(new File(args[1]));
			
			
			Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[2]);
			Kmer kmer=new Kmer();
			kmer.addIntervals(regions);
			int binResolution=Integer.parseInt(args[3]);
			
			Map<SingleInterval, Double> distances1=getDistanceToNuclearBody(data1, kmer, binResolution);
			Map<SingleInterval, Double> distances2=getDistanceToNuclearBody(data2, kmer, binResolution);
			Collection<SingleInterval> genes=BEDFileIO.loadSingleIntervalFromFileFromGTF(args[4], true);
			
			for(SingleInterval gene: genes) {
				String name=gene.getName();
				Collection<SingleInterval> bins=gene.allBins(binResolution);
				for(SingleInterval bin: bins) {
					double score1=getScore(distances1, bin);
					double score2=getScore(distances2, bin);
					System.out.println(name+"\t"+bin.toUCSC()+"\t"+score1+"\t"+score2);
				}
				
			}
			
			
			
			
			System.err.println("done");
		}
		else {
			System.err.println(usage);
		}
	}
	

	


	static String usage=" args[0]=data1 \n args[1]=data2 \n args[2]=kmer \n args[3]=bin resolution \n args[4]=genes";
}

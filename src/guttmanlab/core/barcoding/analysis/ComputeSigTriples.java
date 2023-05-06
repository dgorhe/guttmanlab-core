package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.math.Statistics;

public class ComputeSigTriples {

	int numPerm=100;
	int resolution=1000000;
	int totalNumberOfRegions;
	Map<String, Integer> totalByChromosome;
	Map<String, Integer> clusterSize;
	Map<String, Map<String, Integer>> clusterSizeByChr;
	//Map<String, Cluster> clusters;
	
	public ComputeSigTriples(BarcodingDataStreaming data, Collection<Cluster> triples, String save) throws IOException{
		
		
		this.clusterSize=new TreeMap<String, Integer>();
		this.clusterSizeByChr=new TreeMap<String, Map<String, Integer>>();
		//this.clusters=new TreeMap<String, Cluster>();
		Map<SingleInterval, Collection<String>> regionsToCluster=getRegionsToCluster(data);
		
		this.totalByChromosome=new TreeMap<String, Integer>();
		this.totalNumberOfRegions=getTotalRegions(regionsToCluster);
		
		Map<SingleInterval, Integer> frequencyOfSingles=getFrequencyOfAllRegions(regionsToCluster);
		
		int count=0;
		//FileWriter writer=new FileWriter(save);
		for(Cluster triple: triples){
			//double p=score(regionsToCluster, frequencyOfSingles, triple);
			double z=zscore(regionsToCluster, frequencyOfSingles, triple, save);
			//System.out.println(triple.toStringNoName()+"\t"+p+"\t"+z);
			/*if(p<0.01){
				writer.write("c"+count+"\t"+triple.toString(false)+"\n");
				//write(save+"."+triple.toFileName()+".clusters", triple, regionsToCluster);
				count++;	
			}*/
		}
		//writer.close();
		//System.err.println(count+" "+triples.size());
	}
	
	/*private void write(String save, Cluster triple, Map<SingleInterval, Collection<String>> regionsToCluster) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<String> barcodes=this.intersect(regionsToCluster, triple.getAllIntervals());
		
		for(String barcode: barcodes){
			writer.write(this.clusters.get(barcode)+"\n");
		}
		
		writer.close();
	}*/

	private int getTotalRegions(Map<SingleInterval, Collection<String>> regionsToCluster) {
		int totalNumberOfRegions=0;
		
		for(SingleInterval region: regionsToCluster.keySet()){
			String chr=region.getReferenceName();
			int count=0;
			if(this.totalByChromosome.containsKey(chr)){count=this.totalByChromosome.get(chr);}
			count+=regionsToCluster.get(region).size();
			this.totalByChromosome.put(chr, count);
			totalNumberOfRegions+=regionsToCluster.get(region).size();
		}
		
		return totalNumberOfRegions;
	}

	private Map<SingleInterval, Collection<String>> getRegionsToCluster(BarcodingDataStreaming data) {
		Map<SingleInterval, Collection<String>> rtrn=new TreeMap<SingleInterval, Collection<String>>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			//this.clusters.put(c.getBarcode(), c);
			c=c.bin(resolution);
			for(SingleInterval region: c.getAllIntervals()){
				Collection<String> list=new TreeSet<String>();
				if(rtrn.containsKey(region)){list=rtrn.get(region);}
				list.add(c.getBarcode());
				rtrn.put(region, list);
			}
			Map<String, Integer> chrCount=getChrCount(c);
			this.clusterSizeByChr.put(c.getBarcode(), chrCount);
			this.clusterSize.put(c.getBarcode(), c.getClusterSize());
		}
		
		return rtrn;
	}

	private Map<String, Integer> getChrCount(Cluster c) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(SingleInterval region: c.getAllIntervals()){
			int count=0;
			String chr=region.getReferenceName();
			if(rtrn.containsKey(chr)){
				count=rtrn.get(chr);
			}
			count++;
			rtrn.put(chr, count);
		}
		
		return rtrn;
	}

	private Map<SingleInterval, Integer> getFrequencyOfAllRegions(Map<SingleInterval, Collection<String>> regionsToClusters) {
		
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(SingleInterval region: regionsToClusters.keySet()){
			int count=regionsToClusters.get(region).size();
			rtrn.put(region, count);
		}
		
		return rtrn;	
	}

	private int getTotal(Map<SingleInterval, Collection<String>> allRegions) {
		int total=0;
		
		for(SingleInterval region: allRegions.keySet()){
			total+=allRegions.get(region).size();
		}
		
		return total;
	}

	private double score(Map<SingleInterval, Collection<String>> regionsToBarcodes, Map<SingleInterval, Integer> frequencyOfSingles, Cluster triple) throws IOException {
		//for each triple get all clusters with AB
		//Count all C
		
		int actual=intersect(regionsToBarcodes, triple.getAllIntervals()).size();
		
		Map<SingleInterval, Collection<SingleInterval>> regions=getRegions(triple);
		
		double maxP=-1.0;
		//double maxP2=-1.0;
		double maxAvg=-1.0;
		
		for(SingleInterval region: regions.keySet()){
			double globalFreq=(double)frequencyOfSingles.get(region)/(double)this.totalNumberOfRegions; //TODO
			double localFreq=(double)frequencyOfSingles.get(region)/(double)this.totalByChromosome.get(region.getReferenceName()); //TODO
				
			Collection<String> pairClusters=intersect(regionsToBarcodes, regions.get(region));
			//For each cluster, what is the prob of having region by chance?
			
			int numGreaterThan=1;
			//int numGreaterThan2=1;
			double[] permCounts=new double[numPerm];
			for(int i=0; i<this.numPerm; i++){
				int sum=0;
				//int sum2=0;
				for(String c: pairClusters){
					double size=this.clusterSize.get(c)-regions.get(region).size();
					//double sameChrCount=count(this.clusterSizeByChr.get(c), region.getReferenceName(), regions.get(region));
					//double diffChrCount=size-sameChrCount;
					double prob=(globalFreq*size);
					//double prob2=((globalFreq*diffChrCount)+(localFreq*sameChrCount));
					//double prob2=(localFreq*sameChrCount);
					double random=Math.random();
					//System.err.println(prob+" "+prob2+" "+localFreq+" "+globalFreq);
					if(random<prob){
						sum++;
					}
					/*if(random<prob2){
						sum2++;
					}*/
				}
				if(sum>=actual){numGreaterThan++;}
				//if(sum2>=actual){numGreaterThan2++;}
				permCounts[i]=sum;
				
			}
			double p=(double)numGreaterThan/((double)numPerm+1);
			//double p2=(double)numGreaterThan2/((double)numPerm+1);
			maxP=Math.max(maxP, p);
			//maxP2=Math.max(maxP2, p2);
			double avg=Statistics.mean(permCounts);
			maxAvg=Math.max(maxAvg, avg);
				
			//System.err.println(triple.toStringNoName()+" "+region.toUCSC()+" "+numGreaterThan+" "+p+" "+maxP);
			
		}
		
		/*if(maxP<0.1 && maxP2<0.1){
			System.out.println(triple.toStringNoName()+"\t"+actual+"\t"+maxAvg+"\t"+maxP+"\t"+maxP2);
		}*/
		
		return maxP;
	}

	
	private double zscore(Map<SingleInterval, Collection<String>> regionsToBarcodes, Map<SingleInterval, Integer> frequencyOfSingles, Cluster triple, String saveDir) throws IOException {
		//for each triple get all clusters with AB
		//Count all C
		double minZ=Double.MAX_VALUE;
		double minZ2=Double.MAX_VALUE;
		int actual=intersect(regionsToBarcodes, triple.getAllIntervals()).size();
		
		Map<SingleInterval, Collection<SingleInterval>> regions=getRegions(triple);
		
		for(SingleInterval region: regions.keySet()){
			Collection<String> pairClusters=intersect(regionsToBarcodes, regions.get(region));
			FileWriter writer=new FileWriter(saveDir+"/"+filename(regions.get(region))+".bedgraph");
			//For each cluster, what is the prob of having region by chance?
			ArrayList<Double> otherScores=new ArrayList<Double>();
			ArrayList<Double> otherScoresChr=new ArrayList<Double>();
			Map<SingleInterval, Double> randomScores=new TreeMap<SingleInterval, Double>();
			
			//get all other regions
			for(SingleInterval otherRegion: regionsToBarcodes.keySet()){
				Cluster randomTriple=new Cluster("random");
				randomTriple.addRead(otherRegion);
				double randomScore=intersect(pairClusters, regionsToBarcodes.get(otherRegion)).size();
				otherScores.add(randomScore);
				randomScores.put(otherRegion, randomScore);
				if(otherRegion.getReferenceName().equalsIgnoreCase(region.getReferenceName())){otherScoresChr.add(randomScore);}
			}
			
			//double avg=Statistics.mean(otherScores);
			//double sem=Statistics.sem(otherScores, avg);
			double z=zscore(actual, otherScores);
			double z2=zscore(actual, otherScoresChr);
			minZ=Math.min(z, minZ);
			minZ2=Math.min(z2, minZ2);
			
			for(SingleInterval otherRegion: randomScores.keySet()){
				double score=randomScores.get(otherRegion);
				double z1=zscore(score, otherScoresChr);
				writer.write(otherRegion.getReferenceName()+"\t"+otherRegion.getReferenceStartPosition()+"\t"+otherRegion.getReferenceEndPosition()+"\t"+z1+"\n");
			}
			writer.close();
			
		}
		//System.err.println(triple.toStringNoName()+" "+minZ+" "+minZ2);
		return minZ;
	}
	
	
	private String filename(Collection<SingleInterval> regions) {
		String rtrn="";
		for(SingleInterval region: regions){
			rtrn+=region.getFileName();
		}
		return rtrn;
	}

	private double zscore(double actual, ArrayList<Double> otherScores) {
		double z=Statistics.zScore(actual, otherScores);
		return z;
	}

	private double count(Map<String, Integer> map, String chr, Collection<SingleInterval> list) {
		int count=0;
		if(map.containsKey(chr)){count=map.get(chr);}
		
		for(SingleInterval region: list){
			if(region.getReferenceName().equals(chr)){count=count-1;}
		}
		
		return Math.max(count, 0);
	}

	private Collection<String> intersect(Map<SingleInterval, Collection<String>> regionsToBarcodes, Collection<SingleInterval> allIntervals) {
		Collection<String> rtrn=new TreeSet<String>();
		boolean started=false;
		
		for(SingleInterval region: allIntervals){
			Collection<String> temp=regionsToBarcodes.get(region);
			rtrn=intersect(temp, rtrn, started);
			started=true;
		}
		return rtrn;
	}

	
	private Collection<String> intersect(Collection<String> list1, Collection<String> list2) {
		return intersect(list1, list2, true);
	}
	private Collection<String> intersect(Collection<String> list1, Collection<String> list2, boolean started) {
		if(!started){return list1;}
		Collection<String> rtrn=new TreeSet<String>();
		Collection<String> smaller=getSmaller(list1, list2);
		Collection<String> larger=getLarger(list1, list2);
		
		for(String barcode: smaller){
			if(larger.contains(barcode)){rtrn.add(barcode);}
		}
		
		return rtrn;
	}

	private Collection<String> getSmaller(Collection<String> list1, Collection<String> list2) {
		if(list1.size()<list2.size()){return list1;}
		return list2;
	}
	
	private Collection<String> getLarger(Collection<String> list1, Collection<String> list2) {
		if(list1.size()<list2.size()){return list2;}
		return list1;
	}

	private Map<SingleInterval, Collection<SingleInterval>> getRegions(Cluster triple) {
		Map<SingleInterval, Collection<SingleInterval>> rtrn=new TreeMap<SingleInterval, Collection<SingleInterval>>();
		for(SingleInterval region: triple.getAllIntervals()){
			Collection<SingleInterval> remainder=new TreeSet<SingleInterval>();
			remainder.addAll(triple.getAllIntervals());
			remainder.remove(region);
			rtrn.put(region, remainder);
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		Collection<Cluster> triples=new BarcodingData(new File(args[1])).getUniqueClusters();
		String save=args[2];
		System.err.println(triples.size());
		
		/*Collection<Cluster> rtrn=new TreeSet<Cluster>();
		Cluster triple=triples.iterator().next();
		Collection<Cluster> clusters=data.getClustersOverlappingRegion(triple.getAllIntervals());
		for(Cluster c: clusters){
			boolean overlapsAll=true;
			for(SingleInterval region: triple.getAllIntervals()){
				if(!c.overlapsInterval(region)){overlapsAll=false;}
			}
			if(overlapsAll){
				System.out.println(c.getBarcode()+"\t"+c.getClusterSize());
				rtrn.add(c);}
		}
		
		System.err.println(rtrn.size());*/
		
		
		new ComputeSigTriples(data, triples, save);
	}
	
}

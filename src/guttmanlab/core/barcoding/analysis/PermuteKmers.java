package guttmanlab.core.barcoding.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.Statistics;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Permute the position of a kmer in genome by size
 * @author mguttman
 *
 */
public class PermuteKmers {
	int maxClusterSize;
	int resolution=1;
	Map<String, IntervalTree<File>> barcodeFileTree;
	Map<String, Integer> chrSizes;
	Collection<Annotation> regionsToExclude;
	
	
	public PermuteKmers(Map<String, IntervalTree<File>> barcodeFileTree, Map<String, Integer> chrSizes, int maxClusterSize) throws IOException{
		this.maxClusterSize=maxClusterSize;
		this.barcodeFileTree=barcodeFileTree;
		this.chrSizes=chrSizes;
		this.regionsToExclude=new TreeSet<>();
	}
	
	public PermuteKmers(Map<String, IntervalTree<File>> barcodeFileTree, Map<String, Integer> chrSizes, int maxClusterSize, int resolution) throws IOException{
		this.maxClusterSize=maxClusterSize;
		this.barcodeFileTree=barcodeFileTree;
		this.chrSizes=chrSizes;
		this.regionsToExclude=new TreeSet<>();
		this.resolution=resolution;
	}
	
	

	/**
	 * Given a set of regions, count the unique number of any k-mer>=k
	 * @param regions the list of regions (i.e. NORs)
	 * @param k the size
	 * @param numPerm (number of permutations to perform)
	 * @return
	 * @throws IOException 
	 */
	public Map<Cluster, Double> getRandomizedCountsForAllUniqueObservations(Collection<SingleInterval> regions, int k, int numPerm, String save, boolean onlyInterchromosomal) throws IOException {
		Map<Cluster, Double> rtrn=new TreeMap<Cluster, Double>();
		FileWriter writer=new FileWriter(save);
		
		Cluster c=new Cluster("Observed", regions);
		System.err.println("Observed\t"+c.toStringNoName());
		
		double observed=countUniqueSPRITEClusters(regions, k, onlyInterchromosomal); //countUniqueSPRITEClusters
		rtrn.put(c, observed);
		writer.write("Observed\t"+c.toStringNoName()+"\t"+observed+"\n");
		System.err.println("Observed\t"+c.toStringNoName()+"\t"+observed);
		
		for(int i=0; i<numPerm; i++){
			//Cluster randomCluster=randomCluster(c, chrSizes);
			Cluster randomCluster=c.getPermutedCluster(chrSizes, resolution);
			double randomScore=countUniqueSPRITEClusters(randomCluster.getAllIntervals(), k, onlyInterchromosomal);
			rtrn.put(randomCluster, randomScore);
			writer.write("Random\t"+randomCluster.toStringNoName()+"\t"+randomScore+"\n");
			System.err.println("Random"+i+"\t"+randomCluster.toStringNoName()+"\t"+randomScore);
		}
		writer.close();
		return rtrn;
	}
	
	public void getRandomizedCountsForMaxKmer(Collection<SingleInterval> regions, int k, int numPerm, String save, boolean onlyInterchromosomal) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Cluster c=new Cluster("Observed", regions);
		System.err.println("Observed\t"+c.toStringNoName());
		
		Cluster maxObserved=countSPRITEClustersForMaxKmer(regions, k, onlyInterchromosomal); //countUniqueSPRITEClusters
		
		writer.write("Observed\t"+c.toStringNoName()+"\t"+maxObserved.getScore()+"\n");
		System.err.println("Observed\t"+c.toStringNoName()+"\t"+maxObserved.toStringNoName()+"\t"+maxObserved.getScore());
		
		for(int i=0; i<numPerm; i++){
			//Cluster randomCluster=randomCluster(c, chrSizes);
			Cluster randomCluster=maxObserved.getPermutedCluster(chrSizes, resolution);
			double randomScore=weightedQuantify(randomCluster);
			//double randomScore=countSPRITEClustersForMaxKmer(randomCluster.getAllIntervals(), k, onlyInterchromosomal);
			writer.write("Random\t"+randomCluster.toStringNoName()+"\t"+randomScore+"\n");
			System.err.println("Random"+i+"\t"+randomCluster.toStringNoName()+"\t"+randomScore);
		}
		writer.close();
	}
	
	public Map<Cluster, Double> getRandomizedCountsForAllUniqueObservations(Collection<SingleInterval> regions, int k, int numPerm, String save) throws IOException {
		return getRandomizedCountsForAllUniqueObservations(regions, k, numPerm, save, false);
	}
	
	private double countUniqueSPRITEClusters(Collection<SingleInterval> regions, int k) throws IOException {
		return countUniqueSPRITEClusters(regions, k, false);
	}
	
	private double countUniqueSPRITEClusters(Collection<SingleInterval> regions, int k, boolean onlyInterChromosomal) throws IOException {
		//TODO We can load all data clusters having regions and score off of these
		BarcodingData data=getData(regions);
		System.err.println("loaded data");
		Collection<Cluster> clusters=getUniqueSPRITEClusters(regions, k, onlyInterChromosomal);
		double count=getAllUniqueSPRITEClusters(data, clusters).size();
		return count;
	}
	
	private Cluster countSPRITEClustersForMaxKmer(Collection<SingleInterval> regions, int k, boolean onlyInterChromosomal) throws IOException {
		//TODO We can load all data clusters having regions and score off of these
		BarcodingData data=getData(regions);
		System.err.println("loaded data");
		Collection<Cluster> clusters=getUniqueSPRITEClusters(regions, k, onlyInterChromosomal);
		return getMaxKmerScore(data, clusters);
	}

	
	Collection<Cluster> getAllUniqueSPRITEClusters(Collection<SingleInterval> regions, int k, boolean onlyInterChromosomal) throws IOException {
		BarcodingData data=getData(regions);
		System.err.println("loaded data");
		Collection<Cluster> clusters=getUniqueSPRITEClusters(regions, k, onlyInterChromosomal);
		return getAllUniqueSPRITEClusters(data, clusters);
	}
	
	Collection<Cluster> getAllUniqueSPRITEClusters(Collection<SingleInterval> regions) throws IOException {
		BarcodingData data=getData(regions);
		System.err.println("loaded data");
		//Collection<Cluster> clusters=getUniqueSPRITEClusters(regions);
		
		Cluster cluster=new Cluster("c");
		cluster.addReads(regions);
		
		return getAllUniqueSPRITEClusters(data, cluster);
	}
	
	Collection<Cluster> getAllUniqueSPRITEClusters(SingleInterval region) throws IOException {
		BarcodingData data=getData(region);
		System.err.println("loaded data");
		
		Cluster cluster=new Cluster("c");
		cluster.addRead(region);
		
		return getAllUniqueSPRITEClusters(data, cluster);
	}
	
	private BarcodingData getData(SingleInterval region) throws IOException {
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		regions.add(region);
		
		BarcodingData data=loadData(this.barcodeFileTree, regions, resolution);
		return data;
	}
	
	private BarcodingData getData(Collection<SingleInterval> regions) throws IOException {
		BarcodingData data=loadData(this.barcodeFileTree, regions, resolution);
		return data;
	}
	
	private BarcodingData getData(Collection<SingleInterval> regions, int resolution) throws IOException {
		BarcodingData data=loadData(this.barcodeFileTree, regions, resolution);
		return data;
	}
	
	private Collection<Cluster> getAllUniqueSPRITEClusters(BarcodingData data, Cluster cluster) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		Collection<Cluster> allSpriteMatches=quantify(data, cluster);
		System.err.println(cluster.toStringNoName()+" "+allSpriteMatches.size());
		rtrn.addAll(allSpriteMatches);
		
		return rtrn;
	}

	private Collection<Cluster> getAllUniqueSPRITEClusters(BarcodingData data, Collection<Cluster> clusters) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(Cluster c: clusters){
			Collection<Cluster> allSpriteMatches=quantify(data, c);
			System.err.println(c.toStringNoName()+" "+allSpriteMatches.size());
			rtrn.addAll(allSpriteMatches);
		}
		return rtrn;
	}
	
	private Cluster getMaxKmerScore(BarcodingData data, Collection<Cluster> clusters) throws IOException {
		double max=0;
		Cluster maxCluster=null;
		
		for(Cluster c: clusters){
			double score=weightedQuantify(c);
			System.err.println(c.toStringNoName()+" "+score+" max so far "+max);
			if(score>max){
				maxCluster=c;
				maxCluster.setScore(score);
				max=score;
			}
		}
		return maxCluster;
	}
	
	

	Collection<Cluster> getUniqueSPRITEClusters(Collection<SingleInterval> regions, int k, boolean onlyInter){
		Collection<Cluster> clusters=new TreeSet<Cluster>();
		for(SingleInterval region1: regions){
			Cluster cluster=new Cluster(region1.getName());
			cluster.addRead(region1);
			clusters.add(cluster);
			for(int i=2; i<=k; i++){
				clusters=add(regions, clusters, i);
			}
		}
		
		clusters=getClusters(clusters, k);
		
		if(onlyInter){
			//Remove clusters that are intra
			clusters=removeIntrachromosomalClusters(clusters);
		}
		
		return clusters;
	}
	
	Collection<Cluster> getUniqueSPRITEClusters(Collection<SingleInterval> regions){
		Collection<Cluster> clusters=new TreeSet<Cluster>();
		int k=regions.size();
		
		//For each region
		for(SingleInterval region1: regions){
			Cluster cluster=new Cluster(region1.getName());
			cluster.addRead(region1);
			clusters.add(cluster);
			for(int i=2; i<=k; i++){
				clusters=add(regions, clusters, i);
			}
		}
		
		clusters=getClusters(clusters, k);
		
		
		
		return clusters;
	}
	
	private Collection<Cluster> getClusters(Collection<Cluster> clusters, int k) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Cluster c: clusters){
			if(c.getAllIntervals().size() ==k){rtrn.add(c);}
		}
		
		return rtrn;
	}

	private Collection<Cluster> removeIntrachromosomalClusters(Collection<Cluster> clusters) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Cluster c: clusters){
			Set<String> chromosomes=new TreeSet<String>();
			for(SingleInterval region: c.getAllIntervals()){chromosomes.add(region.getReferenceName());}
			if(chromosomes.size() == c.getAllIntervals().size()){
				rtrn.add(c);
			}
			
		}
		
		return rtrn;
	}

	private Collection<Cluster> add(Collection<SingleInterval> intervals, Collection<Cluster> clusters, int k) {
		//Iterate through each cluster and add each interval
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Cluster cluster: clusters){
			for(SingleInterval interval: intervals){
				Cluster newCluster=new Cluster(cluster);
				newCluster.addRead(interval);
				if(newCluster.getSize()==k){
					rtrn.add(newCluster);
				}
			}
		}
		
		return rtrn;
	}
	

	public PermuteKmers(Map<String, IntervalTree<File>> barcodeFileTree, BarcodingDataStreaming clusters, Map<String, Integer> chrSizes, String save, int maxClusterSize, int numPerm, int resolution) throws IOException{
		this.maxClusterSize=maxClusterSize;
		
		FileWriter writer=new FileWriter(save);
		writer.write("cluster name"+"\t"+"Observed count\tO/E"+"\t Average \t Max \t"+"Percent greater than random\n");
		
		
		Map<Distance, double[]> randomScores=new TreeMap<Distance, double[]>();
		
		int counter=0;
		while(clusters.hasNext()){
			Cluster cluster=clusters.next();
			/*if(!cluster.isInterchromosomal()){ //TODO This is more efficient for many of same size
				Distance d=getDistance(cluster);//TODO Extrapolate Distance to include multiple chromosomes
				this.resolution=d.getResolution(); //TODO This is not safe
				double[] scores;
				if(randomScores.containsKey(d)){
					scores=randomScores.get(d);
				}
				else{
					scores=random(d, numPerm, chrSizes, barcodeFileTree);
					randomScores.put(d, scores);
				}
				
				//double observed=quantify(barcodeFileTree, cluster).size();
				double observed=getObservedFromCluster(cluster);
				double mean=Statistics.mean(scores);
				double percentLessThan=Statistics.percentLessThan(observed, scores);
				double max=Statistics.max(scores);
				double enrichment=observed/mean;
				writer.write(cluster.toStringNoName()+"\t"+observed+"\t"+enrichment+"\t"+mean+"\t"+max+"\t"+percentLessThan+"\n");
				counter++;
				if(counter%1 ==0){System.err.println(counter +" "+observed+" "+mean +" "+randomScores.size());}
			}
			else{
				this.resolution=1000000;
				Cluster random=randomCluster(cluster, chrSizes);
				System.err.println("inter "+cluster.getReadString()+" "+cluster.toStringNoName()+" \n Random: "+random.toStringNoName());
			}*/
				
			
			/****
			 * Delete once merged
			 */
			this.resolution=1000000;
			double[] scores=random(cluster, numPerm, chrSizes, barcodeFileTree);
			double observed=getObservedFromCluster(cluster);
			double mean=Statistics.mean(scores);
			double percentLessThan=Statistics.percentLessThan(observed, scores);
			double max=Statistics.max(scores);
			double enrichment=observed/mean;
			writer.write(cluster.toStringNoName()+"\t"+observed+"\t"+enrichment+"\t"+mean+"\t"+max+"\t"+percentLessThan+"\n");
			counter++;
			if(counter%1 ==0){System.err.println(cluster.toStringNoName() +" "+observed+" "+mean+" "+max);}
		}
		
		clusters.close();
		writer.close();
	}
	
	void computePermutationsForSpecificRegion(BarcodingDataStreaming clusters, int numPerm, String save, int resolution) throws IOException{
		FileWriter writer=new FileWriter(save);
		writer.write("cluster name"+"\t"+"Observed count\tO/E"+"\t Average \t Max \t"+"Percent greater than random\n");

		this.resolution=resolution;
		
		int counter=0;
		while(clusters.hasNext()){
			Cluster cluster=clusters.next();
			double[] scores=random(cluster, numPerm, chrSizes, barcodeFileTree);
			double observed=getObservedFromCluster(cluster);
			double mean=Statistics.mean(scores);
			double percentLessThan=Statistics.percentLessThan(observed, scores);
			double max=Statistics.max(scores);
			double enrichment=observed/mean;
			writer.write(cluster.toStringNoName()+"\t"+observed+"\t"+enrichment+"\t"+mean+"\t"+max+"\t"+percentLessThan+"\n");
			counter++;
			if(counter%1 ==0){System.err.println(cluster.toString() +" "+observed+" "+mean+" "+max);}
		}
		
		clusters.close();
		writer.close();
	}
	
	
	/**
	 * Go through each cluster and compute the intrachromosomal distance, then generate permutations and write out observed and expected
	 * @param clusters
	 * @param numPerm
	 * @param save
	 * @param minCount 
	 * @throws IOException
	 */
	void computePermutationsByDistance(CloseableIterator<Cluster> clusters, int numPerm, String save, int minCount) throws IOException{
		FileWriter writer=new FileWriter(save);
		writer.write("cluster name"+"\t"+"Observed count\tNormalized Observed\tO/E"+"\t Average \t Max \t"+"Percent greater than random\n");
		
		
		Map<Distance, double[]> randomScores=new TreeMap<Distance, double[]>();
		
		int counter=0;
		while(clusters.hasNext()){
			Cluster cluster=clusters.next();
			double observedCount=getObservedFromCluster(cluster);
			if(observedCount>=minCount && !cluster.isInterchromosomal()){
				Distance d=getDistance(cluster);//TODO Extrapolate Distance to include multiple chromosomes
				this.resolution=d.getResolution(); //TODO This is not safe
				double[] scores;
				if(randomScores.containsKey(d)){
					scores=randomScores.get(d);
				}
				else{
					scores=randomWeighted(d, numPerm, chrSizes, barcodeFileTree);
					randomScores.put(d, scores);
				}
				double observed=weightedQuantify(cluster);
				double mean=Statistics.mean(scores);
				double percentLessThan=Statistics.percentLessThan(observed, scores);
				double max=Statistics.max(scores);
				double enrichment=observed/mean;
				writer.write(cluster.toStringNoName()+"\t"+observedCount+"\t"+observed+"\t"+enrichment+"\t"+mean+"\t"+max+"\t"+percentLessThan+"\n");
				counter++;
				if(counter%1 ==0){System.err.println(counter +" "+observed+" "+mean +" "+randomScores.size());}
			}
			/*else{
				this.resolution=1000000;
				Cluster random=randomCluster(cluster, chrSizes);
				System.err.println("inter "+cluster.getReadString()+" "+cluster.toStringNoName()+" \n Random: "+random.toStringNoName());
			}*/
		}
		
		clusters.close();
		writer.close();
	}

	
	private double getObservedFromCluster(Cluster cluster) throws IOException {
		String[] tokens=cluster.getReadString().split("\t");
		double observed;
		if(tokens[tokens.length-1].contains("chr")){
			observed=quantify(cluster).size();
		}
		else{
			observed=new Double(tokens[tokens.length-1]);
		}
		return observed;
	}


	private double[] random(Distance d, int numPerm, Map<String, Integer> chrSizes, Map<String, IntervalTree<File>> barcodeFileTree) throws IOException {
		double[] rtrn=new double[numPerm];
		for(int i=0; i<numPerm; i++){
			rtrn[i]=random(d, chrSizes, barcodeFileTree);
		}
		return rtrn;
	}
	
	private double[] randomWeighted(Distance d, int numPerm, Map<String, Integer> chrSizes, Map<String, IntervalTree<File>> barcodeFileTree) throws IOException {
		double[] rtrn=new double[numPerm];
		for(int i=0; i<numPerm; i++){
			rtrn[i]=randomWeighted(d, chrSizes, barcodeFileTree);
		}
		return rtrn;
	}
	
	private double[] random(Cluster c, int numPerm, Map<String, Integer> chrSizes, Map<String, IntervalTree<File>> barcodeFileTree) throws IOException {
		double[] rtrn=new double[numPerm];
		for(int i=0; i<numPerm; i++){
			rtrn[i]=random(c, chrSizes, barcodeFileTree);
		}
		return rtrn;
	}
	
	double[] randomWeighted(Cluster c, int numPerm, Map<String, Integer> chrSizes, Map<String, IntervalTree<File>> barcodeFileTree) throws IOException {
		double[] rtrn=new double[numPerm];
		for(int i=0; i<numPerm; i++){
			rtrn[i]=randomWeighted(c, chrSizes, barcodeFileTree);
		}
		return rtrn;
	}

	private double randomWeighted(Cluster cluster, Map<String, Integer> chrSizes, Map<String, IntervalTree<File>> barcodeFileTree) throws IOException {
		Cluster c=randomCluster(cluster, chrSizes);
		double random=weightedQuantify(c);
		System.err.println("Random: "+c.toStringNoName()+" "+random);
		return random;
	}

	private double random(Cluster cluster, Map<String, Integer> chrSizes, Map<String, IntervalTree<File>> barcodeFileTree) throws IOException {
		Cluster c=randomCluster(cluster, chrSizes);
		double random=quantify(c).size();
		return random;
	}
	
	private double random(Distance d, Map<String, Integer> chrSizes, Map<String, IntervalTree<File>> barcodeFileTree) throws IOException {
		Cluster c=randomCluster(d, chrSizes);
		double random=quantify(c).size();
		return random;
	}

	private double randomWeighted(Distance d, Map<String, Integer> chrSizes, Map<String, IntervalTree<File>> barcodeFileTree) throws IOException {
		Cluster c=randomCluster(d, chrSizes);
		double random=weightedQuantify(c);
		return random;
	}
	

	void writeIGV(String save, Collection<Cluster> clusters, Collection<SingleInterval> regions, int resolution) throws IOException {
		BarcodingData data=getData(regions, 1);
		Collection<Cluster> rawClusters=getRaw(data, clusters);
		
		writeIGV(save, rawClusters, resolution);
	}
	
	void writeIGV(String save, Collection<Cluster> rawClusters, int resolution) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("chromosome\tstart\tend\tfeature");
		
		Collection<SingleInterval> allReads=new TreeSet<SingleInterval>();
		
		for(Cluster c1: rawClusters){
			Cluster c=c1.bin(resolution);
			allReads.addAll(c.getAllIntervals());
			writer.write("\t"+c.getBarcode());
		}
		writer.write("\n");
		
		for(SingleInterval read: allReads){
			writer.write(read.getReferenceName()+"\t"+read.getReferenceStartPosition()+"\t"+(read.getReferenceEndPosition())+"\t"+read.toUCSC());
			for(Cluster c1: rawClusters){
				Cluster c=c1.bin(resolution);
				int score=0;
				if(c.containsInterval(read)){score=1;}
				writer.write("\t"+score);
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	private Collection<Cluster> getRaw(BarcodingData data, Collection<Cluster> clusters) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(Cluster c: clusters){
			String barcode=c.getBarcode();
			Cluster rawCluster=data.getClusterWithBarcode(barcode);
			rtrn.add(rawCluster);
		}
		return rtrn;
	}

	private void write(String save, Collection<Cluster> observedClusters) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Cluster c: observedClusters){
			writer.write(c.toSPRITEFormat()+"\n");
		}
		
		writer.close();
	}

	
	private Distance getDistance(Cluster cluster) {
		int[] distanceList=new int[cluster.getAllIntervals().size()-1];
		Iterator<SingleInterval> iter=cluster.getAllIntervals().iterator();
		SingleInterval current=iter.next();
		Collection<Integer> resolutions=new TreeSet<Integer>();
		resolutions.add(current.getLength());
		int counter=0;
		while(iter.hasNext()){
			SingleInterval next=iter.next();
			resolutions.add(next.getLength());
			int distance=next.getReferenceStartPosition()-current.getReferenceEndPosition();
			distanceList[counter]=distance;
			current=next;
			counter++;
		}
		if(resolutions.size()!=1){throw new IllegalStateException ("Kmers are not of equal resolution");}
		int resolution=resolutions.iterator().next();
		Distance d=new Distance(distanceList, resolution);
		return d;
	}
	

	private Map<Distance, Collection<Cluster>> getDistances(Collection<Cluster> clusters) {
		Map<Distance, Collection<Cluster>> rtrn=new TreeMap<Distance, Collection<Cluster>>();
		
		//Only if cluster is intrachromosomal
		for(Cluster cluster: clusters){
			if(!cluster.isInterchromosomal()){
				int[] distanceList=new int[cluster.getAllIntervals().size()-1];
				Iterator<SingleInterval> iter=cluster.getAllIntervals().iterator();
				SingleInterval current=iter.next();
				Collection<Integer> resolutions=new TreeSet<Integer>();
				resolutions.add(current.getLength());
				int counter=0;
				while(iter.hasNext()){
					SingleInterval next=iter.next();
					resolutions.add(next.getLength());
					int distance=next.getReferenceStartPosition()-current.getReferenceEndPosition();
					distanceList[counter]=distance;
					current=next;
					counter++;
				}
				if(resolutions.size()!=1){throw new IllegalStateException ("Kmers are not of equal resolution");}
				int resolution=resolutions.iterator().next();
				Distance d=new Distance(distanceList, resolution);
				Collection<Cluster> list=new TreeSet<Cluster>();
				if(rtrn.containsKey(d)){
					list=rtrn.get(d);
				}
				list.add(cluster);
				rtrn.put(d, list);
			}
		}
		return rtrn;
	}


	Collection<Cluster> quantify(Cluster c) throws IOException {
		SingleInterval current=c.getAllIntervals().iterator().next();
		BarcodingData data=loadData(this.barcodeFileTree, current, resolution);
		Collection<Cluster> clusters=data.getClustersOverlappingRegion(current);
		
		for(SingleInterval region: c.getAllIntervals()){
			clusters=overlaps(clusters, region);
		}
		
		return clusters;
	}
	
	double weightedQuantify(Cluster c) throws IOException {
		SingleInterval current=c.getAllIntervals().iterator().next();
		BarcodingData data=loadData(this.barcodeFileTree, current, resolution);
		Collection<Cluster> clusters=data.getClustersOverlappingRegion(current);
		
		for(SingleInterval region: c.getAllIntervals()){
			clusters=overlaps(clusters, region);
		}
		
		double score=weight(clusters);
		return score;
	}
	
	private double weight(Collection<Cluster> clusters) {
		double sum=0;
		for(Cluster c: clusters){
			double n=c.getAllIntervals().size();
			sum+=(1.0/n);
		}
		return sum;
	}

	void writeClusters(BarcodingDataStreaming clusters, String save) throws IOException{
		while(clusters.hasNext()){
			Cluster c=clusters.next();
			Collection<Cluster> rawClusters=quantify(c);
			writeIGV(save+"."+c.getBarcode()+".igv", rawClusters, maxClusterSize);
		}
	}
	
	
	void writeClusters(Collection<Cluster> clusters, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		for(Cluster c: clusters){
			writer.write(c.toSPRITEFormat()+"\n");
		}
		writer.close();
	}
	
	private Collection<Cluster> quantify(BarcodingData data, Cluster c) throws IOException {
		Collection<Cluster> clusters=data.getClusters();
		System.err.println("Overlap 1: "+clusters.size());
		
		for(SingleInterval region: c.getAllIntervals()){
			clusters=overlaps(clusters, region);
			System.err.println(region.toUCSC()+" "+clusters.size());
		}
		
		return clusters;
	}

	private BarcodingData loadData(Map<String, IntervalTree<File>> barcodingDataFiles, SingleInterval region, int resolution) throws IOException {
		
		if(barcodingDataFiles.containsKey(region.getReferenceName())){
			Iterator<File> barcodingData=barcodingDataFiles.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			BarcodingData data=new BarcodingData(barcodingData);
			data=data.bin(resolution, 0, maxClusterSize);
			//System.err.println("loaded data for "+region.toUCSC());
			return data;
		}
		
		return new BarcodingData();
	}
	
	private BarcodingData loadData(Map<String, IntervalTree<File>> barcodingDataFiles, Collection<SingleInterval> regions, int resolution) throws IOException {
		
		Collection<File> files=new TreeSet<File>();
		BarcodingData data=new BarcodingData();
		
		for(SingleInterval region: regions){
			if(barcodingDataFiles.containsKey(region.getReferenceName())){
				Iterator<File> iter=barcodingDataFiles.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
				while(iter.hasNext()){files.add(iter.next());}
			}
		}
		
		data=new BarcodingData(files.iterator());
		data=data.bin(resolution, 0, maxClusterSize);
		return data;
	}
	
	private Collection<Cluster> overlaps(Collection<Cluster> clusters, SingleInterval region) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(Cluster cluster: clusters){
			if(overlaps(cluster, region)){rtrn.add(cluster);}
		}
		return rtrn;
	}
	
	private boolean overlaps(Cluster cluster, SingleInterval region) {
		for(SingleInterval interval: cluster.getAllIntervals()){
			if(interval.overlaps(region)){
				return true;
			}
		}
		return false;
	}

	/**
	 * Randomize the structure of the cluster
	 * @param cluster
	 * @param chrSizes 
	 */
	private Cluster randomIntraCluster(Cluster cluster, Map<String, Integer> chrSizes) {
		//Only if cluster is intrachromosomal
		ArrayList<Integer> distanceList=new ArrayList<Integer>();
		Iterator<SingleInterval> iter=cluster.getAllIntervals().iterator();
		SingleInterval current=iter.next();
		int clusterStart=current.getReferenceStartPosition();
		int clusterEnd=current.getReferenceEndPosition();
		String chr=current.getReferenceName();
		while(iter.hasNext()){
			SingleInterval next=iter.next();
			int distance=next.getReferenceStartPosition()-current.getReferenceEndPosition();
			distanceList.add(distance);
			current=next;
			clusterEnd=Math.max(clusterEnd, current.getReferenceEndPosition());
		}
		
		
		Collection<SingleInterval> clusterRegions=new TreeSet<SingleInterval>();
		
		Iterator<SingleInterval> clusters=cluster.getAllIntervals().iterator();
		
		//Randomly pick a region in the genome
		SingleInterval currentRegion=randomRegion((clusterEnd-clusterStart), chrSizes, clusters.next().getLength());
		clusterRegions.add(currentRegion);
		for(Integer distance: distanceList){
			int clusterLength= clusters.next().getLength();
			int newStart=currentRegion.getReferenceEndPosition()+distance;
			int newEnd=newStart+clusterLength;
			currentRegion=new SingleInterval(currentRegion.getReferenceName(), newStart, newEnd);
			clusterRegions.add(currentRegion);
		}
		
		
		
		Cluster c=new Cluster("random");
		c.addReads(clusterRegions);
		//System.err.println("before "+c.toString(false));
		c=c.bin(resolution);
		//System.err.println("after "+c.toString(false));
		return c;
	}
	
	/**
	 * Randomize the structure of the cluster
	 * @param Distance
	 * @param chrSizes 
	 */
	private Cluster randomCluster(Distance d, Map<String, Integer> chrSizes) {
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		int totalGenomicLength=d.getTotalGenomeLength();
		
		//Randomly pick a region in the genome
		SingleInterval currentRegion=randomRegion(totalGenomicLength, chrSizes, d.getResolution());
		regions.add(currentRegion);
		
		int[] distances=d.getPermutedDistanceList();
		for(Integer distance: distances){
			int newStart=currentRegion.getReferenceEndPosition()+distance;
			int newEnd=newStart+d.getResolution();
			currentRegion=new SingleInterval(currentRegion.getReferenceName(), newStart, newEnd);
			regions.add(currentRegion);
		}

		Cluster c=new Cluster("random");
		c.addReads(regions);
		c=c.bin(resolution);
		return c;
	}
	
	/**
	 * Randomize the structure of the cluster -- including interchromosomal
	 * @param Distance
	 * @param chrSizes 
	 */
	private Cluster randomCluster(Cluster c, Map<String, Integer> chrSizes) {
		return c.getPermutedCluster(chrSizes, resolution);
		
		/*Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		
		Map<String, Cluster> subClustersByChr=getSubclusters(c);
		
		Cluster random=new Cluster("random");
		
		for(String chr: subClustersByChr.keySet()){
			//for each subcluster get random
			Cluster sub=subClustersByChr.get(chr);
			Cluster randomSub=randomIntraCluster(sub, chrSizes);
			random.addReads(randomSub.getAllIntervals());
		}
		
		return random;*/
	}

	private Map<String, Cluster> getSubclusters(Cluster c) {
		Map<String, Cluster> rtrn=new TreeMap<String, Cluster>();
		
		for(SingleInterval interval: c.getAllIntervals()){
			String chr=interval.getReferenceName();
			Cluster temp=new Cluster(chr);
			if(rtrn.containsKey(chr)){
				temp=rtrn.get(chr);
			}
			temp.addRead(interval);
			rtrn.put(chr, temp);
		}
		
		return rtrn;
	}


	private SingleInterval randomRegion(int totalClusterLength, Map<String, Integer> chrSizes, int clusterLength) {
		ArrayList<String> chromosomes=new ArrayList<String>();
		
		for(String chr: chrSizes.keySet()){
			int size=chrSizes.get(chr);
			if(size>totalClusterLength){chromosomes.add(chr);}
		}
	
		int random=new Double(Math.random()*chromosomes.size()).intValue();
		String chr= chromosomes.get(random);
		return randomRegion(totalClusterLength, chr, chrSizes.get(chr), clusterLength);
	}

	private SingleInterval randomRegion(int clusterGenomeLength, String chr, int chrSize, int resolution) {
		//Pick a random region
		int sizeToSelectFrom=chrSize-clusterGenomeLength;
		int startPosition=new Double(Math.random()*sizeToSelectFrom).intValue();
		int newStart=startPosition;
		int newEnd=startPosition+resolution;
		return new SingleInterval(chr, newStart, newEnd);
	}
	
	
	static Map<String, IntervalTree<File>> makeFileTree(File[] files) {
		Map<String, IntervalTree<File>> rtrn=new TreeMap<String, IntervalTree<File>>();
		for(int i=0; i<files.length; i++){
			String name=files[i].getName();
			SingleInterval region=parse(name);
			IntervalTree<File> tree=new IntervalTree<File>();
			if(rtrn.containsKey(region.getReferenceName())){
				tree=rtrn.get(region.getReferenceName());
			}
			tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), files[i]);
			rtrn.put(region.getReferenceName(), tree);
		}
		return rtrn;
	}

	private static SingleInterval parse(String name) {
		String chr=name.split("_")[0];
		int start=new Integer(name.split("_")[1]);
		int end=new Integer(name.split("_")[2]);
		return new SingleInterval(chr, start, end);
	}
	
	private static Collection<Cluster> parseKmers(String file) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		String nextLine;
		int counter=0;
		while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
			String[] tokens=nextLine.split("\t");
			Cluster c=new Cluster(tokens[0]);
			for(int i=1; i<tokens.length-1; i++){
				SingleInterval region=parseName(tokens[i]);
				c.addRead(region);
			}
			rtrn.add(c);
			counter++;
			if(counter %10000000 ==0 ){System.err.println(counter);}
		}
		
		reader.close();
		return rtrn;
	}
	
	private static SingleInterval parseName(String line) {
		String chr=line.split(":")[0];
		int start=new Integer(line.split(":")[1].split("-")[0]);
		int end=new Integer(line.split(":")[1].split("-")[1]);
		return new SingleInterval(chr, start, end);
	}
	
	public static void main(String[] args) throws IOException{
		/*if(args.length>4){
		
		
		//Collection<Cluster> clusters=parseKmers(args[0]);	
			
		
		
		
		int resolution=1000000;
		if(args.length>6){
			resolution=new Integer(args[6]);
		}
		
		new PermuteKmers(barcodeFileTree, clusters, , save, maxClusterSize, numPerm, resolution);*/
		
		if(args.length>4){
			BarcodingDataStreaming clusters=new BarcodingDataStreaming(new File(args[0]));
			File[] files=new File(args[1]).listFiles();
			Map<String, IntervalTree<File>> barcodeFileTree=makeFileTree(files);
			String sizes=args[2];
			String save=args[3];
			int maxClusterSize=new Integer(args[4]);
			Map<String, Integer> chrSizes=CoordinateSpace.getRefSeqLengthsFromTable(sizes);
			
			int numPerm=100;
			if(args.length>5){
				numPerm=new Integer(args[5]);
			}
			
			boolean writeIGV=new Boolean(args[6]);
			
			int resolution=new Integer(args[7]);
			
		
			PermuteKmers permute=new PermuteKmers(barcodeFileTree, chrSizes, maxClusterSize);
			
			permute.computePermutationsForSpecificRegion(clusters, numPerm, save, resolution);
			
			if(writeIGV){
				permute.writeClusters(clusters, save);
			}
			
			System.out.println("Done");
		}
		else{System.err.println(usage);}
	}
	

	static String usage=" args[0]=clusters (cluster file) \n args[1]=barcoding files \n args[2]=sizes \n args[3]=save \n arg[4]=max cluster size \n args[5]=num perm (default=100) \n args[6]=write .igv files (default=false) \n args[7]=resolution";


	public void writeBEDStack(String saveDir, Collection<Cluster> rawSpriteClusters, int resolution) throws IOException {
		int i=0;
		for(Cluster c: rawSpriteClusters){
			FileWriter writer=new FileWriter(saveDir+"/f"+i+".bed");
			
			for(SingleInterval interval: c.getAllIntervals()){
				writer.write(interval.toBED()+"\n");
			}
			
			writer.close();
			i++;
		}
		
	}
	
	
}

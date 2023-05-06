package guttmanlab.core.barcoding.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;

public class Cluster implements Comparable<Cluster>{
		
	private Collection<SingleInterval> intervals;
	private String barcode;
	private String readString;
	private double score;
	private boolean hasScore=false;
		
	public Cluster(String barcode){
			intervals=new TreeSet<SingleInterval>();
			this.barcode=barcode;
		}
	
	public Cluster(String barcode, Collection<SingleInterval> intervals){
		this.intervals=intervals;
		this.barcode=barcode;
	}
		
	public Cluster(Cluster cluster) {
		this.barcode=cluster.barcode;
		this.intervals=new TreeSet<SingleInterval>();
		
		for(SingleInterval interval: cluster.intervals){
			this.intervals.add(interval);
		}
	}

		public int getSize() {
			return intervals.size();
		}

		public void addRead(SingleInterval interval){
			intervals.add(interval);
		}
		
		public void addReads(Collection<SingleInterval> reads){
			for(SingleInterval read: reads){addRead(read);}
		}
		
		public Collection<SingleInterval> getAllIntervals(){
			return this.intervals;
		}
		
		public String toString(boolean withBarcode){
			String rtrn="";
			if(withBarcode){rtrn+=(this.barcode+"\t");}
			for(SingleInterval interval: this.getAllIntervals()){
				rtrn+=interval.getReferenceName()+":"+interval.getReferenceStartPosition()+"-"+interval.getReferenceEndPosition()+"\t";
			}
			return rtrn;
		}
		
		public String toString(){
			String rtrn=this.barcode;
			for(SingleInterval interval: this.getAllIntervals()){
				rtrn+="\t"+interval.getReferenceName()+":"+interval.getReferenceStartPosition()+"-"+interval.getReferenceEndPosition();
			}
			return rtrn;
		}
		
		public String toKmerString(){
			String rtrn=this.barcode;
			for(SingleInterval interval: this.getAllIntervals()){
				rtrn+="\t"+interval.getReferenceName()+":"+interval.getReferenceStartPosition()+"-"+interval.getReferenceEndPosition();
			}
			return rtrn;
		}
		
		public String getBarcode(){
			return this.barcode;
		}
		
		//TODO Implement equals (order doesn't matter)
		@Override
		public boolean equals(Object o){
			Cluster other=(Cluster)o;
			
			//if(!this.barcode.equalsIgnoreCase(other.getBarcode())){return false;} //TODO remove
			
			if(this.getSize()!=other.getSize()){return false;}
			
			for(SingleInterval read: this.intervals){
				if(!other.getAllIntervals().contains(read)){return false;}
			}
			
		return true;
		}

		@Override
		public int compareTo(Cluster cluster2) {
			
			if(this.getSize()!=cluster2.getSize()){return getSize()-cluster2.getSize();}
			
			Iterator<SingleInterval> iter1=getAllIntervals().iterator();
			Iterator<SingleInterval> iter2=cluster2.getAllIntervals().iterator();
			
			while(iter1.hasNext()){
				SingleInterval interval1=iter1.next();
				SingleInterval interval2=iter2.next();
				if(!interval1.equals(interval2)){return interval1.compareTo(interval2);}
			}
			
			return this.getBarcode().compareTo(cluster2.getBarcode());
		}

		/*
		 * Make a new Cluster that is a clone of the current one excluding interval
		 */
		public Cluster removeRead(SingleInterval interval) {
			Cluster c=new Cluster(this);
			c.intervals.remove(interval);
			return c;
		}

		public boolean containsInterval(SingleInterval interval) {
			return this.intervals.contains(interval);
		}

		public String toStringNoName() {
			String rtrn="";
			for(SingleInterval interval: this.getAllIntervals()){
				rtrn+=interval.toUCSC()+";";
			}
			return rtrn;
		}
		
		public Cluster getPermutedCluster(Map<String, Integer> chrSizes){
			return getPermutedCluster(chrSizes, 1);
		}
		
		public Cluster getPermutedCluster(Map<String, Integer> chrSizes, int resolution){
			Cluster c=randomCluster(chrSizes, resolution);
			return c;
		}
		
		/**
		 * Randomize the structure of the cluster -- including interchromosomal
		 * @param Distance
		 * @param chrSizes 
		 */
		private Cluster randomCluster(Map<String, Integer> chrSizes, int resolution) {
			Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
			
			Map<String, Cluster> subClustersByChr=getSubclusters();
			
			Map<String, Integer> updatedChrSizes=new TreeMap<String, Integer>(chrSizes);
			updatedChrSizes.remove("chrY");
			
			Cluster random=new Cluster("random");
			
			for(String chr: subClustersByChr.keySet()){
				//for each subcluster get random
				Cluster sub=subClustersByChr.get(chr);
				Cluster randomSub=randomIntraCluster(sub, updatedChrSizes);
				random.addReads(randomSub.getAllIntervals());
				
				for(String chrToRemove: randomSub.getAllChromosomes()){
					updatedChrSizes.remove(chrToRemove);
				}
				
			}
			
			random=random.bin(resolution);
			return random;
		}
		
		
		public Collection<String> getAllChromosomes() {
			Collection<String> rtrn=new TreeSet<String>();
			
			for(SingleInterval interval: getAllIntervals()){
				rtrn.add(interval.getReferenceName());
			}
			
			return rtrn;
		}

		private Map<String, Cluster> getSubclusters() {
			Map<String, Cluster> rtrn=new TreeMap<String, Cluster>();
			
			for(SingleInterval interval: getAllIntervals()){
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
			return c;
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

		public int getClusterSize() {
			return this.getAllIntervals().size();
		}

		public boolean isInterchromosomal() {
			Set<String> chrSet=new TreeSet<String>();
			for(SingleInterval interval:getAllIntervals()){
				chrSet.add(interval.getReferenceName());
			}
			return chrSet.size()>1;
		}

		public String toSPRITEFormat() {
			String rtrn=this.barcode;
			for(SingleInterval region: this.getAllIntervals()){
				rtrn+="\t"+region.getReferenceName()+":"+region.getReferenceStartPosition();
			}
			return rtrn;
		}

		public void setReadString(String nextLine) {
			this.readString=nextLine;
		}

		public String getReadString(){
			return this.readString;
		}

		//TODO This needs to be fixed
		public Cluster bin(int resolution) {
			Cluster rtrn=new Cluster(getBarcode());
			for(SingleInterval interval: getAllIntervals()){
				int startIndex=interval.getReferenceStartPosition()/resolution;
				int newStart=startIndex*resolution;
				int newEnd=newStart+Math.max(interval.getLength(), resolution);
				SingleInterval newInterval=new SingleInterval(interval.getReferenceName(), newStart, newEnd);
				rtrn.addRead(newInterval);	
			}
			return rtrn;
		}

		
		public static void main(String[] args){
			SingleInterval i1=new SingleInterval("chr10", 3000000, 10000000);
			SingleInterval i2=new SingleInterval("chr12", 5000000, 16000000);
			SingleInterval i3=new SingleInterval("chr12", 25000000, 31000000);
			
			Cluster c=new Cluster("Test");
			c.addRead(i1);
			c.addRead(i2);
			c.addRead(i3);
			
			System.err.println(c.toString(true));
			System.err.println(i1.size());
			System.err.println(i2.size());
			System.err.println(i3.size());
			
			Cluster random=c.randomCluster(CoordinateSpace.MM9.getRefSizes(), 1000000);
			
			System.err.println(random.toString(true));
			
			for(SingleInterval r1: random.getAllIntervals()){
				System.err.println(r1.size());
			}
			
		}

		/**
		 * Set a score for the cluster
		 * @param score a numeric score
		 */
		public void setScore(double score) {
			this.hasScore=true;
			this.score=score;
		}
		
		/**
		 * Get the numeric score
		 * @return the score
		 */
		public double getScore(){
			return this.score;
		}
		
		/**
		 * Checks whether a score was set for this cluster
		 * @return true if there is a numeric score
		 */
		public boolean hasScore(){
			return this.hasScore;
		}

		public double getObservedFromCluster() throws IOException {
			String[] tokens=getReadString().split("\t");
			double observed=new Double(tokens[tokens.length-1]);
			return observed;
		}

		public double computeInterchromosomalDistance(Annotation geneWindow, Collection<Cluster> clusters) {
			Collection<Double> rtrn=new ArrayList<Double>();
			for(SingleInterval ah: getAllIntervals()){
				if(!ah.getReferenceName().equalsIgnoreCase(geneWindow.getReferenceName())){
					//double count=BarcodingData.quantify(clusters, ah);
					double count=BarcodingData.quantify(clusters, ah, geneWindow);
					rtrn.add(count);
				}
			}
			return Statistics.mean(rtrn);
		}
		
		public Map<Cluster, Double> computeInterchromosomalDistance(BarcodingDataStreaming data) {
			Collection<Cluster> clusters=data.getClustersOverlappingRegion(getAllIntervals());
			
			Map<Cluster,Double> rtrn=new TreeMap<Cluster, Double>();
			
			int counter=0;
			for(SingleInterval region1: getAllIntervals()){
				Collection<Cluster> sc=BarcodingData.getClustersOverlappingRegion(clusters, region1);
				for(SingleInterval region2: getAllIntervals()){
					if(region1!=region2 && !region1.getReferenceName().equalsIgnoreCase(region2.getReferenceName())){
						double count=BarcodingData.quantify(sc, region1, region2);
						Cluster c=new Cluster("c"+counter);
						c.addRead(region2);
						c.addRead(region1);
						rtrn.put(c, count);
						counter++;
					}
				}
				System.err.println(region1.toUCSC());
			}
			return rtrn;
		}
		
		public Collection<Cluster> countInterchromosomalPairs(BarcodingDataStreaming data) {
			
			
			//Map<Cluster, Double> pairs=new TreeMap<Cluster, Double>();
			
			Collection<Cluster> rtrn=new TreeSet<Cluster>();
			
			for(SingleInterval region1: getAllIntervals()){
				System.err.println(region1);
				Collection<Cluster> clusters=data.getClustersOverlappingRegion(region1);
				for(SingleInterval region2: getAllIntervals()){
					//System.err.println(region1.toUCSC()+" "+region2.toUCSC());
					if(region1!=region2 && !region1.getReferenceName().equalsIgnoreCase(region2.getReferenceName())){
						Cluster c=new Cluster("");
						c.addRead(region2);
						c.addRead(region1);
						
						rtrn.addAll(BarcodingData.overlaps(clusters, region1, region2));
					}
				}
				//System.err.println(region1.toUCSC()+" "+rtrn.size()+" "+clusters.size());
			}
			
			//System.err.println("Fraction "+rtrn.size()+" "+clusters.size());
			return rtrn;
		}

		public boolean overlapsInterval(Annotation region) {
			for(SingleInterval interval: intervals){
				if(region.overlaps(interval)){return true;}
			}
			return false;
		}

		public Collection<Cluster> enumerateKmers(int k) {
			Collection<SingleInterval> intervals=getAllIntervals();
				
			Collection<Cluster> clusters=new TreeSet<Cluster>();
			for(SingleInterval interval: intervals){
				Cluster cluster=new Cluster(getBarcode());
				cluster.addRead(interval);
				clusters.add(cluster);
			}
				
			for(int i=2; i<=k; i++){
				clusters=add(intervals, clusters, i);
			}
				
			return clusters;
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

		public boolean hasAdjacentRegions() {
			Iterator<SingleInterval> iter=getAllIntervals().iterator();
			SingleInterval current=iter.next();
			int counter=0;
			while(iter.hasNext()){
				SingleInterval next=iter.next();
				int distance=next.getReferenceStartPosition()-current.getReferenceEndPosition();
				if(current.getReferenceName().equals(next.getReferenceName()) && distance==0){return true;}
				current=next;
				counter++;
			}
			return false;
		}
		
		public int minDistance() {
			Iterator<SingleInterval> iter=getAllIntervals().iterator();
			SingleInterval current=iter.next();
			int minDistance=Integer.MAX_VALUE;
			int counter=0;
			while(iter.hasNext()){
				SingleInterval next=iter.next();
				if(current.getReferenceName().equals(next.getReferenceName())){
					int distance=next.getReferenceStartPosition()-current.getReferenceEndPosition();
					if(distance<minDistance){minDistance=distance;}
				}
				current=next;
				counter++;
			}
			return minDistance;
		}

		public Map<String, Collection<SingleInterval>> getAllIntervalsByChr() {
			Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
			
			for(SingleInterval interval: getAllIntervals()){
				Collection<SingleInterval> temp=new TreeSet<SingleInterval>();
				if(rtrn.containsKey(interval.getReferenceName())){
					temp=rtrn.get(interval.getReferenceName());
				}
				temp.add(interval);
				rtrn.put(interval.getReferenceName(), temp);
			}
			
			return rtrn;
		}

		public boolean containsAllIntervals(Collection<SingleInterval> allIntervals) {
			for(SingleInterval interval: allIntervals){
				boolean overlaps=overlapsInterval(interval);
				if(!overlaps){return false;}
			}
			return true;
		}

		//Only works for intrachromosome
		public Annotation toAnnotation(String name) {
					
			Collection<SingleInterval> blocks=new TreeSet<SingleInterval>();
			Collection<String> chromosomes=new TreeSet<String>();
			for(SingleInterval region: this.getAllIntervals()){
				blocks.add(region);
			}
					
			Annotation spliced=new BlockedAnnotation(blocks, name);
			return spliced;		
		}
		
		
		//Only works for intrachromosome
		public String toBED(String name) {
			
			Collection<SingleInterval> blocks=new TreeSet<SingleInterval>();
			Collection<String> chromosomes=new TreeSet<String>();
			for(SingleInterval region: this.getAllIntervals()){
				blocks.add(region);
			}
			
			Annotation spliced=new BlockedAnnotation(blocks, name);
			return spliced.toBED();		
		}
		
		

		public String toFileName() {
			String rtrn="";
			for(SingleInterval region: this.getAllIntervals()){
				rtrn+=region.getReferenceName()+"_"+region.getReferenceStartPosition();
			}
			return rtrn;
		}

		public boolean containsChr(String chr) {
			return this.getAllChromosomes().contains(chr);
		}

		
		
	}
	


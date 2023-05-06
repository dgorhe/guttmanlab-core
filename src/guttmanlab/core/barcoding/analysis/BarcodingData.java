package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;

/**
 * 
 * @author mitch
 * A data structure to represent the barcoding data
 */
public class BarcodingData{

	private Collection<SingleInterval> allPositions;
	private Map<String, Collection<SingleInterval>> barcodeToPositions;
	private Map<SingleInterval, Collection<String>> interactionMap; //indexed by chromosome
	private Collection<File> barcodeFiles;
	
	public BarcodingData(File barcodeFile) throws IOException{
		this(BarcodeFileParser.getBarcodeToPositions(barcodeFile), barcodeFile);
		//allPositions=BarcodeFileParser.getAllGenomicPositions(barcodeFile);
		//barcodeToPositions=BarcodeFileParser.getBarcodeToPositions(barcodeFile);
	}
	
	public BarcodingData(Iterator<File> barcodeFiles) throws IOException{
		this(BarcodeFileParser.getBarcodeToPositions(barcodeFiles), barcodeFiles);
	}
	
	public BarcodingData(File barcodeFile, String chr) throws IOException{
		this(BarcodeFileParser.getBarcodeToPositions(barcodeFile, chr), barcodeFile);
	}
	
	public BarcodingData(File barcodeFile, SingleInterval region) throws IOException{
		this(BarcodeFileParser.getBarcodeToPositions(barcodeFile, region), barcodeFile);
	}
	
	public Collection<File> getBarcodeFiles(){return this.barcodeFiles;}
	
	public BarcodingData(Map<String, Collection<SingleInterval>> data, Collection<File> barcodeFiles) {
		this.barcodeToPositions=data;
		this.barcodeFiles=barcodeFiles;
	}
	
	public BarcodingData(Map<String, Collection<SingleInterval>> data, Iterator<File> barcodeFiles) {
		this.barcodeToPositions=data;
		this.barcodeFiles=new ArrayList<File>();
		while(barcodeFiles.hasNext()){this.barcodeFiles.add(barcodeFiles.next());}
	}
	
	public BarcodingData(Map<String, Collection<SingleInterval>> data, File barcodeFile) {
		this.barcodeToPositions=data;
		this.barcodeFiles=new ArrayList<File>();
		this.barcodeFiles.add(barcodeFile);
	}
	
	public BarcodingData(Collection<Cluster> data) {
		this.barcodeToPositions=new TreeMap<String, Collection<SingleInterval>>();
		for(Cluster c: data){
			this.barcodeToPositions.put(c.getBarcode(), c.getAllIntervals());
		}
		
	}
	
	public BarcodingData() {
		this.barcodeToPositions=new TreeMap<String, Collection<SingleInterval>>();
		this.barcodeFiles=new ArrayList<File>();
	}

	public Collection<String> getBarcodes() {
		return barcodeToPositions.keySet();
	}

	public Collection<SingleInterval> getPositionsWithBarcode(String barcode) {
		return barcodeToPositions.get(barcode);
	}
	
	public Cluster getClusterWithBarcode(String barcode){
		Cluster rtrn=new Cluster(barcode);
		if(this.barcodeToPositions.containsKey(barcode)){
			rtrn.addReads(this.barcodeToPositions.get(barcode));
		}
		return rtrn;
	}
	
	private Cluster bin(Cluster cluster, int resolution) {
		if(resolution<=1){return cluster;}
		Cluster rtrn=new Cluster(cluster.getBarcode());
		for(SingleInterval interval: cluster.getAllIntervals()){
			int startIndex=interval.getReferenceStartPosition()/resolution;
			int newStart=startIndex*resolution;
			int newEnd=(startIndex+1)*resolution;
			SingleInterval newInterval=new SingleInterval(interval.getReferenceName(), newStart, newEnd);
			rtrn.addRead(newInterval);	
		}
		return rtrn;
	}

	public Collection<SingleInterval> getAllGenomicPositions() {
		if(this.allPositions==null || this.allPositions.isEmpty()){
			this.allPositions=new TreeSet<SingleInterval>();
			for(String barcode: this.barcodeToPositions.keySet()){
				Collection<SingleInterval> regions=barcodeToPositions.get(barcode);
				this.allPositions.addAll(regions);
			}
		}
		return allPositions;
	}

	public Collection<String> getBarcodesOverlappingRegion(SingleInterval region) throws NumberFormatException, IOException {
		if(interactionMap==null || interactionMap.isEmpty()){
			this.interactionMap=new TreeMap<SingleInterval, Collection<String>>();
			for(String barcode: this.barcodeToPositions.keySet()){
				Collection<SingleInterval> regions=this.barcodeToPositions.get(barcode);
				for(SingleInterval r: regions){
					Collection<String> list=new TreeSet<String>();
					if(this.interactionMap.containsKey(r)){list=interactionMap.get(r);}
					list.add(barcode);
					this.interactionMap.put(r, list);
				}
			}
		}
		
		return this.interactionMap.get(region);
	}
	
	
	public List<SingleInterval> getPositionsOnChromosome(String chr) {
		List<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		for(SingleInterval pos: allPositions){
			if(pos.getReferenceName().equalsIgnoreCase(chr)){
				rtrn.add(pos);
			}
		}
		return rtrn;
	}
	
	
	public int quantify(Annotation region1, Annotation region2){
		Collection<Cluster> clusters=getClustersOverlappingRegion(region1);
		clusters=overlaps(clusters, region2);
		return clusters.size();
	}
	
	public static int quantify(Collection<Cluster> clusters, Annotation region2){
		return overlaps(clusters, region2).size();
	}
	
	public static int quantify(Collection<Cluster> clusters, Annotation region1, Annotation region2){
		Collection<Cluster> clusters1=overlaps(clusters, region1);
		Collection<Cluster> clusters2=overlaps(clusters1, region2);
		return clusters2.size();
	}
	
	public static Collection<Cluster> overlaps(Collection<Cluster> clusters, Annotation region1, Annotation region2){
		Collection<Cluster> clusters1=overlaps(clusters, region1);
		Collection<Cluster> clusters2=overlaps(clusters1, region2);
		return clusters2;
	}
	
	public int quantify(Cluster cluster) {
		SingleInterval current=cluster.getAllIntervals().iterator().next();
		Collection<Cluster> clusters=getClustersOverlappingRegion(current);
		
		for(SingleInterval region: cluster.getAllIntervals()){
			clusters=overlaps(clusters, region);
		}
		
		return clusters.size();
		
	}
	
	private static Collection<Cluster> overlaps(Collection<Cluster> clusters, Annotation region) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(Cluster cluster: clusters){
			if(overlaps(cluster, region)){rtrn.add(cluster);}
		}
		return rtrn;
	}
	
	private static boolean overlaps(Cluster cluster, Annotation region) {
		for(SingleInterval interval: cluster.getAllIntervals()){
			if(interval.overlaps(region)){
				return true;
			}
		}
		return false;
	}
	
	public InteractionScore getInteractionScore(SingleInterval region1, SingleInterval region2) throws NumberFormatException, IOException{
		if(region1.equals(region2)){return InteractionScore.self;}
		
		Set<String> sharedBarcodes=new HashSet<String>();
		
		Collection<String> barcodes1=getBarcodesOverlappingRegion(region1);
		Collection<String> barcodes2=getBarcodesOverlappingRegion(region2);
		
		for(String barcode1: barcodes1){
			if(barcodes2.contains(barcode1)){sharedBarcodes.add(barcode1);}
		}
		
		InteractionScore score=new InteractionScore(sharedBarcodes, barcodes1, barcodes2);
		return score;
	}
	
	
	
	
	
	
	private Collection<Annotation> getRegions(String chr1, int resolution, CoordinateSpace genome) {
		int chr1Length=genome.getRefSeqLengths().get(chr1);
		Collection<Annotation> rtrn=new ArrayList<Annotation>();
		
		for(int i=0; i<chr1Length; i+=resolution){
			SingleInterval interval=new SingleInterval(chr1, i, i+=resolution);
			rtrn.add(interval);
		}
		return rtrn;
	}
	
	private Collection<Annotation> getRegions(int resolution, CoordinateSpace genome) {
		Collection<Annotation> rtrn=new ArrayList<Annotation>();
		
		for(String chr: genome.getRefSeqLengths().keySet()){
			rtrn.addAll(getRegions(chr, resolution, genome));
		}
		return rtrn;
	}
	
	

	
	
	public BarcodingData bin(int resolution){
		return bin(resolution, 1, Integer.MAX_VALUE);
	}
	
	public BarcodingData bin(int resolution, int minClusterSize, int maxClusterSize) {
		Map<String, Collection<SingleInterval>> data=new TreeMap<String, Collection<SingleInterval>> ();
		
		for(String barcode: getBarcodes()){
			Collection<SingleInterval> list=getPositionsWithBarcode(barcode);
			List<SingleInterval> unique=bin(list, resolution);
			if(unique.size()>minClusterSize && unique.size()<maxClusterSize){
				data.put(barcode, unique);
			}
			/*else if(unique.size()>=maxClusterSize){
				System.err.println("Excluding "+barcode+" "+unique.size());
			}*/
		}
		
		
		BarcodingData rtrn=new BarcodingData(data, barcodeFiles);
		return rtrn;
	}

	private List<SingleInterval> bin(Collection<SingleInterval> list, int resolution) {
		Set<SingleInterval> set=new HashSet<SingleInterval>();
		for(SingleInterval interval: list){
			int startIndex=interval.getReferenceStartPosition()/resolution;
			int newStart=startIndex*resolution;
			int newEnd=(startIndex+1)*resolution;
			SingleInterval newInterval=new SingleInterval(interval.getReferenceName(), newStart, newEnd);
			set.add(newInterval);
		}
		return new ArrayList<SingleInterval>(set);
	}

	/**
	 * Return clusters that are of at least size n
	 * @param min The size of minimum clusters
	 * @param max 
	 * @return barcodes for clusters bigger than n
	 */
	public Collection<String> getClusters(int min, int max) {
		Collection<String> rtrn=new ArrayList<String>();
		Collection<String> barcodes=this.getBarcodes();
		for(String barcode: barcodes){
			Collection<SingleInterval> intervals=this.getPositionsWithBarcode(barcode);
			int count=intervals.size();
			if(count>min && count<max){rtrn.add(barcode);}
		}
		return rtrn;
	}

	public Collection<Cluster> getClusters() {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(String barcode: this.barcodeToPositions.keySet()){
			Collection<SingleInterval> list=this.barcodeToPositions.get(barcode);
			Cluster c=new Cluster(barcode);
			c.addReads(list);
			rtrn.add(c);
		}
		return rtrn;
	}

	
	public Collection<Cluster> getUniqueClusters() {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(String barcode: this.barcodeToPositions.keySet()){
			Collection<SingleInterval> list=this.barcodeToPositions.get(barcode);
			Cluster c=new Cluster("c");
			c.addReads(list);
			rtrn.add(c);
		}
		return rtrn;
	}
	
	/**
	 * Get reads normalized to resolution
	 * @param cluster
	 * @param resolution
	 * @return
	 */
	public Collection<SingleInterval> getPositionsWithBarcode(String cluster, int resolution) {
		Collection<SingleInterval> reads=getPositionsWithBarcode(cluster);
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval read: reads){
			int newStart=read.getReferenceStartPosition()/resolution;
			int newEnd=newStart+1;
			SingleInterval newInterval=new SingleInterval(read.getReferenceName(), newStart, newEnd);
			rtrn.add(newInterval);
		}
		
		return rtrn;
	}

	
	public Map<Cluster, Set<String>> enumerateAllKmers(int k, int maxClusterSize){
		Map<Cluster, Set<String>> rtrn=new TreeMap<Cluster, Set<String>>(); //Cluster and the barcodes that have it
		for(String barcode: this.barcodeToPositions.keySet()){
			int clusterSize=barcodeToPositions.get(barcode).size();
			if(clusterSize>=k && clusterSize<maxClusterSize){
				Collection<Cluster> clusters=enumerateAllKmers(barcode, k);
				for(Cluster cluster: clusters){
					Set<String> barcodeSet=new HashSet<String>();
					int clusterCount=0;
					if(rtrn.containsKey(cluster)){
						barcodeSet=rtrn.get(cluster);
					}
					barcodeSet.add(barcode);
					rtrn.put(cluster, barcodeSet);
				}
			}
		}
		return rtrn;
	}
	
	/**
	 * Enumerate all k-mers
	 * @param k 
	 */
	public Collection<Cluster> enumerateAllKmers(String barcode, int k){
		Collection<SingleInterval> intervals=barcodeToPositions.get(barcode);
		
		Collection<Cluster> clusters=new TreeSet<Cluster>();
		for(SingleInterval interval: intervals){
			Cluster cluster=new Cluster(barcode);
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
	
	public Collection<Cluster> getClustersOverlappingRegion(Annotation region){
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		Map<String, Collection<SingleInterval>> barcodes=this.barcodeToPositions;
		for(String barcode: barcodes.keySet()){
			Collection<SingleInterval> intervals=barcodes.get(barcode);
			if(overlaps1(intervals, region)){
				Cluster newCluster=new Cluster(barcode);
				newCluster.addReads(intervals);
				rtrn.add(newCluster);
			}
		}
		return rtrn;
	}
	
	
	public static Collection<Cluster> getClustersOverlappingRegion(Collection<Cluster> clusters, Annotation region) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Cluster c: clusters){
			if(overlaps(c, region)){
				rtrn.add(c);
			}
		}
		return rtrn;
	}
	
	public static Collection<Cluster> getClustersOverlappingRegion(Collection<Cluster> clusters, Collection<Annotation> regions) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Annotation region: regions){
			rtrn.addAll(getClustersOverlappingRegion(clusters, region));
		}
		return rtrn;
	}
	
	public Map<String, Cluster> getClustersOverlappingRegionByBarcodeName(Annotation region){
		Map<String, Cluster> rtrn=new TreeMap<String, Cluster>();
		Map<String, Collection<SingleInterval>> barcodes=this.barcodeToPositions;
		for(String barcode: barcodes.keySet()){
			Collection<SingleInterval> intervals=barcodes.get(barcode);
			if(overlaps1(intervals, region)){
				Cluster newCluster=new Cluster(barcode);
				newCluster.addReads(intervals);
				rtrn.put(barcode, newCluster);
			}
		}
		return rtrn;
	}

	private boolean overlaps1(Collection<SingleInterval> intervals, Annotation region) {
		for(SingleInterval interval: intervals){
			if(interval.overlaps(region)){return true;}
		}
		return false;
	}

	
	
	private void write(String save, Collection<Cluster> clusters) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Cluster c: clusters){
			writer.write(c.toString()+"\n");
		}
		
		writer.close();
		
	}


	public Collection<Cluster> getClustersWithBarcodes(Collection<String> names) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(String barcode: names){
			rtrn.add(this.getClusterWithBarcode(barcode));
		}
		return rtrn;
	}

	public void write(String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: this.barcodeToPositions.keySet()){
			Cluster c=this.getClusterWithBarcode(barcode);
			writer.write(c.toString()+"\n");
		}
		
		writer.close();
	}
	
	
	
	/*public static void main(String[] args) throws IOException{
		if(args.length>2){
			File dataFile=new File(args[0]);
			Collection<Gene> regions=BEDFileIO.loadRegionsFromFile(args[1]);
			String save=args[2];
			BarcodingData data=new BarcodingData(dataFile);
			data.writeTable(regions, regions, save);
		}
		else{System.err.println(usage);}
	}*/
	static String usage=" args[0]=barcoding data file \n args[1]=regions \n args[2]=save";

	public void computePairwiseInteractions(String save, String chr) throws NumberFormatException, IOException {
		computePairwiseInteractions(save, new TreeSet<SingleInterval>(), chr);
	}
	
	
	
	public void computePairwiseInteractions(String save, Collection<SingleInterval> excluded) throws NumberFormatException, IOException {
		FileWriter writer=new FileWriter(save);
		//writer.write("Name");
		for(SingleInterval interval2:this.getAllGenomicPositions()){
			writer.write(interval2.toUCSC()+"\t");
		}
		writer.write("\n");
		
		int counter=0;
		for(SingleInterval interval1: this.getAllGenomicPositions()){
			writer.write(interval1.toUCSC());
			boolean overlap1=overlaps(interval1, excluded);
			
			Collection<Cluster> clusters=getClustersOverlappingRegion(interval1);
			System.err.println(interval1.toUCSC()+" "+clusters.size());
			
			for(SingleInterval interval2: this.getAllGenomicPositions()){
				boolean overlap2=overlaps(interval2, excluded);
				double score;
				if(overlap1 || overlap2){score=-1.0;}
				else{
					Collection<Cluster> clusters2=getClustersOverlappingRegion(clusters, interval2);
					score=clusters2.size();
				}
						
				writer.write(score+"\t");
				counter++;
				if(counter % 1000000 ==0){System.err.println(counter+" "+interval1.toUCSC());}
			}
			writer.write("\n");
			
		}
		writer.close();
	}
	
	public void computePairwiseInteractions(String save, Collection<SingleInterval> excluded, String chr) throws NumberFormatException, IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		
		//writer.write("Name");
		for(SingleInterval interval2:this.getAllGenomicPositions()){
			if(interval2.getReferenceName().equalsIgnoreCase(chr)){regions.add(interval2);}
		}
		
		System.err.println(regions.size());
		
		for(SingleInterval region: regions){
			writer.write(region.toUCSC()+"\t");
		}
		
		writer.write("\n");
		
		int counter=0;
		for(SingleInterval interval1: regions){
			writer.write(interval1.toUCSC());
			boolean overlap1=overlaps(interval1, excluded);
			for(SingleInterval interval2: regions){
				boolean overlap2=overlaps(interval2, excluded);
				InteractionScore score;
				if(overlap1 || overlap2){score=InteractionScore.self;}
				else{score=getInteractionScore(interval1, interval2);}
						
				writer.write(score.getScore()+"\t");
				
				
			}
			writer.write("\n");
			counter++;
			if(counter % 1 ==0){System.err.println(counter+" "+interval1.toUCSC());}
		}
		writer.close();
	}

	private boolean overlaps(SingleInterval interval1, Collection<SingleInterval> excluded) {
		for(SingleInterval region: excluded){
			if(region.overlaps(interval1)){return true;}
		}
		return false;
	}

	private int distance(SingleInterval interval1, SingleInterval interval2) {
		return Math.max(interval1.getReferenceStartPosition(), interval2.getReferenceStartPosition())-Math.min(interval1.getReferenceEndPosition(), interval2.getReferenceEndPosition());
	}

	

	

	

	
}

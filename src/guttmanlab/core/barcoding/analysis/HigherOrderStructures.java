package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import guttmanlab.core.annotation.SingleInterval;

/**
 * Identify higher order compartments and their observed frequency
 * @author mitch
 *
 */
public class HigherOrderStructures {

	private Map<Cluster, Collection<String>> clusters;
	private int resolution;
	/**
	 * 
	 * @param data A barcoding dataset
	 * @param min the minimum cluster to consider
	 * @param max the maximum compartment to consider
	 */
	public HigherOrderStructures(BarcodingData data, int min, int max, int resolution){
		this.resolution=resolution;
		//get all subclusters
		this.clusters=enumerateSubclusters(data, min, max);
	}

	/**
	 * 
	 * @param data
	 * @param min
	 * @param max
	 * @param resolution
	 * @return Clusters by barcode
	 */
	private Map<Cluster, Collection<String>> enumerateSubclusters(BarcodingData data, int min, int max) {
		Map<Cluster, Collection<String>> rtrn=new TreeMap<Cluster, Collection<String>>();
		
		//Get clusters
		Collection<String> clusters=data.getClusters(2,10);
		
		//Go through each cluster and enumerate subclusters
		for(String cluster: clusters){
			Collection<SingleInterval> reads=data.getPositionsWithBarcode(cluster, resolution);
			System.err.println("barcode: "+cluster+ "\nNumber of positions "+ reads.size());
			Collection<Cluster> subsetsOfClusters=subsets(reads, cluster, min, max);
			rtrn=merge(subsetsOfClusters, rtrn);
		}
		return rtrn;
		
	}


	private Map<Cluster, Collection<String>> merge(Collection<Cluster> subsetsOfClusters, Map<Cluster, Collection<String>> rtrn) {
		for(Cluster cluster: subsetsOfClusters){
			Collection<String> barcodes=new TreeSet<String>();
			if(rtrn.containsKey(cluster)){
				barcodes=rtrn.get(cluster);
			}
			barcodes.add(cluster.getBarcode());
			rtrn.put(cluster, barcodes);
		}
		return rtrn;
	}

	private Collection<Cluster> subsets(Collection<SingleInterval> reads, String barcode, int min, int max) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		Collection<Cluster> updatedClusters=getPairwise(reads, barcode);
		System.err.println("Pairwise: "+updatedClusters.size());
		for(int i=3; i<=Math.min(max, reads.size()); i++){
			updatedClusters=update(updatedClusters, reads, barcode, i);
			if(i>=min){
				rtrn.addAll(updatedClusters);
			}
			System.err.println(i+" "+updatedClusters.size()+" "+rtrn.size());
		}
		return rtrn;
	}
	
	private Collection<Cluster> update(Collection<Cluster> updatedClusters, Collection<SingleInterval> reads, String barcode, int clusterSize) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(SingleInterval read1: reads){
			for(Cluster cluster2: updatedClusters){
				Cluster c=new Cluster(barcode);
				c.addRead(read1);
				c.addReads(cluster2.getAllIntervals());
				if(c.getSize()==clusterSize){
					rtrn.add(c);
				}
			}
		}
		return rtrn;
	}

	private Collection<Cluster> getPairwise(Collection<SingleInterval> reads, String barcode) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(SingleInterval read1: reads){
			for(SingleInterval read2: reads){
				Cluster c=new Cluster(barcode);
				c.addRead(read1);
				c.addRead(read2);
				if(c.getSize()==2){
					rtrn.add(c);
				}
			}
		}
		return rtrn;
	}
	
	/**
	 * Get all cluster of a specified size
	 * @param size
	 * @return
	 */
	/*public Collection<Cluster> getAllClusters(int size){
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(Cluster c: this.clusters){
			if(c.getSize()==size){
				rtrn.add(c);
			}
		}
		return rtrn;
	}*/

	
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File file=new File(args[0]);
			int min=new Integer(args[1]);
			int max=new Integer(args[2]);
			int resolution=100000;
			BarcodingData data=new BarcodingData(file);
			HigherOrderStructures hoStructures=new HigherOrderStructures(data, min, max, resolution);
			hoStructures.writeClusters("C:/Data/Barcoding/ClusterCount.txt", 3);
		}
		else{
			System.err.println("not enough params");
		}
	}

	private void writeClusters(String save, int minCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(Cluster c: clusters.keySet()){
			if(clusters.get(c).size()>minCounts){
				writer.write(c+"\t"+clusters.get(c).size()+"\t"+clusters.get(c)+"\n");
			}
		}
		writer.close();
	}
	

	
}

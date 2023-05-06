package guttmanlab.core.barcoding.analysis;

import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class ClusterWithScore implements Comparable<ClusterWithScore>{

	Cluster cluster;
	Collection<String> barcodes;
	
	/*public ClusterWithScore(Cluster cluster, Collection<String> barcodes) {
		this.cluster=cluster;
		this.barcodes=barcodes;
	}*/
	
	public ClusterWithScore(Cluster cluster, Collection<Cluster> barcodes) {
		this.cluster=cluster;
		this.barcodes=new TreeSet<String>();
		for(Cluster c: barcodes){
			this.barcodes.add(c.getBarcode());
		}
	}
	
	public String toString() {
		String rtrn="";
		for(SingleInterval interval: cluster.getAllIntervals()){
			rtrn+=interval.toUCSC()+"\t";
		}
		//rtrn+=count+"\t"+cluster.getBarcode();
		rtrn+=getCount();
		return rtrn;
	}

	@Override
	public int compareTo(ClusterWithScore cluster2) {
		int diff=cluster2.getCount()-getCount();
		if(diff!=0){return diff;}
		return cluster.compareTo(cluster2.cluster);
	}

	public int getCount() {
		return this.barcodes.size();
	}
	
	

	public boolean containsMultipleChr() {
		Collection<String> chromosomes=new TreeSet<String>();
		for(SingleInterval interval:cluster.getAllIntervals()){
			chromosomes.add(interval.getReferenceName());
		}
		if(chromosomes.size()>1){return true;}
		return false;
	}

	public boolean isInterchromosomal() {
		return this.cluster.isInterchromosomal();
	}

	public Collection<String> getAllBarcodes() {
		return this.barcodes;
	}
	
}

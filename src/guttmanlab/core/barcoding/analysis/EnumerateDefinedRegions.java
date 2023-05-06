package guttmanlab.core.barcoding.analysis;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class EnumerateDefinedRegions {
	
	Collection<SingleInterval> intervals;

	public EnumerateDefinedRegions(Collection<SingleInterval> intervals, String save) {
		this.intervals=intervals;
		
	}
	

	


	

	private Collection<Cluster> enumerate(Collection<SingleInterval> regions, int k) {
		Collection<Cluster> clusters=new TreeSet<Cluster>();
		for(SingleInterval region1: regions){
			Cluster cluster=new Cluster(region1.getName());
			cluster.addRead(region1);
			clusters.add(cluster);
			for(int i=2; i<=k; i++){
				clusters=add(regions, clusters, i);
			}
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

	private void write(FileWriter writer, Collection<Cluster> clusters) throws IOException {
		for(Cluster c: clusters){
			writer.write(c.toString(true)+"\n");
		}
	}

	public static void main(String[] args) throws IOException{
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[0]);
		String save=args[1];
		new EnumerateDefinedRegions(regions, save);
	}
}

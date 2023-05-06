package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import Jama.Matrix;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.simulation.AssemblySize;

public class HubStrength {
	
	int numPerm=10;
	int resolution=1000000;

	public HubStrength(BarcodingDataStreaming data, Cluster hub, String save) throws IOException{
		
		Collection<Cluster> clusters=hub.countInterchromosomalPairs(data); 
		
		//cluster(clusters);
		
		write(save+".clusters", clusters);
		
		System.err.println("Observed "+clusters.size());
		
		/*for(int i=0;i< numPerm; i++){
			Cluster random=hub.getPermutedCluster(space.getRefSeqLengths());
			int randomVal=random.countInterchromosomalPairs(data);
			System.err.println("Random "+i+" "+randomVal);
		}*/
		
		
		
	}

	private void cluster(Collection<Cluster> clusters) {
		Map<Double, Collection<Cluster>> distanceClusters=new TreeMap<Double, Collection<Cluster>>();
		
		//Start with individual
		for(Cluster cluster1: clusters){
			System.err.println(cluster1.getBarcode());
			for(Cluster cluster2: clusters){
				double distance=distance(cluster1, cluster2);
				Collection<Cluster> list=new TreeSet<Cluster>();
				if(distanceClusters.containsKey(distance)){list=distanceClusters.get(distance);}
				list.add(cluster1);
				list.add(cluster2);
				distanceClusters.put(distance, list);
			}	
		}
		
		//Take minimum
		double minDistance=distanceClusters.keySet().iterator().next();
		Collection<Cluster> minClusters=distanceClusters.get(minDistance);
		
		System.err.println(minDistance+" "+minClusters);
		
		
	}

	private double distance(Cluster c1, Cluster c2) {
		Cluster cluster1=c1.bin(resolution);
		Cluster cluster2=c2.bin(resolution);
		
		Collection<SingleInterval> allRegions=new TreeSet<SingleInterval>();
		allRegions.addAll(cluster1.getAllIntervals());
		allRegions.addAll(cluster2.getAllIntervals());
		
		double distance=0;
		for(SingleInterval region: allRegions){
			double s1=0.0;
			double s2=0.0;
			if(cluster1.containsInterval(region)){s1=1.0;}
			if(cluster2.containsInterval(region)){s2=1.0;}
			distance+=Math.pow(s2-s1,2);
		}
		distance=Math.sqrt(distance);
		return distance;
	}

	private void write(String save, Collection<Cluster> clusters) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Cluster c: clusters){
			writer.write(c.toSPRITEFormat()+"\n");
		}
		
		writer.close();
	}

	private void write(String save, Map<Cluster, Double> map, Map<Cluster, Double>[] randomVals) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Cluster c: map.keySet()){
			writer.write(map.get(c)+"\t");
		}
		writer.write("\n");
		
		for(int i=0; i<randomVals.length; i++){
			for(Cluster c: randomVals[i].keySet()){
				writer.write(randomVals[i].get(c)+"\t");
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			System.err.println("Loaded data");
			
			
			Collection<SingleInterval> activeHubRegions=BEDFileIO.loadSingleIntervalFromFile(args[1]);
			
			Cluster activeHub=new Cluster("Active", activeHubRegions);
			
			System.err.println(activeHub.getAllIntervals().size()+" "+activeHubRegions.size());
			String save=args[2];
			
			
			
			new HubStrength(data, activeHub, save);
		
		}
		else{
			System.err.println(usage);
		}
	}

	static String usage=" args[0]=barcoding data (Files) \n args[1]=speckle regions \n args[2]=save";
	
	
	
}

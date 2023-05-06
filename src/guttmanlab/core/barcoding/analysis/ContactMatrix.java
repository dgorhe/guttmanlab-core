package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.plaf.synth.Region;

import Jama.Matrix;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.simulation.AssemblySize;

public class ContactMatrix {
	
	int numPerm=10;
	int resolution=1000000;
	int maxSize=1000;

	public ContactMatrix(BarcodingDataStreaming data, Collection<SingleInterval> hub,  Map<String, IntervalTree<SingleInterval>> regions, String save) throws IOException{
		
		List<String> rows=new ArrayList<String>();
		for(SingleInterval region: hub) {rows.add(region.toUCSC());}
		
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, rows);
		
		for(String row: rows) {
			for(String column: rows) {
				if(row.split(":")[0].equals(column.split(":")[0])) {mwh.set(row, column, -1);}
			}
		}
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			System.out.println(c.getBarcode()+"\t"+c.getClusterSize());
			if(c.getClusterSize()<maxSize) {
			for(SingleInterval region1: c.getAllIntervals()) {
				String row=getRow(regions, region1);
				if(row!=null) {
					for(SingleInterval region2: c.getAllIntervals()) {
						String column=getRow(regions, region2);
						if(!region1.equals(region2) && column!=null) {
							if(!region1.getReferenceName().equals(region2.getReferenceName())) {
							double increment=1.0/(double)c.getClusterSize();
							mwh.incrementCount(row, column, increment);
							}
							//else {mwh.set(row, column, -1);}
						}
					}
			}
			}
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		
		
		
		mwh.write(save);
	}

	private String getRow(Map<String, IntervalTree<SingleInterval>> regions, SingleInterval region2) {
		String rtrn=null;
		
		if(regions.containsKey(region2.getReferenceName())) {
			IntervalTree<SingleInterval> tree=regions.get(region2.getReferenceName());
			Iterator<SingleInterval> iter=tree.overlappingValueIterator(region2.getReferenceStartPosition(), region2.getReferenceEndPosition());
			if(iter.hasNext()) {return iter.next().toUCSC();}
		}
		
		return rtrn;
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
			
			Map<String, IntervalTree<SingleInterval>> tree=BEDFileIO.loadSingleIntervalTree(args[1]);
			
			//System.err.println(activeHub.getAllIntervals().size()+" "+activeHubRegions.size());
			String save=args[2];
			
			
			
			new ContactMatrix(data, activeHubRegions, tree, save);
		
		}
		else{
			System.err.println(usage);
		}
	}

	static String usage=" args[0]=barcoding data (Files) \n args[1]=speckle regions \n args[2]=save";
	
	
	
}

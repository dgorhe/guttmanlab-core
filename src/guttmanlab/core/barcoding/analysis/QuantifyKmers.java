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
import guttmanlab.core.datastructures.IntervalTree;

public class QuantifyKmers {
	int resolution=1000000;
	int maxClusterSize;
	int minClusterSize;
	BarcodingData data;

	public QuantifyKmers(File[] files, BarcodingDataStreaming clusters, int minClusterSize, int maxClusterSize, String save) throws IOException{
		Map<String, IntervalTree<File>> barcodeFileTree=makeFileTree(files);
		
		FileWriter writer=new FileWriter(save);
		
		this.maxClusterSize=maxClusterSize;
		this.minClusterSize=minClusterSize;
		
		data=null;
		SingleInterval currentRegion=null;
		
		while(clusters.hasNext()){
			Cluster c=clusters.next();
			int count=quantify(data, barcodeFileTree, c, currentRegion).size();
			currentRegion=c.getAllIntervals().iterator().next();
			writer.write(c.toStringNoName()+"\t"+count+"\n");
			System.err.println(c.toStringNoName()+" "+count);
		}
		
		writer.close();
	}
	
	private BarcodingData loadData(BarcodingData currentData, Map<String, IntervalTree<File>> barcodingDataFiles, SingleInterval region, int resolution, SingleInterval currentRegion) throws IOException {
		if(currentRegion!=null && currentRegion.overlaps(region)){return currentData;}
		
		
		if(barcodingDataFiles.containsKey(region.getReferenceName())){
			Iterator<File> barcodingData=barcodingDataFiles.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			BarcodingData data=new BarcodingData(barcodingData);
			data=data.bin(resolution, minClusterSize, maxClusterSize);
			//System.err.println("loaded data for "+region.toUCSC());
			return data;
		}
		
		return new BarcodingData();
	}
	
	private Collection<Cluster> quantify(BarcodingData currentData, Map<String, IntervalTree<File>> barcodingDataFiles, Cluster c, SingleInterval currentRegion) throws IOException {
		SingleInterval current=c.getAllIntervals().iterator().next();
		this.data=loadData(currentData, barcodingDataFiles, current, resolution, currentRegion);
		Collection<Cluster> clusters=data.getClustersOverlappingRegion(current);
		
		for(SingleInterval region: c.getAllIntervals()){
			clusters=overlaps(clusters, region);
		}
		
		return clusters;
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
	
	private static Map<String, IntervalTree<File>> makeFileTree(File[] files) {
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
	
	public static void main(String[] args) throws IOException{
		if(args.length>4){
		new QuantifyKmers(new File(args[0]).listFiles(), new BarcodingDataStreaming(new File(args[1])), new Integer(args[2]), new Integer(args[3]), args[4]);
		}
		else{System.err.println(usage);}
		}
	
	
	static String usage=" args[0]=barcode file directory \n args[1]=clusters \n args[2]=minClusterSize \n args[3]=maxClusterSize \n args[4]=save";
	
}

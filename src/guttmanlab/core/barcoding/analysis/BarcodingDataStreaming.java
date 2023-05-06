package guttmanlab.core.barcoding.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.Statistics;
import htsjdk.samtools.util.CloseableIterator;

/**
 * 
 * @author mitch
 * A data structure to represent the barcoding data
 */
public class BarcodingDataStreaming implements SPRITEData{

	BufferedReader reader;
	String nextLine;
	int resolution;
	File barcodeFile;
	Collection<Predicate<Cluster>> filters; //TODO Use filters to control the iterator
	int minClusterSize=1;
	int maxClusterSize=10000;
	
	public BarcodingDataStreaming(File barcodeFile, int resolution) throws IOException{
		this.barcodeFile=barcodeFile;
		this.reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		this.resolution=resolution;
	}
	
	public BarcodingDataStreaming(File barcodeFile) throws IOException{
		this.barcodeFile=barcodeFile;
		this.reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		this.resolution=1;
	}

	/**
	 * Enumerate all k-mers
	 * @param k 
	 */
	public Collection<Cluster> enumerateAllKmers(Cluster c, int k){
		Collection<SingleInterval> intervals=c.getAllIntervals();
		
		Collection<Cluster> clusters=new TreeSet<Cluster>();
		for(SingleInterval interval: intervals){
			Cluster cluster=new Cluster(c.getBarcode());
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


	@Override
	public boolean hasNext() {
		try {
			this.nextLine=reader.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return nextLine!=null && (nextLine.trim().length()>0);
	}


	@Override
	public Cluster next() {
		Cluster rtrn;
		
		if(nextLine.contains("DNA_") || nextLine.contains("DNA[-]_") || nextLine.contains("DNA[+]_")){
			String[] tokens=nextLine.split("\t");
			
			String name=tokens[0];
			
			rtrn=new Cluster(name);
			rtrn.setReadString(nextLine);
			
			for(int i=1; i<tokens.length; i++) {
				SingleInterval region=new SingleInterval(tokens[i].split("_")[1]);
				rtrn.addRead(region);
			}
			
			//return rtrn;
			
			/*String[] tokens=nextLine.split(";");
			rtrn=new Cluster("");
			rtrn.setReadString(nextLine);
			for(int i=0; i<tokens.length-1; i++){
				if(tokens[i].contains(":")){
					String chr=tokens[i].split(":")[0];
					int start;
					int end;	
					//check if has span
					if(tokens[i].split(":")[1].contains("-")){
						start=new Integer(tokens[i].split(":")[1].split("-")[0]);
						end=new Integer(tokens[i].split(":")[1].split("-")[1]);
					}
					
					else{
						start=new Integer(tokens[i].split(":")[1]);
						end=start+1;
					}
					SingleInterval interval=new SingleInterval(chr, start, end);
					rtrn.addRead(interval);
				}
			}*/
		}
		
		else if(nextLine.contains(";")){
			String[] tokens=nextLine.split(";");
			rtrn=new Cluster("");
			rtrn.setReadString(nextLine);
			for(int i=0; i<tokens.length-1; i++){
				if(tokens[i].contains(":")){
					String chr=tokens[i].split(":")[0];
					int start;
					int end;	
					//check if has span
					if(tokens[i].split(":")[1].contains("-")){
						start=new Integer(tokens[i].split(":")[1].split("-")[0]);
						end=new Integer(tokens[i].split(":")[1].split("-")[1]);
					}
					
					else{
						start=new Integer(tokens[i].split(":")[1]);
						end=start+1;
					}
					SingleInterval interval=new SingleInterval(chr, start, end);
					rtrn.addRead(interval);
				}
			}	
		}
		
		else{
			String[] tokens=nextLine.split("\t");
			String barcode=tokens[0];
			rtrn=new Cluster(barcode);
			rtrn.setReadString(nextLine);
			
			for(int i=1; i<tokens.length; i++){
				if(tokens[i].contains(":")){
					String chr=tokens[i].split(":")[0];
					int start;
					int end;
								
					//check if has span
					if(tokens[i].split(":")[1].contains("-")){
						start=new Integer(tokens[i].split(":")[1].split("-")[0]);
						end=new Integer(tokens[i].split(":")[1].split("-")[1]);
					}
					
					else{
						start=new Integer(tokens[i].split(":")[1]);
						end=start+1;
					}
					
					SingleInterval interval=new SingleInterval(chr, start, end);
					
					rtrn.addRead(interval);
				}
			}
		}
		return bin(rtrn, this.resolution);
	}


	@Override
	public void close() {
		try {
			reader.close();
			this.reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public Map<Cluster, Set<String>> enumerateAllKmers(int k, int maxClusterSize){
		Map<Cluster, Set<String>> rtrn=new TreeMap<Cluster, Set<String>>(); //Cluster and the barcodes that have it
		int counter=0;
		while(this.hasNext()){
			Cluster cluster=this.next();
			int clusterSize=cluster.getSize();
			if(clusterSize>=k && clusterSize<maxClusterSize){
				Collection<Cluster> clusters=enumerateAllKmers(cluster, k);
				for(Cluster cluster2: clusters){
					Set<String> barcodeSet=new HashSet<String>();
					int clusterCount=0;
					if(rtrn.containsKey(cluster2)){
						barcodeSet=rtrn.get(cluster2);
					}
					barcodeSet.add(cluster.getBarcode());
					rtrn.put(cluster2, barcodeSet);
				}
			}
			counter++;
			//if(counter%100000==0){System.err.println(counter);}
		}
		this.close();
		return rtrn;
	}
	
	public Collection<Cluster> getClustersOverlappingRegion(Annotation region, int minClusterSize, int maxClusterSize){
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()>minClusterSize && c.getClusterSize()<maxClusterSize){
				if(overlaps(c, region)){rtrn.add(c);}
			}
			else if(c.getClusterSize()>maxClusterSize){System.err.println("skipped "+c.getBarcode()+" "+c.getClusterSize());}
			counter++;
			if(counter% 1000000 ==0){System.err.println(counter+" retained "+rtrn.size());}
		}
		
		this.close();
		return rtrn;
	}
	
	public Collection<Cluster> getClustersOverlappingRegion(Annotation region){
		return getClustersOverlappingRegion(region, minClusterSize, maxClusterSize);
	}
	
	
	public Collection<Cluster> getClustersOverlappingRegion(Collection<? extends Annotation> regions, int minClusterSize, int maxClusterSize){
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()>minClusterSize && c.getClusterSize()<maxClusterSize){
				if(overlaps(c, regions)){rtrn.add(c);}
			}
			counter++;
			//if(counter% 1000000 ==0){System.err.println(counter+" retained "+rtrn.size());}
		}
		
		this.close();
		return rtrn;
	}
	
	public Collection<Cluster> getClustersOverlappingRegion(Collection<? extends Annotation> regions){
		return getClustersOverlappingRegion(regions, minClusterSize, maxClusterSize);
	}
	
	private boolean overlaps(Cluster cluster, Collection<? extends Annotation> regions) {
		for(Annotation region: regions){
			if(overlaps(cluster, region)){return true;}
		}
		return false;
	}
	
	private boolean overlaps(Cluster cluster, Annotation region) {
		for(SingleInterval interval: cluster.getAllIntervals()){
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

	
	public static void main(String[] args) throws IOException{
		File dataFile=new File(args[0]);
		BarcodingDataStreaming data=new BarcodingDataStreaming(dataFile);
		Annotation region=new SingleInterval(args[1], new Integer(args[2]), new Integer(args[3]));
		String save=args[4];
		Collection<Cluster> clusters=data.getClustersOverlappingRegion(region, 3, 100);
		data.write(save, clusters);
		//data.close();
	}

	public Map<SingleInterval, Double> distance(Cluster hub) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Collection<Cluster> nbClusters=getClustersOverlappingRegion(hub.getAllIntervals());
		System.err.println("NB "+nbClusters.size());
		
		
		for(SingleInterval ah1: hub.getAllIntervals()){
			List<Double> vals=new ArrayList<Double>();
			Collection<Cluster> sc=BarcodingData.getClustersOverlappingRegion(nbClusters, ah1);
			System.err.println("Region "+sc.size());
			for(Annotation ah2: hub.getAllIntervals()){
				if(!ah1.equals(ah2) && !ah1.getReferenceName().equalsIgnoreCase(ah2.getReferenceName())){
					double count=BarcodingData.quantify(sc, ah1, ah2);
					vals.add(count);
				}
			}
			double avg=Statistics.mean(vals);
			System.err.println(ah1.toUCSC()+" "+avg);
			rtrn.put(ah1, avg);
		}
		return rtrn;
	}

	public int getReferenceLength(String rna1) {
		int max=0;
		while(hasNext()) {
			Cluster c=next();
			for(SingleInterval region: c.getAllIntervals()) {
				if(region.getReferenceName().equalsIgnoreCase(rna1)) {
					max=Math.max(max, region.getReferenceEndPosition());
				}
			}
		}
		close();
		return max;
	}
	
}

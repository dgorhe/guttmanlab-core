package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.math.ScanStat;

public class Enumerate2 {

	int minCount=3;
	int resolution=1000000;
	Collection<SingleInterval> allPossible;
	int totalNumClusters;
	String save="/Users/mguttman/Desktop/";
	
	public Enumerate2(BarcodingData data, SingleInterval regionLoaded, int minCount, int resolution) throws IOException{
		this.minCount=minCount;
		this.resolution=resolution;
		
		//Step 0: Store regions to clusters
		Map<SingleInterval, Collection<String>> regionsToCluster=getRegionsToCluster(data);
		
		
		//Step 1: Get all regions contained with regionLoaded based on resolution
		Collection<SingleInterval> regions=getRegions(regionLoaded);
		
		for(SingleInterval region: regions){
			System.err.println(region.toUCSC()+" "+resolution);
			//Step 1: For each region, Initialize by enumerating all pairs
			Map<Cluster, Collection<String>> currentKmers=enumeratePairs(data, region); //Store kmer and its associated list of barcodes
			
			System.err.println("pairs done "+currentKmers.size());
			
			//Step 2: Iterate
			for(int i=3; i<4; i++){
				if(currentKmers.isEmpty()){break;}
				//Step 3: Extend pairs to 3 by adding all possible, for each 
				Map<Cluster, Collection<String>> nextKmers=iterate(currentKmers, data, regionsToCluster, region);
				write(save+i+"mer", nextKmers); //Output all kmers to a file
				currentKmers=nextKmers;
				//TODO Merge k-mers that have exactly the same barcode collections into maximum intersection
			}
		}
	}

	private Collection<SingleInterval> getRegions(SingleInterval regionLoaded) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		Iterator<? extends Annotation> iter=regionLoaded.getGenomicWindows(resolution, resolution).sortedIterator();
		while(iter.hasNext()){
			Annotation a=iter.next();
			rtrn.add(new SingleInterval(a.getReferenceName(), a.getReferenceStartPosition(), a.getReferenceEndPosition()));
		}
		return rtrn;
	}

	//This will merge clusters with same barcodes
	private Map<Cluster, Collection<String>> cleanup(Map<Cluster, Collection<String>> kmers) {
		Map<Cluster, Collection<String>> rtrn=new TreeMap<Cluster, Collection<String>>();
		
		Map<Barcodes, Collection<Cluster>> tmp=new TreeMap<Barcodes, Collection<Cluster>>();
		for(Cluster node: kmers.keySet()){
			Collection<String> barcodeCollection=kmers.get(node);
			Barcodes barcode=new Barcodes(barcodeCollection);
			Collection<Cluster> list=new TreeSet<Cluster>();
			if(tmp.containsKey(barcode)){list=tmp.get(barcode);}
			list.add(node);
			tmp.put(barcode, list);
		}
		
		System.err.println(kmers.size()+" "+tmp.size());
		
		//TODO Collapse tmp
		
		return rtrn;
	}

	private Map<Cluster, Collection<String>> iterate(Map<Cluster, Collection<String>> pairs, BarcodingData data, Map<SingleInterval, Collection<String>> regionsToCluster, SingleInterval referenceRegion) {
		
		Map<Cluster, Collection<String>> triples=new TreeMap<Cluster, Collection<String>>();
		
		double sum=0;
		double total=0;
		int counter=0;
		
		int countTotal=0;
		int countSig=0;
		
		//For each current kmer
		for(Cluster pair: pairs.keySet()){
			//Get possible extensions and the barcodes that support them
			Map<SingleInterval, Collection<String>> possibleExtensions=getPossible(data, pairs.get(pair), regionsToCluster, pair, referenceRegion); //TODO Get regions that are in at least minCount clusters
			
			//For each extension
			for(SingleInterval extension: possibleExtensions.keySet()){
				//Get barcodes supporting it
				Collection<String> barcodes=possibleExtensions.get(extension);
				
				double ratio1=(double)barcodes.size()/(double)pairs.get(pair).size();
				double ratio2=(double)regionsToCluster.get(extension).size()/(double)totalNumClusters;
				/*Collection<String> test=new TreeSet<String>();
				for(SingleInterval region: triple.getAllIntervals()){
					test.addAll(regionsToCluster.get(region)); //TODO Frequency of A, B, and C
				}*/
				
				//If more than the minimum number of barcodes support extension, then extend
				//if(barcodes.size()>this.minCount){
				
				
				//TODO Compute Binomial
				double p=ratio2;
				int n=pairs.get(pair).size();
				double lambda=n*p;
				//double pval=ScanStat.poisson(barcodes.size(), lambda);
				//double pvalCDF=1-ScanStat.poissonCDF(barcodes.size(), lambda);
				double binomialP=ScanStat.getBinomialPValue(regionsToCluster.get(extension).size(), barcodes.size(), 100000000, pairs.get(pair).size()).doubleValue();
				//double scanP=ScanStat.getPValue(controlCount, sampleCount, controlTotal, sampleTotal, winSize, chrSize)
				
				//System.err.println(barcodes.size()+" "+pairs.get(pair).size()+" "+regionsToCluster.get(extension).size()+" "+totalNumClusters+" "+binomialP+" "+ratio2 +" "+ratio1);
				
				
				
				countTotal++;
				if(binomialP<0.01){countSig++;}
				if(barcodes.size()>this.minCount){	
					//if(binomialP<0.0001){
						Cluster triple=new Cluster("triple"+countTotal);
						triple.addReads(pair.getAllIntervals());
						triple.addRead(extension);
						
						//Don't retain a kmer if it has 2 regions within some defined distance of each other
						if(triple.minDistance()>5*resolution){
							triples.put(triple, barcodes); //TODO Don't retain those that don't pass threshold
						}
					//}
				}
				sum+=barcodes.size();
				total++;
			}
			
			counter++;
			
			if(counter%100 ==0){
				System.err.println(counter+" "+pairs.size()+" "+possibleExtensions.size()+" "+triples.size()+" "+countSig+" "+countTotal);}
		}
		
		System.err.println("Average "+sum+" "+total+" "+(sum/total));
		
		
		
		return triples;
	}

	private boolean adjacent(Cluster triple) {
		return triple.hasAdjacentRegions();
		
	}

	private boolean greaterThanThreshold(Collection<String> list1, Collection<String> list2) {
		Collection<String> smaller=getSmaller(list1, list2);
		
		int counter=0;
		
		for(String barcode: smaller){
			if(list1.contains(barcode) && list2.contains(barcode)){
				counter++;
				if(counter>this.minCount){return true;}
			}
		}
		//System.err.println(list1.size()+" "+list2.size()+" "+rtrn.size());
		return false;
	}
	
	private Collection<String> intersection(Collection<String> list1, Collection<String> list2) {
		Collection<String> rtrn=new TreeSet<String>();
		
		Collection<String> smaller=getSmaller(list1, list2);
		
		for(String barcode: smaller){
			if(list1.contains(barcode) && list2.contains(barcode)){rtrn.add(barcode);}
		}
		//System.err.println(list1.size()+" "+list2.size()+" "+rtrn.size());
		return rtrn;
	}

	private Collection<String> getSmaller(Collection<String> list1, Collection<String> list2) {
		if(list1.size()<list2.size()){return list1;}
		return list2;
	}

	private Map<SingleInterval, Collection<String>> getPossible(BarcodingData data, Collection<String> barcodes, Map<SingleInterval, Collection<String>> regionsToBarcodes, Cluster pair, SingleInterval referenceRegion) {
		//TODO We only need regions that are within minimum number of clusters
		
		Map<SingleInterval, Collection<String>> counts=new TreeMap<SingleInterval, Collection<String>>();
		
		for(Cluster c:data.getClustersWithBarcodes(barcodes)){
			c=c.bin(resolution);
			
			for(SingleInterval region: c.getAllIntervals()){
				if(!pair.getAllIntervals().contains(region) && (region.compareTo(referenceRegion)>0)){ //Only retain regions that are after the current one)
					Collection<String> list=new TreeSet<String>();
					if(counts.containsKey(region)){list=counts.get(region);}
					list.add(c.getBarcode());
					counts.put(region, list);
				}
				
			}
			
		}
				
		return counts;
	}

	private Map<SingleInterval, Collection<String>> getRegionsToCluster(BarcodingData data) {
		Map<SingleInterval, Collection<String>> rtrn=new TreeMap<SingleInterval, Collection<String>>();
		Collection<Cluster> clusters=data.getClusters();
		
		totalNumClusters=0;
		for(Cluster c: clusters){
			c=c.bin(resolution);
			for(SingleInterval region: c.getAllIntervals()){
				Collection<String> list=new TreeSet<String>();
				if(rtrn.containsKey(region)){list=rtrn.get(region);}
				list.add(c.getBarcode());
				rtrn.put(region, list);
			}
			totalNumClusters++;
		}
		
		return rtrn;
	}

	private void write(String save, Map<Cluster, Collection<String>> kmers) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Cluster c: kmers.keySet()){
			writer.write(c.toString(true)+"\t"+kmers.get(c).size()+"\n");
		}
		
		writer.close();
	}

	private Map<Cluster, Collection<String>> enumeratePairs(BarcodingData data, SingleInterval region) {
		Map<Cluster, Collection<String>> rtrn=new TreeMap<Cluster, Collection<String>>();
		this.allPossible=new TreeSet<SingleInterval>();
		Collection<Cluster> clusters=data.getClusters();
		
		for(Cluster c: clusters){
			c=c.bin(resolution);
			if(c.getClusterSize()<100){
				Collection<SingleInterval> regions=c.getAllIntervals();
				allPossible.addAll(regions);
				for(SingleInterval region2: regions){
					if(region2.compareTo(region)>0){ //Only retain regions that are after the current one
						Cluster r=new Cluster("c");
						r.addRead(region);
						r.addRead(region2);;
						Collection<String> list=new TreeSet<String>();
						if(rtrn.containsKey(r)){list=rtrn.get(r);}
						list.add(c.getBarcode());
						rtrn.put(r, list);
					}
				}
			}
		}
		return rtrn;
	}

	private void write(FileWriter writer, Collection<Cluster> pairs) throws IOException {
		for(Cluster c: pairs){
			writer.write(c.toKmerString()+"\n");
		}
		
	}

	private void write(Collection<Cluster> kmers) throws IOException {
		FileWriter writer=new FileWriter("kmers.txt");
		
		for(Cluster c: kmers){writer.write(c.toKmerString()+"\n");}
		
		writer.close();
	}

	private Collection<Cluster> enumeratePairs(SingleInterval region, Collection<SingleInterval> regions) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(SingleInterval region2: regions){
			Cluster c=new Cluster("c");
			c.addRead(region);
			c.addRead(region2);
			if(c.getClusterSize()==2){rtrn.add(c);}
		}
		
		return rtrn;
	}

	private Collection<Cluster> extend(Collection<Cluster> pairs, Collection<SingleInterval> regions) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
				
		for(Cluster c: pairs){
			for(SingleInterval region: regions){
				Cluster extended=new Cluster(c);
				extended.addRead(region); //TODO Make sure it isn't already in it
				if(extended.getClusterSize()> c.getClusterSize()){
					rtrn.add(extended);
				}
			}
		}
		
		return rtrn;
	}

	private Collection<Cluster> filter(Collection<Cluster> pairs, BarcodingData data) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(Cluster c: pairs){
			int count=data.quantify(c);
			if(count>minCount){
				rtrn.add(c);
				//System.err.println(c.toKmerString()+" "+count);
			}
		}
		return rtrn;
	}
	
	/*
	public class Node implements Comparable<Node>{
		Cluster kmer;
		Collection<Cluster> parentKmers;
		Collection<SingleInterval> extensionRegions;
		Collection<SingleInterval> edgesToExclude;
		
		public Node(Cluster kmer){
			this.kmer=kmer;
			this.parentKmers=new TreeSet<Cluster>();
			this.extensionRegions=new TreeSet<SingleInterval>();
			this.edgesToExclude=new TreeSet<SingleInterval>();
		}
		
		public void addIncomingExtensionRegion(SingleInterval extension) {
			this.extensionRegions.add(extension);
			
		}

		public void addExcludedExtensions(Collection<SingleInterval> edgesToExclude) {
			this.edgesToExclude.addAll(edgesToExclude);
		}

		public Collection<Cluster> getParents() {
			return this.parentKmers;
		}


		public Collection<SingleInterval> getExtensionEdges() {
			return extensionRegions;
		}

		public Collection<SingleInterval> getRegionsToExclude(){
			return this.edgesToExclude;
		}
		
		public void addIncomingExtensionRegions(Collection<SingleInterval> incomingList) {
			for(SingleInterval region: incomingList){
				extensionRegions.add(region);
			}
			
		}

		public void addParentKmers(Collection<Cluster> parentNodes) {
			for(Cluster c: parentNodes){
				addParentKmer(c);
			}
			
		}

		public void addParentKmer(Cluster parent){
			this.parentKmers.add(parent);
		}
		

		@Override
		public int compareTo(Node o) {
			return kmer.compareTo(o.kmer);
		}
		
		@Override
		public boolean equals(Object o){
			Node other=(Node)o;
			return compareTo(other)==0;
		}
	}
	*/
	
	
	
	class Barcodes implements Comparable<Barcodes>{
		Collection<String> barcodes;
		public Barcodes(){
			this.barcodes=new TreeSet<String>();
		}
		
		
		
		public boolean equals(Barcodes other){
			return compareTo(other)==0;
		}


		public boolean contains(String barcode) {
			return barcodes.contains(barcode);
		}

		public Barcodes(Collection<String> barcodes) {
			this.barcodes=barcodes;
		}

		public Barcodes(String barcode1, String barcode2) {
			this.barcodes=new TreeSet<String>();
			this.barcodes.add(barcode1);
			this.barcodes.add(barcode2);
		}

		public Collection<String> getBarcodes() {
			return barcodes;
		}

		public int size() {return barcodes.size();}
		
		public String toString(){
			/*String rtrn="";
			for(String b: barcodes){
				rtrn+=b+"_";
			}
			rtrn+=barcodes.size();
			return rtrn;*/
			return barcodes.toString();
		}

		public void addBarcode(String barcode){barcodes.add(barcode);}

		@Override
		public int compareTo(Barcodes o) {
			int compare=o.getBarcodes().size()-this.getBarcodes().size();
			if(compare==0){compare=toString().compareTo(o.toString());}
			return compare;
		}
	}
	
	private static SingleInterval parseFromFileName(String fileName) {
		String[] tokens=fileName.split("_");
		String chr=tokens[0];
		int start=new Integer(tokens[1]);
		int end=new Integer(tokens[2]);
		return new SingleInterval(chr, start, end);
	}
	
	
	public static void main(String[] args) throws IOException{
		long start=System.currentTimeMillis();
		BarcodingData data=new BarcodingData(new File(args[0]));
		int minCount=new Integer(args[1]);
		//SingleInterval region=new SingleInterval("chr3", 35000000, 36000000);
		SingleInterval region=parseFromFileName(new File(args[0]).getName());
		int resolution=new Integer(args[2]);
		new Enumerate2(data, region, minCount, resolution);
		long end=System.currentTimeMillis();
		System.err.println("Done in "+(end-start)/1000.0);
	}

	


	
}

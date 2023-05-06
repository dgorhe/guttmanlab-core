package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;
import net.sf.samtools.util.CloseableIterator;

public class EnumerateAllClustersRewriteEfficient {
	private static int maxClusterSize;
	private SingleInterval loadedRegion;
	private int binLoadSize=100000000;
	private int minClusterSize=2;
	Collection<String> chromosomesVisited;
	private int minNumberToExtend;
	private IntervalTree<SingleInterval> regionsVisited;
	
	
	public EnumerateAllClustersRewriteEfficient(Map<String, IntervalTree<File>> barcodingData, Map<String, Integer> chrSizes, String save, int resolution, Annotation regionToUse, boolean testPairwise) throws IOException{
		// Get all barcodes that are within a region
		//enumerate all pairs of barcodes  and take the intersection
		
		
		
		this.minNumberToExtend=5; //TODO Make this a variable
		this.chromosomesVisited=new TreeSet<String>();
		FileWriter writer=new FileWriter(save);
		SingleInterval loadedRegion=null;
		BarcodingData data=null;
		
		for(String chr: chrSizes.keySet()){
			if(regionToUse==null || chr.equalsIgnoreCase(regionToUse.getReferenceName())){
				//BarcodingData data=loadChrData(barcodingData, chr, resolution);
				Collection<SingleInterval> bins=getBins(chr, chrSizes.get(chr), resolution);
				this.regionsVisited=new IntervalTree<SingleInterval>();
				
				for(SingleInterval region: bins){
					if(regionToUse==null || region.overlaps(regionToUse)){
						data=loadData(barcodingData, region, resolution, loadedRegion, data);
						loadedRegion=region;
						
						long start=System.currentTimeMillis();
						//get all clusters overlapping this region
						Map<String, Cluster> clusters=getBarcodes(data, region);
						Map<SingleInterval, Barcodes> uniqueBarcodeSets=makeReadsToBarcode(clusters, region); //Pairwise
						
						//Get Triples
						Collection<Barcodes> triples=getTriples(uniqueBarcodeSets, region, clusters); //TODO As we add combos, we can check if already have superstring - even EXACT match
						
						//Enumerate all pairwise
						Map<Cluster, Collection<String>> list=enumeratePairwise(triples, clusters);
						
						//Map<Cluster, Integer> scores=quantify(list, data);
						
						//Write
						int clusterNum=0;
						for(Cluster cws: list.keySet()){
							int score=list.get(cws).size();
							if(score>1){writer.write(clusterNum+"\t"+cws.toStringNoName()+"\t"+score+"\n");}
							clusterNum++;
						}
						
						writeAsBED(save+".bed", list);
						
						long end=System.currentTimeMillis();
						System.err.println(region.toUCSC()+" reads in region "+clusters.size()+" triple sets "+triples.size()+" number of kmers "+list.size()+" "+" time: "+(end-start)/1000.0);
						regionsVisited.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), region);
					}
				}
				chromosomesVisited.add(chr);
			}
		}
		writer.close();
	}
	
	
	private Map<Cluster, Integer> quantify(Map<Cluster, Collection<String>> list, BarcodingData data) throws NumberFormatException, IOException {
		Map<Cluster, Integer> rtrn=new TreeMap<Cluster, Integer>();
		for(Cluster c: list.keySet()){
			int score=QuantifyFrequencyOfKmers.quantify(c, data);
			rtrn.put(c, score);
		}
		return rtrn;
	}

	private void writeAsBED(String save, Map<Cluster, Collection<String>> list) throws IOException {
		FileWriter writer=new FileWriter(save, true);
		
		for(Cluster cluster: list.keySet()){
			int count=list.get(cluster).size();
			//On the same chromosome
			//Observed at least 5 times?
			if(!cluster.isInterchromosomal() && count>=5){
				Annotation a=new BlockedAnnotation(cluster.getAllIntervals(), "C="+list.get(cluster).size());
				writer.write(a.toBED()+"\n");
			}
		}
		
		writer.close();
	}


	private void writeKmers(String save, Map<Cluster, Collection<String>> list) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int clusterNum=0;
		for(Cluster cws: list.keySet()){
			int score=list.get(cws).size();
			if(score>1){writer.write(clusterNum+"\t"+cws.toStringNoName()+"\t"+score+"\n");}
			clusterNum++;
		}
		
		writer.close();
	}


	private Collection<Barcodes> getTriples(Map<SingleInterval, Barcodes> uniqueBarcodeSets, SingleInterval region, Map<String, Cluster> clusters) {
		Collection<Barcodes> rtrn=new TreeSet<Barcodes>();
		int counter=0;
		for(SingleInterval region2: uniqueBarcodeSets.keySet()){
			Barcodes barcodes=uniqueBarcodeSets.get(region2);
			Collection<SingleInterval> reads=getAllReads(clusters, barcodes, uniqueBarcodeSets);
			//System.err.println(region.toUCSC()+" "+region2.toUCSC()+" Read count: "+reads.size()+" Barcode count: "+barcodes.getBarcodes().size());
			for(SingleInterval read: reads){
				if(!read.equals(region2) && !passed(region, read)){
					Barcodes overlappingBarcodes=uniqueBarcodeSets.get(read);
					Barcodes intersection=intersect(barcodes, overlappingBarcodes);
					//System.err.println(read.toUCSC()+" Barcodes: "+overlappingBarcodes.getBarcodes().size()+" intersection: "+intersection.getBarcodes().size());
					if(intersection.size()>=this.minNumberToExtend){rtrn.add(intersection);} //TODO See if there is an easy way to compare substrings to a set we save
				}
			}
			counter++;
			if(counter%1000==0){System.err.println(counter+" "+uniqueBarcodeSets.size()+" "+rtrn.size());}
		}
		
		return rtrn;
	}


	private boolean passed(SingleInterval currentInterval, SingleInterval read) {
		if(!read.getReferenceName().equalsIgnoreCase(currentInterval.getReferenceName())){
			if(this.chromosomesVisited.contains(read.getReferenceName())){return true;}
			return false;
		}
		if(this.regionsVisited.hasOverlappers(read.getReferenceStartPosition(), read.getReferenceEndPosition())){return true;}
		return false;
	}


	private boolean overlaps(Barcodes intersection, Collection<Barcodes> barcodes) {
		for(Barcodes b: barcodes){
			if(b.isSubset(intersection)){return true;}
		}
		return false;
	}


	private Collection<SingleInterval> getAllReads(Map<String, Cluster> clusters, Barcodes barcodes, Map<SingleInterval, Barcodes> uniqueBarcodeSets) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(String barcode: barcodes.getBarcodes()){
			for(SingleInterval read: clusters.get(barcode).getAllIntervals()){
				if(uniqueBarcodeSets.containsKey(read)){
					rtrn.add(read);
				}
			}
		}
			
		return rtrn;
	}


	private Barcodes intersect(Barcodes barcodes, Barcodes overlappingBarcodes) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String barcode: barcodes.getBarcodes()){
			if(overlappingBarcodes.contains(barcode)){rtrn.add(barcode);}
		}
		
		return new Barcodes(rtrn);
	}


	//Working rewrite
	/*public EnumerateAllClustersRewriteEfficient(File barcodingData, Map<String, Integer> chrSizes, String save, int resolution, String chrToUse, boolean testPairwise) throws IOException{
		System.err.println("Max cluster size "+maxClusterSize);
		int clusterNum=0;
		//Collection<String> visted=new TreeSet<String>();
		visited=new TreeSet<Barcodes>();
		FileWriter writer=new FileWriter(save);
		
		for(String chr: chrSizes.keySet()){
			if(chrToUse.equalsIgnoreCase("all") || chrToUse.equalsIgnoreCase(chr)){
				BarcodingData data=loadChrData(barcodingData, chr, resolution);
				Collection<SingleInterval> bins=getBins(chr, chrSizes.get(chr), resolution);
				
				for(SingleInterval region: bins){
					//get all clusters overlapping this region
					Collection<Cluster> clusters=getBarcodes(data, region);
					Map<SingleInterval, Barcodes> uniqueBarcodeSets=makeReadsToBarcode(clusters, region);
					
					//Enumerate all pairwise
					Map<Cluster, Collection<String>> list=enumeratePairwise(uniqueBarcodeSets, data);
			
					//Collection<Cluster> list=enumeratePairwise(clusters);
					
					int count=0;
					for(Cluster cws: list.keySet()){
						int score=list.get(cws).size();
						if(score>1){writer.write(clusterNum+"\t"+cws.toStringNoName()+"\t"+score+"\n");}
						count++;
						clusterNum++;
					}
					System.err.println(region.toUCSC()+" "+clusters.size());
				}	
			}
		}
		writer.close();
	}*/
	
	private Collection<Cluster> enumeratePairwise(Collection<Cluster> clusters) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
			
			int counter=0;
			
			
				for(Cluster c1: clusters){
					for(Cluster c2: clusters){
						Cluster i=intersect(c1, c2);
						if(i.getClusterSize()>2){rtrn.add(i);}
						counter++;
						if(counter%100000 ==0){System.err.println(counter+" "+clusters.size());}
					}
				}
				
			
			
			return rtrn;
		
	}

	private Map<Cluster, Collection<String>> enumeratePairwise(Collection<Barcodes> uniqueBarcodeSets, Map<String, Cluster> clusters) {
		Map<Cluster, Collection<String>> rtrn=new TreeMap<Cluster, Collection<String>>();
		
		Collection<Barcodes> pairwiseList=new TreeSet<Barcodes>();
		long start=System.currentTimeMillis();
		
		System.err.println("started enumeration");
		
		int counter=0;
		for(Barcodes b: uniqueBarcodeSets){
			Collection<String> barcodes=b.getBarcodes();
			for(String barcode1: barcodes){
				for(String barcode2: barcodes){
					if(!barcode1.equalsIgnoreCase(barcode2)){
						Barcodes pair=new Barcodes(barcode1, barcode2);
						pairwiseList.add(pair);
					}
				}
			}
			counter++;
			if(counter%1000 ==0){System.err.println(counter+" number of sets "+uniqueBarcodeSets.size()+" pairwise combos: "+pairwiseList.size());}
		}
		
		for(Barcodes pair: pairwiseList){
			Iterator<String> iter=pair.getBarcodes().iterator();
			Cluster c1=clusters.get(iter.next());
			Cluster c2=clusters.get(iter.next());
			Cluster i=intersect(c1, c2);
			if(i.getClusterSize()>2){
				Collection<String> list=new TreeSet<String>();
				if(rtrn.containsKey(i)){
					list=rtrn.get(i);
				}
				list.add(c1.getBarcode());
				list.add(c2.getBarcode());
				rtrn.put(i, list);
			}
		}
		
		long end=System.currentTimeMillis();
		System.err.println("Finished making kmers: "+(end-start)/1000.0+" pairwise: "+pairwiseList.size()+" "+(end-start)/1000.0);
		
		return rtrn;
	}

	private boolean checkOverlap(Barcodes b, Collection<Barcodes> combinationsVisited) {
		return overlaps(b, combinationsVisited);
	}


	private Cluster intersect(Cluster c1, Cluster c2) {
		Collection<SingleInterval> temp=new TreeSet<SingleInterval>();
		for(SingleInterval region: c1.getAllIntervals()){
			if(c2.containsInterval(region)){temp.add(region);}
		}
		Cluster rtrn=new Cluster("merged");
		rtrn.addReads(temp);
		return rtrn;
	}

	private void write(String save, Map<SingleInterval, Barcodes> uniqueBarcodeSets, BarcodingData data) throws IOException {
		FileWriter writer=new FileWriter (save);
		
		for(SingleInterval region: uniqueBarcodeSets.keySet()){
			Collection<String> barcodes=uniqueBarcodeSets.get(region).getBarcodes();
			//Collection<Cluster> clusters=data.getClustersWithBarcodes(barcodes);
			//ClusterWithScore intersect=intersection(clusters, uniqueBarcodeSets.get(region));
			writer.write(region.toUCSC()+"\t"+barcodes.size());
			for(String b: barcodes){
				writer.write("\t"+b);
			}
			//writer.write("\t"+intersect);
			writer.write("\n");
		}
		
		writer.close();
	}

	private void writeIGV(String save, ClusterWithScore cws, BarcodingData data) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("chromosome\tstart\tend\tfeature");
		
		Collection<Cluster> rawClusters=data.getClustersWithBarcodes(cws.getAllBarcodes());
		Collection<SingleInterval> allReads=new TreeSet<SingleInterval>();
		
		for(Cluster c: rawClusters){
			allReads.addAll(c.getAllIntervals());
			writer.write("\t"+c.getBarcode());
		}
		writer.write("\n");
		
		for(SingleInterval read: allReads){
			writer.write(read.getReferenceName()+"\t"+read.getReferenceStartPosition()+"\t"+read.getReferenceEndPosition()+"\t"+read.toUCSC());
			for(Cluster c: rawClusters){
				int score=0;
				if(c.containsInterval(read)){score=1;}
				writer.write("\t"+score);
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	private BarcodingData loadData(Map<String, IntervalTree<File>> barcodingDataFiles, SingleInterval region, int resolution, SingleInterval loadedRegion, BarcodingData loadedData) throws IOException {
		if(loadedRegion!=null && loadedData!=null && loadedRegion.contains(region)){
			System.err.println("already loaded "+loadedRegion.toUCSC()+" which contains "+region.toUCSC());
			return loadedData;
		}
		
		if(barcodingDataFiles.containsKey(region.getReferenceName())){
			Iterator<File> barcodingData=barcodingDataFiles.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			BarcodingData data=new BarcodingData(barcodingData);
			data=data.bin(resolution, minClusterSize, maxClusterSize);
			System.err.println("loaded data for "+region.toUCSC());
			return data;
		}
		
		return new BarcodingData();
	}
	
	private BarcodingData loadChrData(File barcodingData, String currentChr, int resolution) throws IOException {
		BarcodingData data=new BarcodingData(barcodingData, currentChr);
		data=data.bin(resolution, minClusterSize, maxClusterSize);
		System.err.println("loaded data for "+currentChr);
		return data;
	}
	
	private BarcodingData loadChrData(File barcodingData, BarcodingData data, int resolution, SingleInterval region) throws IOException {
		if(loadedRegion==null || !region.overlaps(loadedRegion)){
			loadedRegion=new SingleInterval(region.getReferenceName(), region.getReferenceStartPosition(), region.getReferenceStartPosition()+binLoadSize);
			data=new BarcodingData(barcodingData, loadedRegion);
			data=data.bin(resolution);
			System.err.println("loaded data for "+loadedRegion.toUCSC());
		}
		return data;
	}

	private Map<String, Cluster> getBarcodes(BarcodingData data, SingleInterval region) {
		return data.getClustersOverlappingRegionByBarcodeName(region);
	}

	private Collection<SingleInterval> getBins(String currentChr, Integer size, int resolution) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		for(int i=0; i<size; i+=resolution){
			SingleInterval interval=new SingleInterval(currentChr, i, i+resolution);
			rtrn.add(interval);
		}
		return rtrn;
	}

	private void write(String save, ClusterWithScore cws, BarcodingData data) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Cluster> rawClusters=data.getClustersWithBarcodes(cws.getAllBarcodes());
		for(Cluster c: rawClusters){
			writer.write(c.toString()+"\n");
		}
		
		writer.close();
	}

	


	private Collection<ClusterWithScore> intersect(Collection<Barcodes> uniqueBarcodeSets, BarcodingData data) {
		Collection<ClusterWithScore> rtrn=new TreeSet<ClusterWithScore>();
		for(Barcodes barcodes: uniqueBarcodeSets){
			Collection<Cluster> clusters=data.getClustersWithBarcodes(barcodes.getBarcodes());
			ClusterWithScore intersect=intersection(clusters, barcodes);
			rtrn.add(intersect);
		}
		return rtrn;
	}


	private ClusterWithScore intersection(Collection<Cluster> clusters, Barcodes barcodes) {
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		
		for(SingleInterval region: clusters.iterator().next().getAllIntervals()){
			boolean inAll=true;
			for(Cluster c: clusters){
				if(!c.containsInterval(region)){inAll=false;}
			}
			if(inAll){regions.add(region);}
		}
		
		Cluster c=new Cluster(barcodes.toString());
		c.addReads(regions);
		return new ClusterWithScore(c, clusters);
	}


	private Map<SingleInterval, Barcodes> makeReadsToBarcode(Map<String, Cluster> clusters, SingleInterval pulledRegion) {
		//Go through each cluster and assign reads to shared sets of barcodes
		
		Map<SingleInterval, Barcodes> temp=new TreeMap<SingleInterval, Barcodes>();
		//Step 1: Get reads to barcode names
		for(String barcode: clusters.keySet()){
			Cluster c=clusters.get(barcode);
			for(SingleInterval region: c.getAllIntervals()){
				if(!region.equals(pulledRegion) && !this.passed(pulledRegion, region)){
					Barcodes b=new Barcodes();
					if(temp.containsKey(region)){
						b=temp.get(region);
					}
					b.addBarcode(barcode);
					temp.put(region, b);
				}
				//else{System.err.println("pairwise: "+pulledRegion.toUCSC()+" Skipped "+region.toUCSC());}
			}
		}
		
		return temp;
		
		//Go through barcodes and intersect
		/*Collection<Barcodes> rtrn=new TreeSet<Barcodes>();
		
		for(SingleInterval interval: temp.keySet()){
			rtrn.add(temp.get(interval));
		}
		
		return rtrn;*/
	}


	private Set<String> getBarcodes(PopulatedWindow<SAMFragment> region) {
		//Go through the reads and get the name
		Set<String> rtrn=new TreeSet<String>();
		Iterator<SAMFragment> iter=region.getAnnotationsInWindow();
		
		while(iter.hasNext()){
			SAMFragment read=iter.next();
			String barcodeName=read.getName().split("::")[1];
			String newBarcode=barcodeName.replaceAll("\\[", "").replaceAll("\\]", ".");
			newBarcode=newBarcode.substring(0, newBarcode.length()-1);
			rtrn.add(newBarcode);
		}
		return rtrn;
	}





	public static void main(String[] args) throws IOException{
		if(args.length>4){
			File[] files=new File(args[0]).listFiles();
			File sizesFile=new File(args[1]);
			int resolution=new Integer(args[2]);
			String save=args[3];
			String chr=args[4];
			System.err.println("V5.0");
			
			int maxClusterSize=500;
			if(args.length>5){maxClusterSize=new Integer(args[5]);}
			
			EnumerateAllClustersRewriteEfficient.maxClusterSize=maxClusterSize;
			
			Map<String, Integer> sizes=CoordinateSpace.getRefSeqLengthsFromTable(sizesFile.getAbsolutePath());
			
			Annotation region=null;
			if(chr.contains(":")){
				String refName=chr.split(":")[0];
				int start=new Integer(chr.split(":")[1].split("-")[0]);
				int end=new Integer(chr.split(":")[1].split("-")[1]);
				region=new SingleInterval(refName, start, end);
			}
			else if(chr.contains("chr")){
				region=new SingleInterval(chr, 0, sizes.get(chr));
			}
			
			Map<String, IntervalTree<File>> barcodeFileTree=makeFileTree(files);
			
			new EnumerateAllClustersRewriteEfficient(barcodeFileTree, sizes, save, resolution, region, true);
			System.out.println("Done");
		}
		else{System.err.println(usage);}
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

	static String usage=" args[0]=preprocessed barcoding file directory \n args[1]=sizes file \n args[2]=resolution \n args[3]=save \n args[4]=chr (or ALL, or chrA:X-Y) \n args[5]=max cluster size (optional, default=500)";
	
}

class Barcodes implements Comparable<Barcodes>{
	Collection<String> barcodes;
	public Barcodes(){
		this.barcodes=new TreeSet<String>();
	}
	
	/**
	 * Check if one barcode is a complete subset of the other
	 * @param other Second barcode to compare to this
	 * @return true if one is  acomplete subset of the other
	 */
	public boolean isSubset(Barcodes other) {
		if(this.getBarcodes().size()==other.getBarcodes().size()){
			return this.equals(other);
		}
		//go through barcodes and check if one is a complete subset
		
		//iterate over barcodes in small one and make sure they are all in bigger one
		Barcodes smaller=getSmaller(this, other);
		Barcodes larger=getLarger(this, other);
	
		for(String barcode: smaller.getBarcodes()){
			if(!larger.contains(barcode)){return false;}
		}
		return true;
	}
	
	public boolean equals(Barcodes other){
		return compareTo(other)==0;
	}

	private Barcodes getLarger(Barcodes b1, Barcodes b2) {
		if(b1.getBarcodes().size()>b2.getBarcodes().size()){return b1;}
		return b2;
	}
	
	private Barcodes getSmaller(Barcodes b1, Barcodes b2) {
		if(b1.getBarcodes().size()<b2.getBarcodes().size()){return b1;}
		return b2;
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

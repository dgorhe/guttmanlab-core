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

public class SplitClusters {
	int resolution=1000000;
	Map<String, Collection<Cluster>> barcodeToPairs;
	Map<Cluster, Collection<String>> currentKmers;
	
	public SplitClusters(BarcodingData data) throws IOException{
	
		barcodeToPairs=new TreeMap<String, Collection<Cluster>>();
		currentKmers=enumeratePairs(data);
		System.err.println("Done enumerating "+currentKmers.size());
		
		for(Cluster c: data.getClusters()){
			//Make a graph where nodes are region and edges are pairwise distances
			//trim all edges < cutoff
			// enumerate all connected graphs
			
			if(barcodeToPairs.containsKey(c.getBarcode())){
				FileWriter writer=new FileWriter("/Users/mguttman/Desktop/temp/"+c.getBarcode()+".dot");
				writer.write("graph {\n");
				
				//for each kmer, get pairwise scores
				Collection<Cluster> pairs=barcodeToPairs.get(c.getBarcode());
				for(Cluster pair: pairs){
					int score=currentKmers.get(pair).size();
					Iterator<SingleInterval> regions=pair.getAllIntervals().iterator();
					SingleInterval node1=regions.next();
					SingleInterval node2=regions.next();
					writer.write(node1.toNode()+" -- "+node2.toNode()+" "+"[label="+score+",weight="+score+"]\n");
					double average=this.average(pairs);
					//System.err.println(average+" "+stdev+" "+min+" "+max);
					if(currentKmers.get(pair).size()<(average)){
						System.err.println(c.getBarcode());
					}
				}
				writer.write("}");
				writer.close();
			}
		}
		
	}
	
	private Map<Cluster, Collection<String>> enumeratePairs(BarcodingData data) {
		Map<Cluster, Collection<String>> rtrn=new TreeMap<Cluster, Collection<String>>();
		Collection<Cluster> clusters=data.getClusters();
		
		int count=0;
		for(Cluster c: clusters){
			if(c.getClusterSize()<100){
				Collection<Cluster> pairs=new TreeSet<Cluster>();
				c=c.bin(resolution);
				Collection<SingleInterval> regions=c.getAllIntervals();
				
				for(SingleInterval region: regions){
					for(SingleInterval region2: regions){
						if(region.compareTo(region2)>0){
							Cluster pair=new Cluster("pair");
							pair.addRead(region);
							pair.addRead(region2);
							Collection<String> list=new TreeSet<String>();
							if(rtrn.containsKey(pair)){list=rtrn.get(pair);}
							list.add(c.getBarcode());
							rtrn.put(pair, list);
							pairs.add(pair);
						}
					}
				}
				this.barcodeToPairs.put(c.getBarcode(), pairs);
				count++;
				if(count%1000 ==0){
					System.err.println(count+" "+clusters.size() + " "+c.getClusterSize()+" "+pairs.size());
				}
			}
		}
		return rtrn;
	}
	
	private int min(Collection<Cluster> pairs) {
		int min=Integer.MAX_VALUE;
		for(Cluster c: pairs){
			int count=this.currentKmers.get(c).size();
			min=Math.min(min, count);
		}
		return min;
	}
	
	private int max(Collection<Cluster> pairs) {
		int max=0;
		for(Cluster c: pairs){
			int count=this.currentKmers.get(c).size();
			max=Math.max(max, count);
		}
		return max;
	}
	
	private double average(Collection<Cluster> pairs) {
		double sum=0;
		double total=0;
		for(Cluster c: pairs){
			int count=this.currentKmers.get(c).size();
			sum+=count;
			total++;
		}
		return sum/total;
	}
	
	private double stdev(Collection<Cluster> pairs, double avg) {
		double sum=0;
		double total=-1;
		for(Cluster c: pairs){
			int count=this.currentKmers.get(c).size();
			sum+=Math.pow(count-avg,2);
			total++;
		}
		
		
		return Math.sqrt(sum/total);
	}

	public static void main(String[] args) throws IOException{
		BarcodingData data=new BarcodingData(new File(args[0]));
		
		new SplitClusters(data);
	}
	
}

package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.barcoding.analysis.BarcodingDataStreaming;
import guttmanlab.core.barcoding.analysis.Cluster;

public class PullInteractionsWithSpecificRNA {
	int resolution=10000;
	int maxCluster=1000;

	public PullInteractionsWithSpecificRNA(BarcodingDataStreaming data, SingleInterval region, String save) throws IOException{
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		
		Map<SingleInterval, Integer> random=new TreeMap<SingleInterval, Integer>();
		
		System.err.println(region.toUCSC());
		
		Map<Integer, Integer> allClusterSizes=new TreeMap<Integer, Integer>();
		Map<Integer, Integer> retainedClusterSizes=new TreeMap<Integer, Integer>();
		
		Map<SingleInterval, Integer> singles=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		int total=0;
		while(data.hasNext()){
			Cluster c=data.next();
			
			
			//System.err.println(c);
			if(c.getClusterSize()<maxCluster){
				total++;
				int sizeCount=0;
				if(allClusterSizes.containsKey(c.getClusterSize())){sizeCount=allClusterSizes.get(c.getClusterSize());}
				sizeCount++;
				allClusterSizes.put(c.getClusterSize(), sizeCount);
				if(c.overlapsInterval(region)){
					
					if(c.getClusterSize()==1){
						SingleInterval single=c.getAllIntervals().iterator().next();
						int count=0;
						if(singles.containsKey(single)){count=singles.get(single);}
						count++;
						singles.put(single, count); //TODO Write this out
					}
					
					
					
					counter++;
					sizeCount=0;
					if(retainedClusterSizes.containsKey(c.getClusterSize())){sizeCount=retainedClusterSizes.get(c.getClusterSize());}
					sizeCount++;
					retainedClusterSizes.put(c.getClusterSize(), sizeCount);
					//System.err.println(c);
					//write all regions in cluster
					
					int observed=0;
					
					for(SingleInterval interval: c.getAllIntervals()){
						//TODO Count how often the same RNA is in the cluster
						if(interval.getReferenceName().equalsIgnoreCase(region.getReferenceName())){observed++;}
					}
					
					Cluster binned=c.bin(resolution);
					for(SingleInterval interval: binned.getAllIntervals()){
						int count=0;
						if(counts.containsKey(interval)){
							count=counts.get(interval);
						}
						
						
						
						count++;
						counts.put(interval, count);
					}
					if(observed>=2){
						System.out.println(observed+" "+c.getClusterSize()+" "+c.toSPRITEFormat());
					}
				}
			}
		}
		data.close();
		System.err.println("number of clusters with "+region.toUCSC()+": "+counter);
		
		double fraction=(double)counter/(double)total;
		counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()<maxCluster){
				double ran=Math.random();
				if(ran<fraction){
					counter++;
					Cluster binned=c.bin(resolution);
					for(SingleInterval interval: binned.getAllIntervals()){
						int count=0;
						if(random.containsKey(interval)){
							count=random.get(interval);
						}
						count++;
						random.put(interval, count);
					}
				}
			}
		}
		
		data.close();
		
		System.err.println("number of random clusters "+": "+counter);
		
		
		
		
		/*for(Integer clusterSize: allClusterSizes.keySet()){
			System.out.println("ALL\t"+clusterSize+"\t"+allClusterSizes.get(clusterSize));
		}
		
		for(Integer clusterSize: retainedClusterSizes.keySet()){
			System.out.println("Retained\t"+clusterSize+"\t"+retainedClusterSizes.get(clusterSize));
		}*/
		
		
		write(save, counts);
		write(save+".random.bedgraph", random);
		write(save+".singles.bedgraph", singles);
	}

	private void write(String save, Map<SingleInterval, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval interval: counts.keySet()){
			writer.write(interval.getReferenceName()+"\t"+interval.getReferenceStartPosition()+"\t"+interval.getReferenceEndPosition()+"\t"+counts.get(interval)+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		SingleInterval region=new SingleInterval(args[1]);
		String save=args[2];
		new PullInteractionsWithSpecificRNA(data, region, save);
	}
	
}

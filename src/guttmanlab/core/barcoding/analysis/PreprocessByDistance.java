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

public class PreprocessByDistance {

	int freq=500;
	
	public PreprocessByDistance(BarcodingDataStreaming clusters, String save, int frequency, int minNumberOfObservations) throws IOException{
		this.freq=frequency;
		new File(save).mkdir();
				
		Map<Distance, String> distanceToFileWriter=new TreeMap<Distance, String>();
		Map<Cluster, String> cache=new TreeMap<Cluster, String>();
		
		//Go through each and get distances, then write out to individual files
		int counter=0;
		int fileCounter=0;
		int numSkipped=0;
		
		FileWriter currentWriter=null;
		while(clusters.hasNext()){
			Cluster cluster=clusters.next();
			
			double observedCount=cluster.getObservedFromCluster();
			
			if(observedCount>=minNumberOfObservations){
				String currentFile=save+"/f"+(fileCounter/freq);
				if(!cluster.isInterchromosomal()){
					Distance d=getDistance(cluster);
					if(fileCounter%freq ==0){
						if(currentWriter!=null){currentWriter.close();}
						currentWriter=new FileWriter(currentFile);
					}
					
					if(distanceToFileWriter.containsKey(d)){
						System.err.println(distanceToFileWriter.get(d));
						//Add to cache
						cache.put(cluster, distanceToFileWriter.get(d));
					}
					else{
						fileCounter++;
						distanceToFileWriter.put(d, currentFile);
						currentWriter.write(cluster.getReadString()+"\n");
					}
					
					counter++;
					if(counter%100 ==0){System.err.println("record iteration ... "+counter+" "+fileCounter+" "+d.toString());}
									
				}
				
				if(cache.size() ==5000000){
					System.err.println("clearing cache "+counter);
					clearCache(cache);
					cache=new TreeMap<Cluster, String>();
				}
			}
			else{numSkipped++;}
			
		}
		if(currentWriter!=null){currentWriter.close();}
		clearCache(cache);
		
	}
	
	private void clearCache(Map<Cluster, String> cache) throws IOException {
		Map<String, Collection<Cluster>> temp=new TreeMap<String, Collection<Cluster>>();
		
		for(Cluster c: cache.keySet()){
			String name=cache.get(c);
			Collection<Cluster> list=new TreeSet<Cluster>();
			if(temp.containsKey(name)){
				list=temp.get(name);
			}
			list.add(c);
			temp.put(name, list);
		}
		
		for(String fileName: temp.keySet()){
			FileWriter writer=new FileWriter(fileName, true);
			for(Cluster c: temp.get(fileName)){
				writer.write(c.getReadString()+"\n");
			}
			writer.close();
		}
		
	}

	/*public PreprocessByDistance(BarcodingDataStreaming clusters, String save) throws IOException{
		new File(save).mkdir();
		Map<Distance, Collection<Cluster>> distanceToClusters=new TreeMap<Distance, Collection<Cluster>>();
		
		Map<Distance, Integer> distanceCounts=getNumber(clusters);
		System.err.println(distanceCounts.size());
		
		
		Map<Distance, FileWriter> distanceToFileWriter=new TreeMap<Distance, FileWriter>();
		
		//Go through each and get distances, then write out to individual files
		int counter=0;
		int fileCounter=0;
		FileWriter writer=null;
		while(clusters.hasNext()){
			Cluster cluster=clusters.next();
			if(!cluster.isInterchromosomal()){
				Distance d=getDistance(cluster);
				Collection<Cluster> clusterList;
				if(distanceToClusters.containsKey(d)){
					clusterList=distanceToClusters.get(d);
				}
				else{
					clusterList=new TreeSet<Cluster>();
				}
				clusterList.add(cluster);
				if(clusterList.size()== distanceCounts.get(d)){
					if(fileCounter%1000 ==0){
						if(writer!=null){writer.close();}
						writer=new FileWriter(save+"/f"+(fileCounter/1000));
					}
					//Write the clusters
					write(writer, d, clusterList);
					//Remove from map
					distanceToClusters.remove(d);
					distanceCounts.remove(d);
					fileCounter++;
					System.err.println(fileCounter+" removed "+d+" "+clusterList.size());
				}
				else{
					distanceToClusters.put(d, clusterList);
				}
				
			}
			counter++;
			if(counter%1000000 ==0){
				System.err.println("record iteration ... "+counter+" "+distanceToClusters.size());
				//write(save, distanceToClusters);
			}
		}
		if(writer!=null){writer.close();}
		clusters.close();
		
		//write(save, distanceToClusters);
		
		
		//closeAll(distanceToFile);
	}*/
	
	private void write(FileWriter writer, Distance d, Collection<Cluster> clusterList) throws IOException {
		for(Cluster cluster: clusterList){
			writer.write(cluster.toKmerString()+"\n");
		}
	}

	private Map<Distance, Integer> getNumber(BarcodingDataStreaming clusters) {
		Map<Distance, Integer> distanceCounts=new TreeMap<Distance, Integer>();
		int counter=0;
		while(clusters.hasNext()){
			Cluster cluster=clusters.next();
			if(!cluster.isInterchromosomal()){
				Distance d=getDistance(cluster);
				int num=0;
				if(distanceCounts.containsKey(d)){
					num=distanceCounts.get(d);
				}
				num++;
				distanceCounts.put(d, num);
				
				
				/*Collection<Cluster> clusterList;
				if(distanceToClusters.containsKey(d)){
					clusterList=distanceToClusters.get(d);
				}
				else{
					clusterList=new TreeSet<Cluster>();
				}
				clusterList.add(cluster);
				distanceToClusters.put(d, clusterList);*/

			}
			counter++;
			if(counter%1000000 ==0){
				System.err.println(counter+" "+distanceCounts.size());
				//write(save, distanceToClusters);
			}
		}
		clusters.close();
		return distanceCounts;
	}

	private void write(String saveDir, Map<Distance, Collection<Cluster>> distanceToClusters) throws IOException {
		new File(saveDir).mkdir();
		
		int counter=0;
		for(Distance d: distanceToClusters.keySet()){
			FileWriter writer=new FileWriter(saveDir+"/"+d.toFileName(), true);
			
			for(Cluster cluster: distanceToClusters.get(d)){
				writer.write(cluster.toKmerString()+"\n");
			}
			writer.close();
			counter++;
			if(counter %1000 ==0){System.err.println(counter+" "+distanceToClusters.size()+" "+distanceToClusters.get(d).size());}
		}
		
	}

	private void closeAll(Map<Distance, FileWriter> distanceToFile) throws IOException {
		for(Distance d: distanceToFile.keySet()){
			distanceToFile.get(d).close();
		}
	}

	private Distance getDistance(Cluster cluster) {
		int[] distanceList=new int[cluster.getAllIntervals().size()-1];
		Iterator<SingleInterval> iter=cluster.getAllIntervals().iterator();
		SingleInterval current=iter.next();
		Collection<Integer> resolutions=new TreeSet<Integer>();
		resolutions.add(current.getLength());
		int counter=0;
		while(iter.hasNext()){
			SingleInterval next=iter.next();
			resolutions.add(next.getLength());
			int distance=next.getReferenceStartPosition()-current.getReferenceEndPosition();
			distanceList[counter]=distance;
			current=next;
			counter++;
		}
		if(resolutions.size()!=1){throw new IllegalStateException ("Kmers are not of equal resolution");}
		int resolution=resolutions.iterator().next();
		Distance d=new Distance(distanceList, resolution);
		return d;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String saveDir=args[1];
			int frequency=new Integer(args[2]);
			int minNumObservations=new Integer(args[3]);
			System.err.println("V3");
			new PreprocessByDistance(data, saveDir, frequency, minNumObservations);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=kmers \n args[1]=savDir \n args[2]=number of unique distances per file \n args[3]=minimum number of observations to consider";
	
}

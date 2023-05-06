package guttmanlab.core.barcoding.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class MergeTriples {

	public MergeTriples(File[] files, String save) throws IOException{
		//merge all triples into a single Set
		Map<Cluster, Integer> clusters=merge(files);
		
		//write triple frequency
		write(save, clusters);
	}
	
	
	private void write(String save, Map<Cluster, Integer> clusters) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Cluster cluster: clusters.keySet()){
			int count=clusters.get(cluster);
			writer.write(cluster+"\t"+count+"\n");
		}
		
		writer.close();
	}


	private Map<Cluster, Integer> merge(File[] files) throws IOException {
		Map<Cluster, Integer> rtrn=new TreeMap<Cluster, Integer>();
		for(int i=0; i<files.length; i++){
			System.err.println(files[i]);
			Map<Cluster, Integer> clusters=parse(files[i]);
			for(Cluster cluster: clusters.keySet()){
				int count=clusters.get(cluster);
				if(rtrn.containsKey(cluster)){count+=rtrn.get(cluster);}
				rtrn.put(cluster, count);
			}
		}
		return rtrn;
	}


	private Map<Cluster, Integer> parse(File file) throws IOException {
		Map<Cluster, Integer> rtrn=new TreeMap<Cluster, Integer>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
		while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
			String[] tokens=nextLine.split("\t");
			
			String barcode=tokens[0];
			//List<SingleInterval> list=new ArrayList<SingleInterval>();
			Cluster cluster=new Cluster(barcode);
			
			for(int i=1; i<tokens.length-1; i++){
				String chr=tokens[i].split(":")[0];
				int start=new Integer(tokens[i].split(":")[1].split("-")[0]);
				int end=new Integer(tokens[i].split(":")[1].split("-")[1]);;
				SingleInterval interval=new SingleInterval(chr, start, end);
				
				//list.add(interval);
				cluster.addRead(interval);
				
				rtrn.put(cluster, new Integer(tokens[tokens.length-1]));
				
			}
			counter++;
		}
		System.err.println("Read in "+counter +" clusters");
		
		reader.close();
		return rtrn;
		
		
		
	}


	public static void main(String[] args)throws IOException{
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		new MergeTriples(files, save);
	}
	
}

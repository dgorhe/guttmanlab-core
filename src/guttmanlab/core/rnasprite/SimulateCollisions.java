package guttmanlab.core.rnasprite;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

public class SimulateCollisions {
	
	int min=4;
	int max=10000;
	
	public SimulateCollisions(File input, int binNum, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		
		int counter=0;
		ArrayList<String> list=new ArrayList<String>();
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			if(tokens.length>=min && tokens.length<=max) {
				if(counter%binNum ==0) {
					write(writer, counter, binNum, list);
					list=new ArrayList<String>();
				}
				add(list, tokens);
				
				counter++;
				if(counter%1000000 ==0) {System.err.println(counter +" "+binNum);}
			}
		}
		reader.close();
		writer.close();
	}
	

	private void write(FileWriter writer, int counter, int binNum, ArrayList<String> list) throws IOException {
		if(!list.isEmpty()) {
			int bin=counter/binNum;
			writer.write("merged"+bin);
			for(String pos: list) {
				writer.write("\t"+pos);
			}
			writer.write("\n");
		}
	}


	private void add(ArrayList<String> list, String[] tokens) {
		for(int i=2; i<tokens.length; i++) {list.add(tokens[i]);}
	}


	public void randomizeOrder(BarcodingDataStreaming data, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		//go through line by line and assign to merged set
		//Map<Integer, Collection<Cluster>> mergedSet=new TreeMap<Integer, Collection<Cluster>>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			double assignment=Math.random();
			writer.write(assignment+"\t"+c.getLine()+"\n");
			counter++;
			if(counter%1000000==0) {System.err.println(counter+" "+assignment);}
		}
		
		data.close();
		writer.close();
		
	}

	private int getRandom(int numberOfNewClusters) {
		return new Double(Math.random()*numberOfNewClusters).intValue();
	}

	private void write(Map<Integer, Collection<Cluster>> mergedSet, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		for(Integer assignment: mergedSet.keySet()) {
			Collection<Cluster> set=mergedSet.get(assignment);
			Cluster merged=new Cluster("merged"+assignment);
			for(Cluster s: set) {
				merged.addDNAReads(s.getAllDNAIntervals());
				merged.addRNAReads(s.getAllRNARegions());		
			}
			writer.write(merged.toString()+"\n");
			counter++;
			if(counter%1000==0) {System.err.println(counter+" "+set.size()+" "+merged.getClusterSize());}
		}
		
		
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		//BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		File file=new File(args[0]);
		int num=Integer.parseInt(args[1]);
		String save=args[2];
		new SimulateCollisions(file, num, save);
	}
	
}

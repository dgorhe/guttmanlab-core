package guttmanlab.core.clap.updated;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class MakeTable {
	
	private static void write(Map<String, double[]> rtrn, File[] files, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("names");
		for(int i=0; i<files.length; i++) {writer.write("\t"+files[i].getName());}
		writer.write("\n");
		
		for(String key: rtrn.keySet()) {
			writer.write(key);
			double[] vals=rtrn.get(key);
			for(int i=0; i<vals.length; i++) {
				writer.write("\t"+vals[i]);
			}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	
	private static Collection<String> getKeys(SingleInterval region) {
		Collection<SingleInterval> bins=region.allBins(1);
		Collection<String> rtrn=new TreeSet<String>();
		
		//System.err.println(region.toUCSC()+" "+bins.size());
		for(SingleInterval bin: bins) {rtrn.add(bin.toUCSC());}
		
		return rtrn;
	}

	public static void main (String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		
		Map<String, double[]> rtrn=new TreeMap<String, double[]>();
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getAbsolutePath());
			TreeMap<SingleInterval, Double> map=BEDFileIO.loadbedgraph(files[i]);
			for(SingleInterval region: map.keySet()) {
				Collection<String> keys=getKeys(region);
				for(String key: keys) {
					if(!rtrn.containsKey(key)) {rtrn.put(key, new double[files.length]);}
					double[] vals=rtrn.get(key);
					vals[i]=map.get(region);
					rtrn.put(key, vals);
				}
			}
			
			
		}
		
		write(rtrn, files, save);
	}

	/*public static void main (String[] args) throws IOException{
		File file1=new File(args[0]);
		File file2=new File(args[1]);
		
		counts(file1, file2);
	}


	private static void counts(File file1, File file2) throws IOException {
		TreeMap<SingleInterval, Double> map1=BEDFileIO.loadbedgraph(file1);
		TreeMap<SingleInterval, Double> map2=BEDFileIO.loadbedgraph(file2);
		
		Map<String, double[]> rtrn=new TreeMap<String, double[]>();
		
		for(SingleInterval region: map1.keySet()) {
			Collection<String> keys=getKeys(region);
			for(String key: keys) {
				if(!rtrn.containsKey(key)) {rtrn.put(key, new double[2]);}
				double[] vals=rtrn.get(key);
				vals[0]=map1.get(region);
				rtrn.put(key, vals);
			}
		}
			
		for(SingleInterval region: map2.keySet()) {
			Collection<String> keys=getKeys(region);
			for(String key: keys) {
				if(!rtrn.containsKey(key)) {rtrn.put(key, new double[2]);}
				double[] vals=rtrn.get(key);
				vals[1]=map2.get(region);
				rtrn.put(key, vals);
			}
		}
		
	
		int overlap=0;
		int vals1=0;
		int vals2=0;
		for(String key: rtrn.keySet()) {
			double[] vals=rtrn.get(key);
			if(vals[0]>0 && vals[1]>0) {overlap++;}
			if(vals[0]>0) {vals1++;}
			if(vals[1]>0) {vals2++;}
		}
		
		System.err.println(rtrn.keySet().size()+"\t"+vals1+"\t"+vals2+"\t"+overlap);
		
	}*/
		
}
	
	
	

	
	


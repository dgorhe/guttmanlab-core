package guttmanlab.core.test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Cluster;

public class EnumeratePairs {

	public EnumeratePairs(Collection<Cluster> clusters, String save, int min, int max) throws IOException{
		List<String> names=getAllNames(clusters);
		
		clusters=filterClusters(clusters, min, max);
		
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(names, names);
				
		for(Cluster c: clusters){
			MatrixWithHeaders temp=new MatrixWithHeaders(names, names);
			List<String> genes1=c.getRNANameList();
			for(int i=0; i<genes1.size(); i++){
				for(int j=0; j<genes1.size(); j++){
					String gene1=genes1.get(i);
					String gene2=genes1.get(j);
					if(i!=j){
						//System.err.println(gene1+"\t"+gene2);
						double score=temp.get(gene1, gene2);
						score++;
						temp.set(gene1, gene2, score);
					}
					
				}
			}
			
			add(temp, mwh);
			
		}
		
		maskDiagonal(mwh);
		
		
		mwh.write(save);
		
		
	}
	
	private Collection<Cluster> filterClusters(Collection<Cluster> clusters, int min, int max) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		for(Cluster c: clusters){
			if(c.getRNANameList().size()>min && c.getRNANameList().size()<max){rtrn.add(c);}
		}
		
		return rtrn;
	}

	private void maskDiagonal(MatrixWithHeaders mwh) {
		for(String row: mwh.getRowNames()){
			mwh.set(row, row, 0);
		}
		
	}

	private List<String> getAllNames(Collection<Cluster> clusters) {
		List<String> rtrn=new ArrayList<String>();
		
		Set<String> temp=new TreeSet<String>();
		for(Cluster c: clusters){temp.addAll(c.getRNANameList());}
		
		System.err.println(temp.size());
		
		rtrn.addAll(temp);
		return rtrn;
	}

	private void add(MatrixWithHeaders temp, MatrixWithHeaders mwh) {
		for(String row: mwh.getRowNames()){
			for(String column: mwh.getColumnNames()){
				double sum=temp.get(row, column)+mwh.get(row, column);
				mwh.set(row, column, sum);
			}
		}
		
	}

	private static Collection<Cluster> parseClusters(String file, int threshold) throws IOException{
		Collection<String> lines=BEDFileIO.loadLines(file);
		Collection<Cluster> clusters=new ArrayList<Cluster>();
		for(String line: lines){
			Map<String, Integer> counts=new TreeMap<String, Integer>();
			String[] tokens=line.split("\t");
			Cluster c=new Cluster(tokens[0]);
			for(int i=1; i<tokens.length; i++){
				//c.addName(tokens[i]);
				int count=0;
				if(counts.containsKey(tokens[i])){
					count=counts.get(tokens[i]);
				}
				count++;
				counts.put(tokens[i], count);
			}
			
			for(String name: counts.keySet()){
				int count=counts.get(name);
				if(count>threshold){c.addName(name);}
				//if(name.startsWith("DPM")){c.addName(name);}
			}
			
			if(c.getRNANameList().size()>0){
				System.err.println(c.getBarcode()+" "+c.getRNANameList());
			}
			write(tokens[0], counts);
			
			clusters.add(c);
		}
		
		
		
		return clusters;
	}
	
	
	private static void write(String string, Map<String, Integer> counts) {
		System.out.print(string);
		for(String name: counts.keySet()){
			System.out.print("\t"+name+":"+counts.get(name));
		}
		System.out.println();
		
	}

	public static void main(String[] args) throws IOException{
		String file=args[0];
		String save=args[1];
		int threshold=new Integer(args[2]);
		Collection<Cluster>c=parseClusters(file, threshold);
		new EnumeratePairs(c, save, 0, 10000);
	}
	
}

package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class HigherOrderClusters {

	int random=100000;
	
	
	public HigherOrderClusters(BarcodingDataStreaming data, Kmer rnaNames, Map<String, String> rnaClasses, String save) throws IOException{
		
		Collection<Cluster> clusters=data.getRNAClusters(rnaNames, false);
		System.err.println(clusters.size());
		if(clusters.size()>random){
			clusters=getRandom(clusters, random);
		}
		System.err.println("retained "+clusters.size());
		
		MatrixWithHeaders mwh=makeMatrix(clusters, rnaNames, rnaClasses);
		
		mwh.write(save);
	}
	
	
	private Collection<Cluster> getRandom(Collection<Cluster> clusters, int random2) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		double fraction=(double)random2/(double)clusters.size();
		
		for(Cluster c: clusters){
			double rand=Math.random();
			if(rand<=fraction){rtrn.add(c);}
		}
		
		return rtrn;
	}


	private MatrixWithHeaders makeMatrix(Collection<Cluster> clusters, Kmer rnaNames, Map<String, String> rnaClasses) {
		List<String> columns=new ArrayList<String>();
		for(String rna: rnaNames.getRegions()){
			String rnaClass=rnaClasses.get(rna);
			//System.err.println(rna+" "+rnaClass+" "+rna.replaceAll("\"", ""));
			if(!columns.contains(rnaClass)){columns.add(rnaClass);}
		}
		
		List<String> rows=new ArrayList<String>();
		for(Cluster c: clusters){
			rows.add(c.getBarcode());
		}
				
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		
		for(Cluster c: clusters){
			for(String rna: rnaNames.getRegions()){
				if(c.containsRNA(rna) || c.containsRNA(rna.replaceAll("\"", ""))){
					String rnaClass=rnaClasses.get(rna);
					rtrn.set(c.getBarcode(), rnaClass, 1.0);
				}
			}
		}
		
		
		return rtrn;
	}


	private static Kmer makeKmer(String string) throws IOException {
		Kmer rtrn=new Kmer();
		List<String> list=BEDFileIO.loadLines(string);
		for(String line: list){
			rtrn.addRegion(line.split("\t")[0]);
		}
		return rtrn;
	}

	
	private static Map<String, String> parse(String string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			rtrn.put(line.split("\t")[0], line.split("\t")[1]);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Kmer rnaNames=makeKmer(args[1]);
			Map<String, String> classes=parse(args[1]);
			String save=args[2];
			new HigherOrderClusters(data, rnaNames, classes, save);
		}
		else{
			System.err.println(usage);
		}
	}
	
	
	


	static String usage=" args[0]=data \n args[1]=kmer (file) \n args[2]=save";
	
}

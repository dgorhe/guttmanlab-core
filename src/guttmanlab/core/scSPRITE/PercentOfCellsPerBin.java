package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.barcoding.analysis.Cluster;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;

public class PercentOfCellsPerBin {

	int numPerms=1000;
	Map<String, Integer> chrSizes=CoordinateSpace.MM10.getRefSizes();
	int binResolution=1000000;
	
	public PercentOfCellsPerBin(File[] files, String save, int binResolution) throws IOException{
		MatrixWithHeaders counts=null;
		this.binResolution=binResolution;
		
		for(int i=0; i<files.length; i++){
			System.err.println(files[i].getName());
			MatrixWithHeaders mwh=new MatrixWithHeaders(files[i]);
			counts=updateCounts(counts, mwh);
		}
		
		//compute fraction
		MatrixWithHeaders percent=new MatrixWithHeaders(counts.getRowNames(), counts.getColumnNames());
		for(String row: percent.getRowNames()){
			for(String column: percent.getColumnNames()){
				double count=counts.get(row, column);
				double ratio=count/(double)files.length;
				percent.set(row, column, ratio);
			}
		}
		
		percent.write(save);
	}

	
	

	private MatrixWithHeaders updateCounts(MatrixWithHeaders counts, MatrixWithHeaders mwh) {
		
		if(counts==null){
			counts=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		}
		
		for(String row: mwh.getRowNames()){
			if(!counts.containsRow(row)){counts.addRow(row, row);}
			for(String column: mwh.getColumnNames()){
				if(counts.containsColumn(column)){counts.addColumn(column);}
				double score=mwh.get(row, column);
				if(score>0){counts.incrementCount(row, column);}
			}
		}
		return counts;
	}




	private Map<Cluster, Double> score(File[] files, Collection<Cluster> allClusters) throws IOException{
		Map<Cluster, Double> num=new TreeMap<Cluster, Double>();
		Map<Cluster, Double> denom=new TreeMap<Cluster, Double>();
		
		for(int i=0; i<files.length; i++){
			if(i%10==0){System.err.println(i+" "+files.length);}
			MatrixWithHeaders mwh=new MatrixWithHeaders(files[i]);
			
			for(Cluster c: allClusters){
				if(hasKmer(mwh, c)){increment(num, c);}
				increment(denom, c);
			}
		}
		
		System.err.println(num.size()+" "+denom.size());
		
		Map<Cluster, Double> rtrn=new TreeMap<Cluster, Double>();
		for(Cluster c: num.keySet()){
			double numerator=num.get(c);
			double denominator=denom.get(c);
			double fraction=numerator/denominator;
			rtrn.put(c, fraction);
		}
		
		return rtrn;
	}
	
	private void increment(Map<Cluster, Double> num, Cluster c) {
		double val=0;
		if(num.containsKey(c)){val=num.get(c);}
		val=val+1;
		num.put(c, val);
	}


	private double[] fractionOfCells(File[] files, List<Cluster> random, Cluster observed) throws IOException {
		double[] counter=new double[random.size()+1];
		double[] total=new double[random.size()+1];
		
		for(int i=0; i<files.length; i++){
			if(i%10==0){System.err.println(i+" "+files.length);}
			MatrixWithHeaders mwh=new MatrixWithHeaders(files[i]);
			
			if(hasKmer(mwh, observed)){counter[0]++;}
			total[0]++;
			
			for(int j=0; j<random.size(); j++){
				if(hasKmer(mwh, random.get(j))){counter[j+1]++;}
				total[j+1]++;
			}
		}
		
		double[] rtrn=new double[random.size()+1];
		for(int j=0; j<rtrn.length; j++){
			rtrn[j]=counter[j]/total[j];
		}
		return rtrn;
	}

	private double fractionOfCells(File[] files, Cluster kmer) throws IOException {
		double counter=0;
		double total=0;
		for(int i=0; i<files.length; i++){
			MatrixWithHeaders mwh=new MatrixWithHeaders(files[i]);
			if(hasKmer(mwh, kmer)){counter++;}
			total++;
		}
		
		return counter/total;
	}

	private boolean hasKmer(MatrixWithHeaders mwh, Cluster kmer) {
		Iterator<SingleInterval> iter=kmer.getAllIntervals().iterator();
		SingleInterval region1=iter.next();
		SingleInterval region2=iter.next();
		
		Collection<String> bins1=getBins(region1, this.binResolution);
		Collection<String> bins2=getBins(region2, this.binResolution);
		
		for(String bin1: bins1){
			for(String bin2: bins2){
				if(mwh.containsColumn(bin1) && mwh.containsColumn(bin2)){
					double score=mwh.get(bin1, bin2);
					if(score>0){return true;}
				}
			}
		}
		return false;
	}

	private Collection<String> getBins(SingleInterval region1, int binResolution2) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(int i=region1.getReferenceStartPosition(); i<region1.getReferenceEndPosition(); i+=binResolution2){
			SingleInterval region=new SingleInterval(region1.getReferenceName(), i, i+binResolution2);
			rtrn.add(region.toUCSC());
			//System.err.println(region1.toUCSC()+" "+region.toUCSC());
		}
		return rtrn;
	}

	private Collection<Cluster> generatePairs(Collection<SingleInterval> regions) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(SingleInterval region1: regions){
			for(SingleInterval region2: regions){
				if(!region1.getReferenceName().equals(region2.getReferenceName())){
					Cluster c=new Cluster("pairs");
					c.addRead(region1);
					c.addRead(region2);
					rtrn.add(c);
				}
			}
		}
		return rtrn;
	}

	private List<Cluster> getPerms(Cluster kmer, int numPerm) {
		List<Cluster> rtrn=new ArrayList<Cluster>();
		for(int i=0; i<numPerm; i++){
			Cluster randomCluster=kmer.getPermutedCluster(chrSizes, binResolution);
			rtrn.add(randomCluster);
		}
		return rtrn;
	}
	
	
	private static Map<String, MatrixWithHeaders> parse(File[] listFiles) throws IOException {
		Map<String, MatrixWithHeaders> rtrn=new TreeMap<String, MatrixWithHeaders>();
		
		for(int i=0; i<listFiles.length; i++){ //TODO Fix this
			if(i%10==0){System.err.println(i+" "+listFiles.length);}
			rtrn.put(listFiles[i].getName(), new MatrixWithHeaders(listFiles[i]));
		}
		System.err.println("Done loading matrices");
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		//Map<String, MatrixWithHeaders> matrices=parse(new File(args[0]).listFiles());
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		int binResolution=new Integer(args[2]);
		
		
		new PercentOfCellsPerBin(files, save, binResolution);
		
		
		
	}
	
}

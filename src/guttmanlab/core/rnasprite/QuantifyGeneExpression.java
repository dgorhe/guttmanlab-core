package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;

public class QuantifyGeneExpression {
	int binSize=1000000;
	int minClusterCount=100;
	double alpha=0.01;

	
	
	public QuantifyGeneExpression(BarcodingDataStreaming data, String rna, SingleInterval dna, String save) throws IOException{
		
		Collection<String> genes=getGenesInRegion(data, dna);
		
		/*Collection<String> genes=data.getAllRNAs();
		Collection<Cluster> clusters=data.getClustersOverlappingRNA(rna, 1, 1000);
		System.err.println(clusters.size());
		Collection<SingleInterval> dnaRegions=getInterval(clusters);
		
		for(SingleInterval dna: dnaRegions){
			System.err.println(dna.toUCSC());
		}*/
		
		//for(SingleInterval dna: dnaRegions){
			System.err.println(dna.toUCSC());
			Collection<Cluster> dnaClusters=data.getClustersOverlappingRegion(dna,1,1000);
			Pair<Collection<Cluster>> pairs=splitByRNA(dnaClusters, rna);
			for(String gene: genes){
				Pair<Double> txnRate=quantifyTxn(pairs, gene);
				double ratio=txnRate.getValue1()/txnRate.getValue2();
				double logratio=Math.log(ratio)/Math.log(2);
				double p=this.fisherExact(pairs, gene);
				System.err.println(gene+" "+ratio+" "+logratio+" "+p);
			}
		//}
		
		
	}
	
	private Collection<String> getGenesInRegion(BarcodingDataStreaming data, SingleInterval dna) {
		Collection<String> rtrn=new TreeSet<String>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()<1000){
				for(SingleInterval gene: c.getAllRNARegions()){
					if(gene.overlaps(dna)){rtrn.add(gene.getName());}
				}
			}
		}
		
		data.close();
		
		
		return rtrn;
	}

	private Collection<SingleInterval> getInterval(Collection<Cluster> clusters) {
		Map<SingleInterval, Integer> regions=new TreeMap<SingleInterval, Integer>();
		
		for(Cluster c: clusters){
			Cluster binned=c.bin(this.binSize);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				int count=0;
				if(regions.containsKey(region)){
					count=regions.get(region);
				}
				count++;
				regions.put(region, count);
			}
		}
		
		Collection<SingleInterval> merged=filter(regions, this.minClusterCount);
		merged=merge(merged);
		return merged;
	}

	private Collection<SingleInterval> merge(Collection<SingleInterval> merged) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		SingleInterval current=null;
		for(SingleInterval region: merged){
			if(current==null){current=region;}
			else{
				if(current.getReferenceName().equalsIgnoreCase(region.getReferenceName()) && ((region.getReferenceStartPosition()-current.getReferenceStartPosition())<(2*this.binSize))){
					//merge
					SingleInterval temp=new SingleInterval(region.getReferenceName(), Math.min(region.getReferenceStartPosition(), current.getReferenceStartPosition()), Math.max(region.getReferenceEndPosition(), current.getReferenceEndPosition()));
					current=temp;
				}
				else{
					rtrn.add(current);
					current=null;
				}
			}
		}
		
		rtrn.add(current);
		
		return rtrn;
	}

	private Collection<SingleInterval> filter(Map<SingleInterval, Integer> regions, int minCount) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval region: regions.keySet()){
			int count=regions.get(region);
			if(count>minCount){
				System.err.println(region.toUCSC()+"\t"+count);
				rtrn.add(region);}
		}
		
		return rtrn;
	}

	public QuantifyGeneExpression(BarcodingDataStreaming data, String save) throws IOException{
		//MatrixWithHeaders counts=data.getRNADNAContactMatrix(binSize).getValue2();
		
		//bin all genes into 1Mb DNA intervals
		Map<SingleInterval, Collection<String>> genesByBin=getGenesByBin(data);
		Collection<String> allRNA=getAllRNA(genesByBin);
		
		
		MatrixWithHeaders mwh=initializeMatrix(genesByBin, allRNA);
		
		Map<String, Collection<SingleInterval>> genesByChr=getGenesByChr(genesByBin);
		
		for(String chr: genesByChr.keySet()){
			Collection<Cluster> clustersOnChr=data.getClustersOnChr(chr, 1, 1000);
			System.err.println(chr+" "+clustersOnChr.size());
			
			for(SingleInterval genomicBin: genesByChr.get(chr)){
				Collection<String> genes=genesByBin.get(genomicBin);
				Collection<Cluster> clusters=getClustersOverlappingRegion(clustersOnChr, genomicBin);
				//System.err.println(genomicBin.toUCSC()+ " Total clusters "+clusters.size());
				
				if(clusters.size()>minClusterCount){
					
					for(String rna: allRNA){
						Pair<Collection<Cluster>> pairs=splitByRNA(clusters, rna);
						if(pairs.getValue1().size()>minClusterCount){
							for(String gene: genes){
								if(!gene.equals(rna)){
									double p=fisherExact(pairs, gene);
									Pair<Double> txnRate=quantifyTxn(pairs, gene);
									double ratio=txnRate.getValue1()/txnRate.getValue2();
									double logratio=Math.log(ratio)/Math.log(2);
									mwh.set(gene, rna, logratio);
									System.out.println(rna+"\t"+gene+"\t"+pairs.getValue1().size()+"\t"+pairs.getValue2().size()+"\t"+txnRate.getValue1()+"\t"+txnRate.getValue2()+"\t"+ratio+"\t"+p+"\n");	
								}
							}
						}
					}
				
				}
			}
			
			
		}
		
		mwh.write(save);
	}
	
	private Collection<String> getAllRNA(Map<SingleInterval, Collection<String>> genesByBin) {
		Collection<String> allRNA=new ArrayList<String>();
		
		for(SingleInterval region: genesByBin.keySet()){
			allRNA.addAll(genesByBin.get(region));
		}
		
		return allRNA;
	}

	private Map<String, Collection<SingleInterval>> getGenesByChr(Map<SingleInterval, Collection<String>> genesByBin) {
		Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
		
		for(SingleInterval bin: genesByBin.keySet()){
			String chr=bin.getReferenceName();
			Collection<SingleInterval> list=new TreeSet<SingleInterval>();
			if(rtrn.containsKey(chr)){
				list=rtrn.get(chr);
			}
			list.add(bin);
			rtrn.put(chr, list);
		}
		
		return rtrn;
	}

	private MatrixWithHeaders initializeMatrix(Map<SingleInterval, Collection<String>> genesByBin, Collection<String> allRNA) {
		List<String> columns=new ArrayList<String>();
		columns.addAll(allRNA);
		
		List<String> rows=new ArrayList<String>();
		for(SingleInterval bin: genesByBin.keySet()){
			rows.addAll(genesByBin.get(bin));
		}
		
		MatrixWithHeaders mwh= new MatrixWithHeaders(rows, columns);
	
		for(String row: mwh.getRowNames()){
			for(String column: mwh.getColumnNames()){
				mwh.set(row, column, Double.NaN);
			}
		}
		return mwh;
	}

	private Collection<Cluster> getClustersOverlappingRegion(Collection<Cluster> clustersOnChr, SingleInterval genomicBin) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		for(Cluster c: clustersOnChr){
			if(c.containsOverlappingDNA(genomicBin)){rtrn.add(c);}	
		}	
		return rtrn;
	}

	private Map<SingleInterval, Collection<String>> getGenesByBin(BarcodingDataStreaming data) {
		Map<SingleInterval, Collection<String>> rtrn=new TreeMap<SingleInterval, Collection<String>>();
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			Collection<RNAInterval> genes=c.getAllRNARegions();

			for(SingleInterval gene: genes){
				if(!gene.getName().startsWith("Unass")){
					SingleInterval binned=bin(gene, binSize);
					Collection<String> list=new TreeSet<String>();
					if(rtrn.containsKey(binned)){
						list=rtrn.get(binned);
					}
					list.add(gene.getName());
					rtrn.put(binned, list);
				}
			}
			counter++;
			if(counter%100000==0){System.err.println(counter);}
			
		}
		data.close();
		return rtrn;
	}

	private SingleInterval bin(SingleInterval gene, int resolution) {
		int startIndex=gene.getReferenceStartPosition()/resolution;
		int newStart=startIndex*resolution;
		int newEnd=newStart+Math.max(gene.getLength(), resolution);
		SingleInterval newInterval=new SingleInterval(gene.getReferenceName(), newStart, newEnd);
		return newInterval;
	}

	/*public QuantifyGeneExpression(BarcodingDataStreaming data, Collection<Gene> genes, String rna, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		//for each gene, get all clusters containing the DNA
		//split into +lncRNA vs -lncRNA
		//Quantify # of mRNAs
			
		for(Gene gene: genes){
			System.err.println(gene.getName());
			SingleInterval regionToUse=new SingleInterval(gene.getReferenceName(), gene.getReferenceStartPosition()-binSize, gene.getReferenceEndPosition()+binSize);
			Collection<Cluster> clusters=data.getClustersOverlappingRegion(regionToUse, 1, 1000);
			System.err.println("Total clusters "+clusters.size());
			
			Pair<Collection<Cluster>> pairs=splitByRNA(clusters, rna);
			System.err.println("Clusters with RNA "+pairs.getValue1().size()+"\nClusters w/o RNA "+pairs.getValue2().size());
			
			double p=fisherExact(pairs, gene);
			Pair<Double> txnRate=quantifyTxn(pairs, gene);
			double ratio=txnRate.getValue1()/txnRate.getValue2();
			System.out.println(gene.toUCSC()+"\t"+gene.getName()+"\t"+pairs.getValue1().size()+"\t"+pairs.getValue2().size()+"\t"+txnRate.getValue1()+"\t"+txnRate.getValue2()+"\t"+ratio+"\t"+p+"\n");	
		}
		
		writer.close();
	}*/

	private Pair<Collection<Cluster>> splitByRNA(Collection<Cluster> clusters, String rna) {
		Collection<Cluster> withRNA=new ArrayList<Cluster>();
		Collection<Cluster> noRNA=new ArrayList<Cluster>();
		
		for(Cluster c: clusters){
			if(c.containsRNA(rna)){withRNA.add(c);}
			else{noRNA.add(c);}
		}
		
		//noRNA=getTop(noRNA, withRNA.size());
		
		noRNA=sample(noRNA, withRNA);
		
		Pair<Collection<Cluster>> rtrn=new Pair<Collection<Cluster>>(withRNA, noRNA);
		return rtrn;
	}

	private Collection<Cluster> sample(Collection<Cluster> noRNA, Collection<Cluster> withRNA) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		ArrayList<RNAInterval> allRNAs=new ArrayList<RNAInterval>();
		for(Cluster c: noRNA){
			allRNAs.addAll(c.getAllRNARegions());
		}
		
		for(Cluster c: withRNA){
			int rnaClusterSize=c.getAllRNARegions().size();
			Collection<RNAInterval> rnaRegions=get(allRNAs, rnaClusterSize);
			Cluster newCluster=new Cluster("newCluster");
			newCluster.addRNAReads(rnaRegions);
			rtrn.add(newCluster);
		}
		
		return rtrn;
	}

	private Collection<RNAInterval> get(ArrayList<RNAInterval> allRNAs, int rnaClusterSize) {
		Collection<RNAInterval> rtrn=new ArrayList<RNAInterval>();
		for(int i=0; i<rnaClusterSize; i++){
			int index=Double.valueOf(Math.random()*allRNAs.size()).intValue();
			rtrn.add(allRNAs.get(index));
		}
		
		return rtrn;
	}

	private Collection<Cluster> getTop(Collection<Cluster> noRNA, int size) {
		if(noRNA.size()<=size){return noRNA;}
		Map<Integer, Collection<Cluster>> map=new TreeMap<Integer, Collection<Cluster>>();
		
		for(Cluster c: noRNA){
			Collection<Cluster> list=new ArrayList<Cluster>();
			if(map.containsKey(-c.getClusterSize())){list=map.get(-c.getClusterSize());}
			list.add(c);
			map.put(-c.getClusterSize(), list);
		}
		
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		for(Integer s: map.keySet()){
			rtrn.addAll(map.get(s));
			if(rtrn.size()>=size){return rtrn;}
		}
		
		return rtrn;
	}
	
	
	private double fisherExact(Pair<Collection<Cluster>> pairs, String gene) {
		Collection<Cluster> hasRNA=pairs.getValue1();
		Collection<Cluster> noRNA=pairs.getValue2();
		
		Pair<Integer> v1=getCounts(hasRNA, gene);
		Pair<Integer> v2=getCounts(noRNA, gene);
		
		double p1=Statistics.fisherExact(v1.getValue1(), v2.getValue1(), v1.getValue2(), v2.getValue2());
		double p2=Statistics.fisherExact(v1.getValue1(), v1.getValue2(), v2.getValue1(), v2.getValue2());
		
		double ratio1=(double)v1.getValue1()/(double)(v1.getValue1()+v1.getValue2());
		double ratio2=(double)v2.getValue1()/(double)(v2.getValue1()+v2.getValue2());
		
		System.err.println(gene+" "+v1.getValue1()+" "+v2.getValue1()+" "+v1.getValue2()+" "+v2.getValue2()+" "+p1);
		
		return p2;
	}

	private Pair<Integer> getCounts(Collection<Cluster> clusters, String gene) {
		int hasMRNA=0;
		int noMRNA=0;
		for(Cluster c: clusters){
			if(c.containsRNA(gene)){
				hasMRNA++;
			}
			else{
				noMRNA++;
			}
			
		}
		return new Pair<Integer>(hasMRNA, noMRNA);
	}
	

	private Pair<Double> quantifyTxn(Pair<Collection<Cluster>> pairs, String gene) {
		double val1=quantifyTxn(pairs.getValue1(), gene);
		double val2=quantifyTxn(pairs.getValue2(), gene);
		Pair<Double> rtrn=new Pair<Double>(val1, val2);
		return rtrn;
	}

	private double quantifyTxn(Collection<Cluster> clusters, String gene) {
		int count=0;
		int total=0;
		for(Cluster c: clusters){
			if(c.containsRNA(gene)){
				count++;
				//count+=(100000*(2.0/(double)c.getClusterSize()));
			}
			total++;
			//total+=(100000*(2.0/(double)c.getClusterSize()));
		}
		//System.err.println(count+" "+total);
		return (double)count/(double)total;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>1){
			/*BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String save=args[1];
			new QuantifyGeneExpression(data, save);*/
			
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String rna=args[1];
			SingleInterval region=new SingleInterval(args[2]);
			String save=args[3];
			new QuantifyGeneExpression(data, rna,region, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=cluster files \n args[1]=save";
}

package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class RNAContacts {

	int minCount=10;
	double total;
	private double cutoff=2.0;
	//double minRandom=0.1;
	
	public RNAContacts(BarcodingDataStreaming data, String save, Map<String, String> rnas) throws IOException, InterruptedException{
		//TODO Compute input
		
		//Map<String, Double> input=new TreeMap<String, Double>();
		List<String> rnaList=new ArrayList<String>();
		for(String rna: rnas.keySet()){
			String rnaClass=rnas.get(rna);
			if(!rnaList.contains(rnaClass)){rnaList.add(rnaClass);}
		}
		//rnaList.add("input");
		
		MatrixWithHeaders countMatrix=new MatrixWithHeaders(rnaList, rnaList);
		MatrixWithHeaders countMatrixWeighted=new MatrixWithHeaders(rnaList, rnaList);
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()==2){
			Collection<String> genes=getGenes(c);
			for(String gene1: genes){
				//double inputCount=0;
				//if(input.containsKey(gene1)){inputCount=input.get(gene1);}
				//inputCount++;
				//input.put(gene1, inputCount);
				for(String gene2: genes){
					double score=countMatrix.get(gene1, gene2);
					score++;
					countMatrix.set(gene1, gene2, score);
					
					score=countMatrixWeighted.get(gene1, gene2);
					score+=(2.0/c.getClusterSize());
					countMatrixWeighted.set(gene1, gene2, score);
			}
			}
		}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		
		data.close();
		
		for(String class1: countMatrix.getRowNames()){
			countMatrix.set(class1, class1, 0);
			countMatrixWeighted.set(class1, class1, 0);
		}
		
		/*for(String row: countMatrix.getRowNames()){
			if(!row.equals("input")){
				countMatrix.set(row, "input", input.get(row));
				countMatrix.set("input", row, input.get(row));
			}
		}*/
		
		countMatrix.write(save+".all.counts");
		countMatrixWeighted.write(save+".all.weighted");
		
		
		
		Process p=Runtime.getRuntime().exec("python /groups/guttman/SPRITE2/ice_matrix.py "+save+".all.counts");
		p.waitFor();
		
		p=Runtime.getRuntime().exec("python /groups/guttman/SPRITE2/ice_matrix.py "+save+".all.weighted");
		p.waitFor();
	}
	
	
	

	private Collection<String> getGenes(Cluster c) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(SingleInterval region: c.getAllIntervals()){
			rtrn.add(region.getReferenceName());
		}
		
		return rtrn;
	}


	

	private MatrixWithHeaders normalizeByInput(MatrixWithHeaders countMatrix, Map<String, Integer> countsPerGene) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(countMatrix.getRowNames(), countMatrix.getColumnNames());
		for(String row: countMatrix.getRowNames()){
			for(String column: countMatrix.getColumnNames()){
				double score=countMatrix.get(row, column);
				double normFactor=countsPerGene.get(column);
				double norm=score/normFactor;
				rtrn.set(row, column, norm);
			}
		}
		return rtrn;
	}


	


	private void write(String string, Map<String, SingleInterval> rnaRegions) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(String gene: rnaRegions.keySet()){
			SingleInterval region=rnaRegions.get(gene);
			writer.write(gene+"\t"+region.toUCSC()+"\n");
		}
		
		writer.close();
	}


	private void updateRegions(Collection<SingleInterval> allRNARegions, Map<String, SingleInterval> rnaRegions) {
		for(SingleInterval region: allRNARegions){
			//System.err.println(region.getName());
			SingleInterval other=null;
			if(rnaRegions.containsKey(region.getName())){
				other=rnaRegions.get(region.getName());
			}
			other=merge(other, region);
			rnaRegions.put(region.getName(), other);
		}
		
	}


	private SingleInterval merge(SingleInterval other, SingleInterval region) {
		if(other==null){return region;}
		return update(other, region);
	}
	
	private SingleInterval update(SingleInterval rna1, SingleInterval rna2) {
		String chr=rna1.getReferenceName();
		int start=Math.min(rna1.getReferenceStartPosition(), rna2.getReferenceStartPosition());
		int end=Math.max(rna1.getReferenceEndPosition(), rna2.getReferenceEndPosition());
		SingleInterval rtrn=new SingleInterval(chr, start, end);
		rtrn.setName(rna2.getName());
		return rtrn;
	}


	private MatrixWithHeaders filterMatrix(MatrixWithHeaders norm) {
		List<String> list=new ArrayList<String>();
		//filter rows and columns if max<x-fold
		for(String gene: norm.getRowNames()){
			double[] vals=norm.getRow(gene);
			double max=Statistics.max(vals);
			if(max>cutoff){list.add(gene);}
		}
		
		System.err.println("Remaining "+list.size()+" from "+norm.getRowNames().size());
		MatrixWithHeaders temp=new MatrixWithHeaders(list, list);
		
		for(String gene1: list){
			for(String gene2: list){
				temp.set(gene1, gene2, norm.get(gene1, gene2));
			}
		}
				
		return temp;
	}


	private MatrixWithHeaders floorMatrix(MatrixWithHeaders mwh, MatrixWithHeaders countMatrix) {
		for(String row: mwh.getRowNames()){
			for(String column: mwh.getColumnNames()){
				double count=countMatrix.get(row, column);
				if(count<this.minCount){mwh.set(row, column, 0.0);}
			}
		}
		return mwh;
	}


	


	private MatrixWithHeaders getExpected(double sos, List<String> geneNames, Map<String, Integer> countsPerGene) {
		MatrixWithHeaders expectedMatrix=new MatrixWithHeaders(geneNames, geneNames);
		
		for(String gene1: geneNames){
			double p1=(double)countsPerGene.get(gene1)/total;
			for(String gene2: geneNames){
				double p2=((double)countsPerGene.get(gene2)/total);
				double expected=p1*p2*(sos);
				//double expected=getExpected(p1, p2, clusterSizes);
				expectedMatrix.set(gene1, gene2, expected);
			}
		}
		return expectedMatrix;
	}


	private double getExpected(double p1, double p2, ArrayList<Integer> clusterSizes) {
		double rtrn=0;
		
		for(Integer clusterSize: clusterSizes){
			double expected=(p1*clusterSize)*(p2*clusterSize);
			rtrn+=expected;
		}
		
		return rtrn;
	}


	private MatrixWithHeaders normalize(MatrixWithHeaders mwh, MatrixWithHeaders expectedMatrix) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()){
			for(String column: mwh.getColumnNames()){
				if(row==column){rtrn.set(row, column, Double.NaN);}
				else{
					double num=mwh.get(row, column);
					double denom=expectedMatrix.get(row, column);
					double ratio=num/denom;
					rtrn.set(row, column, ratio);
				}
			}
		}
		
		return rtrn;
	}


	private double expectation(Map<String, Integer> countsPerGene, String gene1, String gene2, int clusterSize, int numPerm) {
		int count=0;
		for(int i=0; i<numPerm; i++){
			boolean hit=random(countsPerGene, gene1, gene2, clusterSize);
			if(hit){count++;}
		}
		return (double)count/((double)numPerm);
	}
	
	private double expectation(Map<String, Integer> countsPerGene, String gene1, String gene2, int clusterSize) {
		double p1=clusterSize*((double)countsPerGene.get(gene1)/total);
		double p2=clusterSize*((double)countsPerGene.get(gene2)/total);
		
		double rtrn=p1*p2;
		return rtrn;
	}
	
	private boolean random(Map<String, Integer> countsPerGene, String gene1, String gene2, int clusterSize) {
		double p1=(double)countsPerGene.get(gene1)/total;
		double p2=(double)countsPerGene.get(gene2)/total;
		
		//roll 2 dices
		double random1=Math.random();
		double random2=Math.random();
		
		if(random1<p1*clusterSize && random2<p2*clusterSize){return true;}
		return false;
	}

	private MatrixWithHeaders normalize(MatrixWithHeaders mwh, Map<String, Integer> countsPerGene) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()){
			for(String column: mwh.getColumnNames()){
				double val=mwh.get(row, column);
				int rowCount=countsPerGene.get(row);
				int columnCount=countsPerGene.get(column);
				double rowRatio=(double)rowCount/((double)mwh.getRowNames().size());
				double columnRatio=(double)columnCount/(double)mwh.getColumnNames().size();
				double normscore=val/(columnRatio*rowRatio);
				rtrn.set(row, column, normscore);
			}
		}
		
		return rtrn;
	}

	/*private Map<String, Integer> getCountsPerGene(BarcodingDataStreaming data, String save, String type) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		this.total=0;
		Map<String, SingleInterval> rnaRegions=new TreeMap<String, SingleInterval>();
		
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(this.useCluster(c, type)){
				updateRegions(c.getAllRNARegions(), rnaRegions);
				
				Collection<String> geneNames=c.getRNANames();
				for(String gene: geneNames){
					int count=0;
					if(rtrn.containsKey(gene)){count=rtrn.get(gene);}
					count++;
					total++;
					rtrn.put(gene, count);
				}
				counter++;
				if(counter%100000==0){System.err.println(counter);}
			}
		}
		
		write(save+".chromosomePositions.annotation", rnaRegions);
		data.close();
		return rtrn;
	}*/

	private List<String> getNames(Map<String, Integer> data) {
		List<String> list=new ArrayList<String>();
		
		for(String gene: data.keySet()){
			if(data.get(gene)>this.minCount){list.add(gene);}
		}
		
		return list;
	}

	

	
	private static Map<String, String> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(String line: lines){
			rtrn.put(line.split("\t")[0].replaceAll("\"", ""), line.split("\t")[1]);
		}
		return rtrn;
	}



	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>1){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String save = args[1];
			Map<String, String> genes;
			if(args.length>2){
				genes=parse(args[2]);
			}
			else{genes=makeGenes(data);}
			new RNAContacts(data, save, genes);
		}
		else{System.err.println(usage);}
	}
	
	
	private static Map<String, String> makeGenes(BarcodingDataStreaming data) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			for(SingleInterval region: c.getAllIntervals()){
				String rna=region.getReferenceName();
				rtrn.put(rna, rna);
			}
		}
		data.close();
		
		return rtrn;
	}


	static String usage=" args[0]=clusters \n args[1]=save \n args[2]=list of genes";
	
}

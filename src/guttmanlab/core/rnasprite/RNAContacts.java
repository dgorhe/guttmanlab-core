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

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;

public class RNAContacts {

	int minCount=100;
	double total;
	private double cutoff=2.0;
	//double minRandom=0.1;
	
	
	public RNAContacts(BarcodingDataStreaming data, String save, Map<String, String> rnas, Collection<String> introns, Collection<String> exons) throws IOException, InterruptedException{
		//TODO Compute input
		
		Kmer kmer=new Kmer();
		List<String> rnaList=new ArrayList<String>();
		for(String rna: rnas.keySet()){
			kmer.addRegion(rna);
			String rnaClass=rnas.get(rna);
			if(!rnaList.contains(rnaClass)){rnaList.add(rnaClass);}
		}
		
		/*for(String rna: introns){
			rnas.put(rna+".intron", "intron");
		}
		
		for(String rna: exons){
			rnas.put(rna+".exon", "exon");
		}*/
		
		/*Collection<String> filtered=computeExonIntronCounts(data, mRNAs);
		
		for(String mRNA: filtered){
			rnaList.add(mRNA+".exon");
			rnaList.add(mRNA+".intron");
		}*/
		
		System.err.println(exons.size()+" "+introns.size());
		
		Pair<Collection<String>> filtered=filterExonIntron(introns, exons, data);
		introns=filtered.getValue1();
		exons=filtered.getValue2();
		
		System.err.println("filtered "+exons.size()+" "+introns.size());
		
		/*rnaList.add("intron");
		rnaList.add("exon");*/
		
		for(String exon: exons){
			rnaList.add(exon+".exon");
		}
		for(String intron: introns){
			rnaList.add(intron+".intron");
		}
		
		//rnas.put("exon", "exon");
		//rnas.put("intron", "intron");
		//rnaList.add("intron");
		//rnaList.add("exon");
		
		
		
		MatrixWithHeaders countMatrix=new MatrixWithHeaders(rnaList, rnaList);
		MatrixWithHeaders countMatrixWeighted=new MatrixWithHeaders(rnaList, rnaList);
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			Cluster sub=c.subsetRNA(kmer, introns, exons);
			for(String gene1: sub.getRNANames()){
				String class1=parseClass(gene1, rnas);
				for(String gene2: sub.getRNANames()){
					String class2=parseClass(gene2, rnas);
					double score=countMatrix.get(class1, class2);
					score++;
					countMatrix.set(class1, class2, score);
					
					score=countMatrixWeighted.get(class1, class2);
					score+=(2.0/c.getClusterSize());
					countMatrixWeighted.set(class1, class2, score);
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
		
		MatrixWithHeaders median=average(countMatrix);
		median.write(save+".all.median");
		
		countMatrix.write(save+".all.counts");
		countMatrixWeighted.write(save+".all.weighted");
		
		Process p=Runtime.getRuntime().exec("python /groups/guttman/SPRITE2/ice_matrix.py "+save+".all.counts");
		p.waitFor();
		
		p=Runtime.getRuntime().exec("python /groups/guttman/SPRITE2/ice_matrix.py "+save+".all.weighted");
		p.waitFor();
		
		p=Runtime.getRuntime().exec("python /groups/guttman/SPRITE2/ice_matrix.py "+save+".all.median");
		p.waitFor();
	}
	

	private MatrixWithHeaders average(MatrixWithHeaders countMatrix) {
		List<String> rowNames=new ArrayList<String>();
		
		for(String name: countMatrix.getRowNames()){
			if(!name.contains("intron") && !name.contains("exon")){rowNames.add(name);}
		}
		
		rowNames.add("intron");
		rowNames.add("exon");
		
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rowNames, rowNames);
		
		
		for(String row: rowNames){
			for(String col: rowNames){
				if(countMatrix.containsRow(row) && countMatrix.containsColumn(col)){
					rtrn.set(row, col, countMatrix.get(row, col));
				}
			}
		}
			
		for(String row: countMatrix.getRowNames()){
			if(!row.contains("intron") && !row.contains("exon")){
				double val=getAllVals(countMatrix, row, "intron");
				rtrn.set(row, "intron", val);
				rtrn.set("intron", row, val);
				val=getAllVals(countMatrix, row, "exon");
				rtrn.set(row, "exon", val);
				rtrn.set("exon", row, val);
			}
		}
		
		return rtrn;
	}
	


	private double getAllVals(MatrixWithHeaders countMatrix, String row, String string) {
		ArrayList<Double> vals=new ArrayList<Double>();
		for(String column: countMatrix.getColumnNames()){
			if(column.contains(string)){vals.add(countMatrix.get(row, column));}
		}
		return Statistics.quantile(vals, 0.5);
	}


	private Pair<Collection<String>> filterExonIntron(Collection<String> introns, Collection<String> exons, BarcodingDataStreaming data) {
		Map<String, Integer> intronCounts=new TreeMap<String, Integer>();
		Map<String, Integer> exonCounts=new TreeMap<String, Integer>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			for(RNAInterval gene: c.getAllRNARegions()){
				String name=gene.getName();
				String type=gene.getType();
				if(type.equals("intron")){
					int count=0;
					if(intronCounts.containsKey(name)){count=intronCounts.get(name);}
					count++;
					intronCounts.put(name, count);
				}
				if(type.equals("exon")){
					int count=0;
					if(exonCounts.containsKey(name)){count=exonCounts.get(name);}
					count++;
					exonCounts.put(name, count);
				}
				
			}
		}
		
		data.close();
		
		Collection<String> v1=new TreeSet<String>();
		for(String gene: intronCounts.keySet()){
			int count=intronCounts.get(gene);
			if(count>minCount && introns.contains(gene)){v1.add(gene);}
		}
		
		
		
		
		Collection<String> v2=new TreeSet<String>();
		for(String gene: exonCounts.keySet()){
			int count=exonCounts.get(gene);
			if(count>minCount && exons.contains(gene)){v2.add(gene);}
		}
		
		Pair<Collection<String>> rtrn=new Pair<Collection<String>>(v1, v2);
		return rtrn;
	}
	
	private Collection<String> computeExonIntronCounts(BarcodingDataStreaming data, Collection<String> mRNAs) {
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			for(String gene: c.getRNANames()){
				int count=0;
				if(counts.containsKey(gene)){count=counts.get(gene);}
				count++;
				counts.put(gene, count);
			}
		}
		
		data.close();
		
		Collection<String> rtrn=new TreeSet<String>();
		for(String gene: counts.keySet()){
			int count=counts.get(gene);
			if(count>minCount && mRNAs.contains(gene)){rtrn.add(gene);}
		}
		
		System.err.println(counts.size()+" "+rtrn.size());
		
		return rtrn;
	}


	private String parseClass(String gene1, Map<String, String> rnas) {
		if(rnas.containsKey(gene1)){return rnas.get(gene1);}
		//System.err.println(gene1);
		//return gene1.split("\\.")[1];
		return gene1;
	}


	public RNAContacts(BarcodingDataStreaming data, String save, Map<String, String> rnas) throws IOException, InterruptedException{
		//TODO Compute input
		
		Collection<String> names=new TreeSet<String>();
		
		Kmer kmer=new Kmer();
		List<String> rnaList=new ArrayList<String>();
		for(String rna: rnas.keySet()){
			kmer.addRegion(rna);
			String rnaClass=rnas.get(rna);
			if(!rnaList.contains(rnaClass)){rnaList.add(rnaClass);}
		}
		
		MatrixWithHeaders countMatrix=new MatrixWithHeaders(rnaList, rnaList);
		MatrixWithHeaders countMatrixWeighted=new MatrixWithHeaders(rnaList, rnaList);
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			names.addAll(c.getRNANames());
			Cluster sub=c.subsetRNA(kmer);
			for(RNAInterval gene1: sub.getAllRNARegions()){
				String class1=rnas.get(gene1.getName());
				for(RNAInterval gene2: sub.getAllRNARegions()){
					if(!gene1.equals(gene2)) {
						String class2=rnas.get(gene2.getName());
						double score=countMatrix.get(class1, class2);
						score++;
						countMatrix.set(class1, class2, score);
						
						score=countMatrixWeighted.get(class1, class2);
						score+=(2.0/c.getClusterSize());
						countMatrixWeighted.set(class1, class2, score);
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
		
		countMatrix.write(save+".all.counts");
		countMatrixWeighted.write(save+".all.weighted");
		
		Process p=Runtime.getRuntime().exec("python /groups/guttman/SPRITE2/ice_matrix.py "+save+".all.counts");
		p.waitFor();
		
		p=Runtime.getRuntime().exec("python /groups/guttman/SPRITE2/ice_matrix.py "+save+".all.weighted");
		p.waitFor();
		
		
		
	}
	
	private boolean hasDNA(Cluster c) {
		if(!c.getAllDNAIntervals().isEmpty() && c.getAllDNAIntervals().size()>0){return true;}
		return false;
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


	private boolean useCluster(Cluster c, String type) {
		if(type.equalsIgnoreCase("DNA")){
			//only use if contains DNA
			if(!c.getAllDNAIntervals().isEmpty() && c.getAllDNAIntervals().size()>0){return true;}
			return false;
		}
		else if(type.equalsIgnoreCase("RNA")){
			//only use if does NOT contain DNA (RNA only)
			if(c.getAllDNAIntervals().isEmpty()){return true;}
			return false;
		}
		return true;
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


	private double sumOfSquares(BarcodingDataStreaming data, String type) {
		double rtrn=0;
		
		while(data.hasNext()){
			Cluster c=data.next();
			if(this.useCluster(c, type)){
				rtrn+=c.getClusterSize();
			}
		}
		
		data.close();
		
		return rtrn;
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

	private List<String> getNames(BarcodingDataStreaming data) {
		Collection<String> rtrn=new TreeSet<String>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			rtrn.addAll(c.getRNANames());
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		
		List<String> list=new ArrayList<String>();
		list.addAll(rtrn);
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
			if(args.length>3){
				Collection<String> introns=parseMrna(args[3]);
				Collection<String> exons=parseMrna(args[4]);
				new RNAContacts(data, save, genes, introns, exons);
			}
			else{
				new RNAContacts(data, save, genes);
			}
		}
		else{System.err.println(usage);}
	}
	
	
	private static Collection<String> parseMrna(String string) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		
		List<String> list=BEDFileIO.loadLines(string);
		
		for(String line: list){
			String name=line.split("\t")[0];
			//System.err.println(name);
			rtrn.add(name);
		}
		
		return rtrn;
	}


	private static Map<String, String> makeGenes(BarcodingDataStreaming data) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			for(String rna: c.getRNANames()){rtrn.put(rna, rna);}
		}
		data.close();
		
		return rtrn;
	}


	static String usage=" args[0]=clusters \n args[1]=save \n args[2]=list of genes";
	
}

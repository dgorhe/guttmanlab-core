package guttmanlab.core.sars;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.ScanStat;
import guttmanlab.core.math.Statistics;
import jsc.distributions.Binomial;

public class EnrichmentMatrixForPrecomputedBins {

	Map<String, Integer> counts;
	Map<String, Integer> totalCounts;
	//int totalCounts;
	long totalSize;
	
	static double minCount=10;
	static double minEnrichment=2;
	static double alpha=0.001;
	int binSize;
	
	public EnrichmentMatrixForPrecomputedBins(File bam) throws NumberFormatException, IOException{
		this.counts=score(bam);
		this.binSize=getBinSize(bam);
		this.totalCounts=getTotalCounts(bam);
		this.totalSize=getTotalSize(bam);
	}
	
	private int getBinSize(File bam) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(bam)));
		String nextLine=reader.readLine();
		reader.close();
		String[] tokens=nextLine.split("\t");
		int val=new Integer(tokens[3]);
		return val;
	}

	private Map<String, Integer> getTotalCounts(File bam) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(bam)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String key=tokens[0];
			int val=new Integer(tokens[2]);
			rtrn.put(key, val);
		}
		reader.close();
		return rtrn;
	}

	private Map<String, Integer> score(File bam) throws NumberFormatException, IOException{
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(bam)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String key=tokens[0];
			int val=new Integer(tokens[1]);
			rtrn.put(key, val);
		}
		reader.close();
		return rtrn;
	}
	
	
	private static MatrixWithHeaders makeMatrix(Map<String, Integer>[] maps, Map<String, Integer>[] totalCounts2, File[] files, EnrichmentMatrixForPrecomputedBins input, double minEnrichment) {
		long totalSize=input.totalSize;
		List<String> columns=getColumnNames(files);
		List<String> rows=getRowNames(maps);
		
		int binSize=input.binSize;
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(int i=0; i<files.length; i++){
			String column=files[i].getName();
			for(String row: maps[i].keySet()){
				int totalCount=totalCounts2[i].get(row);
				int geneCount=maps[i].get(row);
				int inputCount=1;
				if(input.counts.containsKey(row)){inputCount=input.counts.get(row);}
				
				
				if(input.totalCounts.containsKey(row)){
					int inputTotalCounts=input.totalCounts.get(row);
					double p=getPValue(geneCount, inputCount, totalCount, inputTotalCounts);
					double ratio=getRatio(geneCount, inputCount, totalCount, inputTotalCounts);
				
				
					if(p<alpha && ratio>minEnrichment && geneCount>minCount){
						double scanP=getScanPValue(geneCount, inputCount, totalCount, inputTotalCounts, binSize, totalSize);
						//System.err.println(row+" "+column+" "+ratio+" "+p+" "+scanP);
						if(scanP<alpha){
							rtrn.set(row, column, ratio);
						}
						
					}
				}
				else{
					System.err.println("skipped "+row);
				}
			}
		}
		
		return rtrn;
	}
	
	
	private static long getTotalSize(File file) throws IOException {
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine=reader.readLine();
		reader.close();
		String[] tokens=nextLine.split("\t");
		int val=new Integer(tokens[4]);
		return val;
	}
	
	
	private static double getRatio(int sampleCounts, int inputCounts, int totalSampleCounts, int totalInputCounts) {
		double ratio=((double)sampleCounts/(double)totalSampleCounts)/((double)inputCounts/(double)totalInputCounts);
		return ratio;
	}

	private static double getPValue(int sampleCounts, int inputCounts, int totalSampleCounts, int totalInputCounts) {
		//Return binomial p
		
		int n=new Double(sampleCounts+inputCounts).intValue();
		if(n>0){
			double p=(double)totalSampleCounts/((double)totalSampleCounts+(double)totalInputCounts);
			Binomial b=new Binomial(n, p);
			return 1-b.cdf(sampleCounts);
		}
		return 1.0;
	}
	
	public static double getScanPValue(int sampleCount, int inputCount, int sampleTotalCount, int inputTotalCounts, int windowSize, long totalSize) {
		return ScanStat.getPValue(inputCount, sampleCount, inputTotalCounts, sampleTotalCount, windowSize, totalSize);
	}
	
	private static List<String> getRowNames(Map<String, Integer>[] maps) {
		Collection<String> genes=new TreeSet<String>();
		
		for(int i=0; i<maps.length; i++){
			genes.addAll(maps[i].keySet());
		}
		
		List<String> rtrn=new ArrayList<String>();
		rtrn.addAll(genes);
		return rtrn;
	}

	private static List<String> getColumnNames(File[] files) {
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<files.length; i++){
			rtrn.add(files[i].getName());
		}
		
		return rtrn;
	}

	
	private static MatrixWithHeaders collapseByGene(MatrixWithHeaders mwh) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(getGenes(mwh.getRowNames()), mwh.getColumnNames());
		
		for(String region: mwh.getRowNames()){
			String gene=region.split(":")[0];
			double[] vals=mwh.getRow(region);
			double[] vals2=rtrn.getRow(gene);
			double[] max=max(vals, vals2);
			rtrn.setRow(gene, max);
		}
		
		return rtrn;
	}
	
	private static List<String> getGenes(List<String> rowNames) {
		Collection<String> names=new TreeSet<String>();
		
		for(String name: rowNames){
			names.add(name.split(":")[0]);
		}
		
		List<String> rtrn=new ArrayList<String>();
		rtrn.addAll(names);
		return rtrn;
	}

	private static double[] max(double[] vals, double[] vals2) {
		double[] rtrn=new double[vals.length];
		
		for(int i=0; i<vals.length; i++){
			rtrn[i]=Math.max(vals[i], vals2[i]);
		}
		
		return rtrn;
	}

	private static MatrixWithHeaders filter(MatrixWithHeaders mwh) {
		List<String> genesToInclude=new ArrayList<String>();
		for(String gene: mwh.getRowNames()){
			double[] vals=mwh.getRow(gene);
			if(Statistics.sum(vals)>0){genesToInclude.add(gene);}
		}
		return mwh.submatrixByRowNames(genesToInclude);
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		File inputFile=new File(args[2]);
		double minEnrichment=new Double(args[3]);	
	
		Map<String, Integer>[] totalCounts=new Map[files.length];
		Map<String, Integer>[] maps=new Map[files.length];
		for(int i=0; i<files.length ; i++){
			System.err.println(files[i].getName());
			EnrichmentMatrixForPrecomputedBins m=new EnrichmentMatrixForPrecomputedBins(files[i]);
			totalCounts[i]=m.totalCounts;
			maps[i]=m.counts;
		}
		
		EnrichmentMatrixForPrecomputedBins input=new EnrichmentMatrixForPrecomputedBins(inputFile);
		MatrixWithHeaders mwh=makeMatrix(maps, totalCounts, files, input, minEnrichment);
		mwh=filter(mwh);
		
		MatrixWithHeaders collapsed=collapseByGene(mwh);
		collapsed.write(save+".collapsed.matrix");
		
		
		mwh.write(save);
		}
		else{System.err.println(usage);}
	}

	

	



	static String usage=" args[0]=all files \n args[1]=save \n args[2]=control file \n args[3]=min enrichment";
	

	
	
	
}

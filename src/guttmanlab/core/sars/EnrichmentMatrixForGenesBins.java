package guttmanlab.core.sars;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.ScanStat;
import guttmanlab.core.math.Statistics;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;

public class EnrichmentMatrixForGenesBins {

	Map<String, Integer> counts;
	int totalCounts;
	
	static double minCount=10;
	static double minEnrichment=2;
	static double alpha=0.001;
	
	public EnrichmentMatrixForGenesBins(File bam, int binSize){
		score(bam, binSize);
	}
	
	private TreeMap<String, Integer> score(File bam, int binSize){
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<String, Integer> positionCount=new TreeMap<String, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				String bin=bin(read, binSize);
				
				int score=0;
				if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
				score++;
				positionCount.put(bin, score);
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		this.totalCounts=totalCount;
		this.counts=positionCount;
		
		return positionCount;
	}
	
	
	
	
	private String bin(SAMRecord read, int binSize) {
		int startIndex=read.getAlignmentStart()/binSize;
		int newStart=startIndex*binSize;
		int newEnd=newStart+binSize;
		SingleInterval newInterval=new SingleInterval(read.getReferenceName(), newStart, newEnd);
		return newInterval.toUCSC();
	}

	private static MatrixWithHeaders makeMatrix(Map<String, Integer>[] maps, int[] totalCounts2, File[] files, int binSize, long totalSize) {
		List<String> columns=getColumnNames(files);
		List<String> rows=getRowNames(maps);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		int sumTotalCounts=sum(totalCounts2);
		Map<String, Integer> sumGeneCounts=sum(maps);
		
		for(int i=0; i<files.length; i++){
			int totalCount=totalCounts2[i];
			String column=files[i].getName();
			for(String row: maps[i].keySet()){
				int geneCount=maps[i].get(row);
				double p=getPValue(geneCount, sumGeneCounts.get(row), totalCount, sumTotalCounts);
				double ratio=getRatio(geneCount, sumGeneCounts.get(row), totalCount, sumTotalCounts);
				if(p<alpha && ratio>minEnrichment && geneCount>minCount){
					double scanP=getScanPValue(geneCount, sumGeneCounts.get(row), totalCount, sumTotalCounts, binSize, totalSize);
					System.err.println(row+" "+column+" "+ratio+" "+p+" "+scanP);
					if(scanP<alpha){
						rtrn.set(row, column, ratio);
					}
				}
			}
		}
		
		
		return rtrn;
	}
	
	private static MatrixWithHeaders makeMatrix(Map<String, Integer>[] maps, int[] totalCounts2, File[] files, int binSize, long totalSize, EnrichmentMatrixForGenesBins input) {
		List<String> columns=getColumnNames(files);
		List<String> rows=getRowNames(maps);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		int sumTotalCounts=sum(totalCounts2);
		Map<String, Integer> sumGeneCounts=sum(maps);
		
		for(int i=0; i<files.length; i++){
			int totalCount=totalCounts2[i];
			String column=files[i].getName();
			for(String row: maps[i].keySet()){
				int geneCount=maps[i].get(row);
				double p=getPValue(geneCount, sumGeneCounts.get(row), totalCount, sumTotalCounts);
				double ratio=getRatio(geneCount, sumGeneCounts.get(row), totalCount, sumTotalCounts);
				if(p<alpha && ratio>minEnrichment && geneCount>minCount){
					double scanP=getScanPValue(geneCount, sumGeneCounts.get(row), totalCount, sumTotalCounts, binSize, totalSize);
					//System.err.println(row+" "+column+" "+ratio+" "+p+" "+scanP);
					if(scanP<alpha){
						int inputCount=1;
						if(input.counts.containsKey(row)){inputCount=input.counts.get(row);}
						double scanP2=getScanPValue(geneCount, inputCount, totalCount, input.totalCounts, binSize, totalSize);
						double p2=getPValue(geneCount, inputCount, totalCount, input.totalCounts);
						double ratio2=getRatio(geneCount, inputCount, totalCount, input.totalCounts);
						if(scanP2<alpha && p2<alpha && ratio2>minEnrichment){
							double minRatio=Math.max(ratio, ratio2);
							rtrn.set(row, column, minRatio);
						}
					}
				}
			}
		}
		
		
		return rtrn;
	}
	
	
	private static long getSize(File file) {
		SAMFileReader inputReader= new SAMFileReader(file);
		
		SAMSequenceDictionary dict=inputReader.getFileHeader().getSequenceDictionary();
					
		inputReader.close();
		return dict.getReferenceLength();
		
	}
	
	private static Map<String, Integer> sum(Map<String, Integer>[] maps) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(int i=0; i<maps.length; i++){
			for(String gene: maps[i].keySet()){
				int sum=0;
				if(rtrn.containsKey(gene)){sum=rtrn.get(gene);}
				sum+=maps[i].get(gene);
				rtrn.put(gene, sum);
			}
		}
		
		return rtrn;
	}

	private static int sum(int[] totalCounts2) {
		int sum=0;
		for(int i=0; i< totalCounts2.length; i++){sum+=totalCounts2[i];}
		return sum;
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
		if(args.length>2){
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		int binSize=new Integer(args[2]);
		File inputFile=new File(args[3]);
		
		long totalSize=getSize(files[0]);
		
		int[] totalCounts=new int[files.length];
		Map<String, Integer>[] maps=new Map[files.length];
		for(int i=0; i<files.length ; i++){
			System.err.println(files[i].getName());
			EnrichmentMatrixForGenesBins m=new EnrichmentMatrixForGenesBins(files[i], binSize);
			totalCounts[i]=m.totalCounts;
			maps[i]=m.counts;
		}
		
		EnrichmentMatrixForGenesBins input=new EnrichmentMatrixForGenesBins(inputFile, binSize);
		MatrixWithHeaders mwh=makeMatrix(maps, totalCounts, files, binSize, totalSize, input);
		mwh=filter(mwh);
		
		MatrixWithHeaders collapsed=collapseByGene(mwh);
		collapsed.write(save+".collapsed.matrix");
		
		//TODO write out intermediates to quickly compute enrichments
		
		
		mwh.write(save);
		}
		else{System.err.println(usage);}
	}

	

	



	static String usage=" args[0]=bam files \n args[1]=save \n args[2]=bin size";
	

	
	
	
}

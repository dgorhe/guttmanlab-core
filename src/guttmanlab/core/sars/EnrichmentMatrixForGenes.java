package guttmanlab.core.sars;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import jsc.distributions.Binomial;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class EnrichmentMatrixForGenes {

	Map<String, Integer> counts;
	int totalCounts;
	
	static double minCount=10;
	static double minEnrichment=2;
	static double alpha=0.001;
	
	public EnrichmentMatrixForGenes(File bam){
		score(bam);
	}
	
	private TreeMap<String, Integer> score(File bam){
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<String, Integer> positionCount=new TreeMap<String, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				String gene=read.getReferenceName();
				
				int score=0;
				if(positionCount.containsKey(gene)){score=positionCount.get(gene);}
				score++;
				positionCount.put(gene, score);
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
	
	private static MatrixWithHeaders makeMatrix(Map<String, Integer>[] maps, int[] totalCounts2, File[] files) {
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
				//rtrn.set(row, column, geneCount);
				double p=getPValue(geneCount, sumGeneCounts.get(row), totalCount, sumTotalCounts);
				double ratio=getRatio(geneCount, sumGeneCounts.get(row), totalCount, sumTotalCounts);
				//if(p<alpha && ratio>minEnrichment && geneCount>minCount){
					rtrn.set(row, column, ratio);
				//}
			}
		}
		
		
		return rtrn;
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

	
	private static MatrixWithHeaders filter(MatrixWithHeaders mwh) {
		List<String> genesToInclude=new ArrayList<String>();
		for(String gene: mwh.getRowNames()){
			double[] vals=mwh.getRow(gene);
			if(Statistics.sum(vals)>0){genesToInclude.add(gene);}
		}
		return mwh.submatrixByRowNames(genesToInclude);
	}
	
	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		
		int[] totalCounts=new int[files.length];
		Map<String, Integer>[] maps=new Map[files.length];
		for(int i=0; i<files.length ; i++){
			System.err.println(files[i].getName());
			EnrichmentMatrixForGenes m=new EnrichmentMatrixForGenes(files[i]);
			totalCounts[i]=m.totalCounts;
			maps[i]=m.counts;
		}
		MatrixWithHeaders mwh=makeMatrix(maps, totalCounts, files);
		mwh=filter(mwh);
		
		mwh.write(save);
	}

	

	
	
	
}

package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;

public class IdentifySignificantRegions {
	
	int windowSize=4;
	private double enrichmentCutoff=5.0;
	private double minCount=2;
	private int minRNAFreq=10;
	private int binResolution=1000000;
	
	
	public IdentifySignificantRegions(BarcodingDataStreaming data, String save, int binResolution, double enrichmentCutoff) throws IOException{
		this.binResolution=binResolution;
		this.enrichmentCutoff=enrichmentCutoff;
		MatrixWithHeaders transSites=identifyTransSites(data);
		transSites.write(save);
	}
	
	private MatrixWithHeaders identifyTransSites(BarcodingDataStreaming data) throws IOException {
		Pair<MatrixWithHeaders> pairedScores=data.getRNADNAContactMatrix(minRNAFreq, binResolution);
		
		MatrixWithHeaders scores=pairedScores.getValue1();
		MatrixWithHeaders counts=pairedScores.getValue2();
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(scores.getRowNames(), scores.getColumnNames());
		
		Collection<String> genes=new TreeSet<String>();
		
		for(String gene: scores.getRowNames()){
			double lambda=mean(scores, gene);
			Collection<String> regions=scanWindows(scores, counts, gene, lambda);
			if(!regions.isEmpty()){
				System.err.println(gene+" "+regions.size());
				update(scores, gene, regions, rtrn, lambda);
				genes.add(gene);
			}
			
			
			
			//else{System.err.println("Skipping "+gene+" "+val);}
		}
		rtrn=rtrn.submatrixByRowNames(genes);
		return rtrn;
	}

	public IdentifySignificantRegions(MatrixWithHeaders scores, MatrixWithHeaders counts, String save) throws IOException{
		
		Collection<String> genes=new TreeSet<String>();
		MatrixWithHeaders mwh=new MatrixWithHeaders(scores.getRowNames(), scores.getColumnNames());
		
		for(String gene: scores.getRowNames()){
			//System.err.println(gene);
			double lambda=mean(scores, gene);
			Collection<String> regions=scanWindows(scores, counts, gene, lambda);
			if(!regions.isEmpty()){
				System.err.println(gene+" "+regions.size());
				update(scores, gene, regions, mwh, lambda);
				genes.add(gene);
			}
			
			
			
			//else{System.err.println("Skipping "+gene+" "+val);}
		}
		mwh=mwh.submatrixByRowNames(genes);
		mwh.write(save);
	}
	
	
	private void update(MatrixWithHeaders scores, String gene, Collection<String> regions, MatrixWithHeaders mwh, double lambda) throws IOException {
		
		for(String genome: regions){
			double val=scores.get(gene, genome);
			double enrichment=val/lambda;
			mwh.set(gene, genome, enrichment);
		}
		
	}
	
	private void write(String save, MatrixWithHeaders scores, String gene, Collection<String> regions) throws IOException {
		FileWriter writer=new FileWriter(save+"."+gene+".bed");
		
		for(String region: regions){
			SingleInterval genomeRegion=new SingleInterval(region);
			writer.write(genomeRegion.getReferenceName()+"\t"+genomeRegion.getReferenceStartPosition()+"\t"+genomeRegion.getReferenceEndPosition()+"\n");
		}
		
		writer.close();
		
		writer=new FileWriter(save+"."+gene+".bedgraph");
		
		for(String genome: regions){
			SingleInterval genomeRegion=new SingleInterval(genome);
			writer.write(genomeRegion.getReferenceName()+"\t"+genomeRegion.getReferenceStartPosition()+"\t"+genomeRegion.getReferenceEndPosition()+"\t"+scores.get(gene, genome)+"\n");
		}
		writer.close();
	}

	private Collection<String> scanWindows(MatrixWithHeaders scores, MatrixWithHeaders counts, String gene, double lambda) {
		Collection<String> rtrn=new TreeSet<String>();
		
		//System.err.println(gene+" "+lambda);
		
		Iterator<Collection<String>> windowIter=getWindows(scores, gene, windowSize);
		
		while(windowIter.hasNext()){
			Collection<String> columns=windowIter.next();
			boolean sig=isSig(scores, counts, gene, columns, lambda);
			if(sig){
				rtrn.addAll(columns);
			}
		}
		return rtrn;
	}

	private double mean(MatrixWithHeaders scores, String gene) {
		List<Double> vals=new ArrayList<Double>();
		
		SingleInterval rnaRegion=new SingleInterval(scores.getPIDToName().get(gene));
		
		for(String column: scores.getColumnNames()){
			SingleInterval genomeRegion=new SingleInterval(column);
			if(genomeRegion.getReferenceName().equals(rnaRegion.getReferenceName())){
				double val=scores.get(gene, column);
				vals.add(val);
			}
		}
		
		return Statistics.mean(vals);
	}

	private Iterator<Collection<String>> getWindows(MatrixWithHeaders scores, String gene, int windowSize2) {
		Collection<Collection<String>> rtrn=new ArrayList<Collection<String>>();
		//System.err.println(gene);
		//System.err.println(scores.getPIDToName().get(gene));
		SingleInterval rnaRegion=new SingleInterval(scores.getPIDToName().get(gene));
		List<String> columnNames=scores.getColumnNames();
		for(int i=0; i<columnNames.size()-windowSize2; i++){
			Collection<String> set=get(columnNames, i, i+windowSize2, rnaRegion);
			if(set!=null){
				rtrn.add(set);
			}
		}
		return rtrn.iterator();
	}

	private Collection<String> get(List<String> columnNames, int start, int end, SingleInterval rnaRegion) {
		Collection<String> rtrn=new ArrayList<String>();
		
		Collection<String> chromosomes=new TreeSet<String>();
		for(int i=start; i<=end; i++){
			String region=columnNames.get(i);
			rtrn.add(region);
			if(region.contains("Input")){chromosomes.add(region);}
			else{
				SingleInterval genome=new SingleInterval(region);
				chromosomes.add(genome.getReferenceName());
				//if(genome.getReferenceName().equals(rnaRegion.getReferenceName())){return null;}
			}
		}
		
		if(chromosomes.size()==1){return rtrn;}
		return null;
	}

	private boolean isSig(MatrixWithHeaders scores, MatrixWithHeaders counts, String gene, Collection<String> columns, double mean) {
		double lambda=mean*windowSize;
		double score=0;
		int numGreater=0;
		int countGreater=0;
		
		for(String column: columns){
			double val=scores.get(gene, column);
			double count=scores.get(gene, column);
			score+=val;
			double localEnrich=val/mean;
			if(localEnrich>enrichmentCutoff){numGreater++;}
			if(count>minCount){countGreater++;}
		}
		
		//double p=ScanStat.poissonCDF(score, lambda);
		double enrichment=score/lambda;
		
		/*if(numGreater>1){
			System.err.println(getRegion(columns).toUCSC()+" "+score+" "+lambda+" "+numGreater+" "+enrichment);
		}*/
		
		//return(p<cutoff && numGreater>1);
		return (enrichment>enrichmentCutoff && countGreater>3);
		
		//return (numGreater>3 && enrichment>enrichmentCutoff);
	}

	private SingleInterval getRegion(Collection<String> columns) {
		String chr="";
		int start=Integer.MAX_VALUE;
		int end=-Integer.MAX_VALUE;
		for(String genome: columns){
			SingleInterval region=new SingleInterval(genome);
			chr=region.getReferenceName();
			start=Math.min(start, region.getReferenceStartPosition());
			end=Math.max(end, region.getReferenceEndPosition());
		}
		
		return new SingleInterval(chr, start, end);
	}

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		int binResolution=new Integer(args[2]);
		double enrichmentCutoff=new Double(args[3]);
		new IdentifySignificantRegions(data, save, binResolution, enrichmentCutoff);
	}
	
}

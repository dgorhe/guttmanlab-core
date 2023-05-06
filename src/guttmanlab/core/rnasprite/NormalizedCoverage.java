package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class NormalizedCoverage {
	
	int minClusterSize=2;
	int maxClusterSize=1000;
	Map<String, SingleInterval> rnaRegions;
	int minRNAFreq=10;
	int binResolution;
	//double cutoff=2.0;
	
	static SingleInterval problemBin1=new SingleInterval("chr2:79490000-79500000");
	static SingleInterval problemBin2=new SingleInterval("chr11:3119270-3192250");
	static SingleInterval problemBin3=new SingleInterval("chr15:99734977-99736026");
	static SingleInterval problemBin4=new SingleInterval("chr3:5173978-5175025");
	static SingleInterval problemBin5=new SingleInterval("chr13:58176952-58178051");
	
	static Collection<SingleInterval> problemBins=new TreeSet<SingleInterval>();
	
	Map<String, Integer> clusterCount;
	double enrichment;
	
	
	
	public NormalizedCoverage(BarcodingDataStreaming data, String save, int binResolution, Collection<String> genesToUse, double enrichment) throws IOException{
		this.enrichment=enrichment;
		this.binResolution=binResolution;
		this.clusterCount=new TreeMap<String, Integer>();
		this.rnaRegions=new TreeMap<String, SingleInterval>();
		System.err.println("Started");
		MatrixWithHeaders mwh= computeInput(data, binResolution);
		
		System.err.println("Completed input "+mwh.getNumberRows()+" "+mwh.getNumberColumns());
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()<maxClusterSize){
				Collection<String> rnas=c.getRNANames();
				score(mwh, rnas, c, binResolution);
				counter++;
			}
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		mwh=filter(mwh, clusterCount);
		mwh=excludeFiltered(mwh, problemBins);
		
		
		//mwh=filter(mwh);
		
		mwh.write(save+".original.scores");
		writeAnnotation(save+".annotation", mwh);
		writeGeneAnnotation(save+".gene.annotation", clusterCount, rnaRegions);
		
		
		
		MatrixWithHeaders subset=mwh.submatrixByRowNames(genesToUse);
		subset.write(save+".subset.scores");
		
		
		
		MatrixWithHeaders normalized=normalize(subset);
		normalized.write(save+".normalized.scores");
		
		MatrixWithHeaders binary=discretize(normalized);
		binary.write(save+".binary.scores");
		
		MatrixWithHeaders normalizedTotal=normalize(subset, clusterCount);
		normalizedTotal.write(save+".normalizedRelativeToTotal.scores");
	}

	
	private MatrixWithHeaders discretize(MatrixWithHeaders mwh) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()){
			for(String column: mwh.getColumnNames()){
				double val=mwh.get(row, column);
				if(val>0){rtrn.set(row, column, 1.0);}
			}
		}
		
		return rtrn;
	}


	private MatrixWithHeaders normalize(MatrixWithHeaders mwh, Map<String, Integer> clusterCount) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()){
			double lambda=clusterCount.get(row);
			for(String column: mwh.getColumnNames()){
				double val=mwh.get(row, column);
				double enrichment=val/lambda;
				if(val>=5){rtrn.set(row, column, enrichment);}
				//if(val>=2 && enrichment>this.enrichment){rtrn.set(row, column, 1.0);}
			}
		}
		
		return rtrn;
	}
	
	private MatrixWithHeaders normalize(MatrixWithHeaders mwh) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()){
			double lambda=Statistics.mean(mwh.getRow(row));
			for(String column: mwh.getColumnNames()){
				double val=mwh.get(row, column);
				double enrichment=val/lambda;
				if(val>=5){rtrn.set(row, column, enrichment);}
				//if(val>=2 && enrichment>this.enrichment){rtrn.set(row, column, 1.0);}
			}
		}
		
		return rtrn;
	}
	
	


	private void writeGeneAnnotation(String save, Map<String, Integer> clusterCount2, Map<String, SingleInterval> rnaRegions2) throws IOException {
		FileWriter writer=new FileWriter(save);
	
		Map<SingleInterval, String> flip=flip(rnaRegions2);
		
		writer.write("ID\tCount\tPosition\tchromosome\torder\n");
		int order=0;
		for(SingleInterval region: flip.keySet()){
			String geneName=flip.get(region);
			int count=clusterCount2.get(geneName);
			writer.write(geneName+"\t"+count+"\t"+region.toUCSC()+"\t"+region.getReferenceName()+"\t"+order+"\n");
			order++;
		}
		
		
		
		
		writer.close();
	}


	private Map<SingleInterval, String> flip(Map<String, SingleInterval> rnaRegions2) {
		Map<SingleInterval, String> rtrn=new TreeMap<SingleInterval, String>();
		
		for(String name: rnaRegions2.keySet()){
			SingleInterval region=rnaRegions2.get(name);
			rtrn.put(region, name);
		}
		
		return rtrn;
	}


	private MatrixWithHeaders excludeFiltered(MatrixWithHeaders mwh, Collection<SingleInterval> problemBins) {
		Collection<String> list=new ArrayList<String>();
		for(SingleInterval problemBin: problemBins){
			String region=problemBin.bin(binResolution).toUCSC();
			list.add(region);
		}
		
		return mwh.excludeByColumnNames(list);
	}


	/*private MatrixWithHeaders getSignificantRegions(MatrixWithHeaders mwh) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String rna: mwh.getRowNames()){
			double percentile=getPercentile(mwh.getRow(rna), 0.9);
			for(String region: mwh.getColumnNames()){
				double val=mwh.get(rna, region);
				double ratio=val/percentile;
				if(ratio>cutoff){
					rtrn.set(rna, region, 1.0);
				}
			}
		}
		
		return rtrn;
	}*/


	private double getPercentile(double[] row, double d) {
		return Statistics.quantile(row, d);
	}


	private void write(String save, MatrixWithHeaders regionsByRNA) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String rna: regionsByRNA.getRowNames()){
			writer.write(rna);
			for(String region: regionsByRNA.getColumnNames()){
				double val=regionsByRNA.get(rna, region);
				if(val==1.0){
					writer.write("\t"+region);
				}
			}
			writer.write("\n");
		}
		
		writer.close();
	}


	private MatrixWithHeaders filter(MatrixWithHeaders mwh, Map<String, Integer> clusterCount) {
		Collection<String> genesToUse=new TreeSet<String>();
		
		for(String gene: mwh.getRowNames()){
			double val=clusterCount.get(gene);
			if(val>this.minRNAFreq){genesToUse.add(gene);}
		}
		
		mwh=mwh.submatrixByRowNames(genesToUse);
		return mwh;
	}


	/*private MatrixWithHeaders excludeCis(MatrixWithHeaders mwh) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String gene: mwh.getRowNames()){
			String chr=this.rnaRegions.get(gene).getReferenceName();
			for(String genome: mwh.getColumnNames()){
				if(genome.contains("Input")){
					rtrn.set(gene, genome, mwh.get(gene, genome));
				}
				else{
					SingleInterval genomeRegion=new SingleInterval(genome);
					if(!genomeRegion.getReferenceName().equals(chr)){
						rtrn.set(gene, genome, mwh.get(gene, genome));
					}
					else{
						rtrn.set(gene, genome, Double.NaN);
					}
				}
			}
		}
		
		return rtrn;
	}*/


	/*private MatrixWithHeaders excludeLocus(MatrixWithHeaders mwh) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String gene: mwh.getRowNames()){
			SingleInterval geneRegion=extend(this.rnaRegions.get(gene));
			for(String genome: mwh.getColumnNames()){
				SingleInterval genomeRegion=new SingleInterval(genome);
				if(!geneRegion.overlaps(genomeRegion)){
					rtrn.set(gene, genome, mwh.get(gene, genome));
				}
				else{
					rtrn.set(gene, genome, Double.NaN);
				}
			}
		}
		
		return rtrn;
	}*/


	private SingleInterval extend(SingleInterval region) {
		int numExtendedBins=2;
		int extension=numExtendedBins*binResolution;
		SingleInterval rtrn=new SingleInterval(region.getReferenceName(), region.getReferenceStartPosition()-extension, region.getReferenceEndPosition()+extension);
		return rtrn;
	}


	private MatrixWithHeaders normalizeByExpression(MatrixWithHeaders mwh) {
		List<String> rows=new ArrayList<String>();
		
		for(String gene: mwh.getRowNames()){
			double inputExpression=mwh.get(gene, "InputChromatin");
			if(inputExpression>this.minRNAFreq){rows.add(gene);}
			//else{System.err.println("Skipping "+gene+" "+inputExpression);}
		}
		
		List<String> columns=new ArrayList<String>();
		columns.addAll(mwh.getColumnNames());
		columns.remove("Input");
		columns.remove("InputChromatin");
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(String gene: rows){
			double inputExpression=mwh.get(gene, "InputChromatin");
			//System.err.println(gene+" "+inputExpression);
			for(String genome: columns){
				double val=mwh.get(gene, genome)/inputExpression;
				rtrn.set(gene, genome, val);
			}
		}
		
		return rtrn;
	}


	/*private MatrixWithHeaders scaleNormalize(MatrixWithHeaders mwh) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String gene: mwh.getRowNames()){
			if(this.rnaRegions.containsKey(gene)){
				SingleInterval geneRegion=this.rnaRegions.get(gene);
				double maxScore=-1;
				for(String genome: mwh.getColumnNames()){
					if(!genome.contains("Input")){
						SingleInterval genomeRegion=new SingleInterval(genome);
						if(geneRegion.overlaps(genomeRegion)){
							maxScore=Math.max(maxScore, mwh.get(gene, genome));
						}
					}
				}
				for(String genome: mwh.getColumnNames()){
					if(maxScore>0){
						double norm=mwh.get(gene, genome)/maxScore;
						rtrn.set(gene, genome, norm);
					}
					//else{System.err.println("Skipping gene "+gene+" no peak to normalize to");}
				}
			}
			//else{System.err.println("skipped "+gene+" because no region associated");}
		}
		return rtrn;
	}*/


	private void writeAnnotation(String save, MatrixWithHeaders mwh) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("ID\tChromosome\n");
		for(String pos: mwh.getColumnNames()){
			writer.write(pos+"\t"+pos.split(":")[0]+"\n");
		}
		
		writer.close();
	}

	private MatrixWithHeaders normalizeByInput(MatrixWithHeaders mwh) {
		return mwh.divideColumnsWithConstants(mwh.getRow("Input"));
	}

	private List<String> getNames(Collection<String> rnas) {
		List<String> rtrn=new ArrayList<String>();
		
		rtrn.addAll(rnas);
		
		return rtrn;
	}

	private void score(MatrixWithHeaders mwh, Collection<String> rnas, Cluster c, int binResolution) {
		Cluster binned=c.bin(binResolution);
		//double score=(2.0/(double)c.getAllDNAIntervals().size());
		double score=1;
		
		for(SingleInterval region: binned.getAllDNAIntervals()){
			for(String rna: rnas){
				double updatedScore=mwh.get(rna, region.toUCSC())+score;
				mwh.set(rna, region.toUCSC(), updatedScore);
			}
		}
		
	}

	private void updateScores(MatrixWithHeaders mwh, Map<SingleInterval, Double> counts, Map<SingleInterval, Double> inputCounts, double percentile, Gene gene) {
		
		for(SingleInterval region: inputCounts.keySet()){
			double sampleCount=0.0;
			if(counts.containsKey(region)){
				sampleCount=counts.get(region);
			}
			double inputCount=inputCounts.get(region);
			double ratio=(sampleCount/inputCount);
			if(inputCount>percentile){
				mwh.set(gene.getName(), region.toUCSC(), ratio);
			}
		}	
	}

	private List<String> getPositions(Map<SingleInterval, Double> inputCounts) {
		List<String> rtrn=new ArrayList<String>();
		
		for(SingleInterval region: inputCounts.keySet()){
			rtrn.add(region.toUCSC());
		}
		
		return rtrn;
	}

	

	private MatrixWithHeaders computeInput(BarcodingDataStreaming data, int binResolution) {
		List<String> geneNames=new ArrayList<String>();
		List<String> genomePosition=new ArrayList<String>();
		Collection<SingleInterval> tempPosition=new TreeSet<SingleInterval>();
		//Collection<RNAInterval> tempRNA=new TreeSet<RNAInterval>();
		Collection<String> rnaSet=new TreeSet<String>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			Cluster binned=c.bin(binResolution);
			tempPosition.addAll(binned.getAllDNAIntervals());
			//rnaSet.addAll(c.getRNANames());
			
			for(String rna: c.getRNANames()){
				rnaSet.add(rna);
				int count=0;
				if(this.clusterCount.containsKey(rna)){count=this.clusterCount.get(rna);}
				count++;
				this.clusterCount.put(rna, count);
			}
			
			for(RNAInterval r: c.getAllRNARegions()){
				String name=r.getName();
				SingleInterval region=merge(name, r);
				this.rnaRegions.put(name, region);
			}
			
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		data.close();
		
		
		for(SingleInterval region: tempPosition){
			genomePosition.add(region.toUCSC());
		}
		
		geneNames.addAll(rnaSet);
		
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(geneNames, genomePosition);
		return mwh;
		
		
		/*
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Map<String, Double> geneExpression=new TreeMap<String, Double>();
		Map<String, Double> geneExpressionOnChromatin=new TreeMap<String, Double>();
		
		Collection<String> genes=new TreeSet<String>();
		
		int totalCount=0;
		int totalChromatinCount=0;
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			Collection<RNAInterval> rnaRegions=c.getAllRNARegions();
			updateCoordinates(rnaRegions);
			genes.addAll(c.getRNANames());
			Collection<String> rnas=c.getRNANames();
			for(String rna: rnas){
				double countTotal=0;
				if(geneExpression.containsKey(rna)){countTotal=geneExpression.get(rna);}
				countTotal++;
				geneExpression.put(rna, countTotal);
			
				if(!c.getAllDNAIntervals().isEmpty()){
					double countOnChromatin=0;
					if(geneExpressionOnChromatin.containsKey(rna)){countOnChromatin=geneExpressionOnChromatin.get(rna);}
					countOnChromatin++;
					totalChromatinCount++;
					geneExpressionOnChromatin.put(rna, countOnChromatin);
				}
				totalCount++;
			}
			
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				double count=0;
				if(rtrn.containsKey(region)){count=rtrn.get(region);}
				count+=(2.0/(double)binned.getAllDNAIntervals().size());
				rtrn.put(region, count);
			}
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		data.close();
		//double percentile=percentile(rtrn, 0.05);
		
		List<String> geneNames=getNames(genes);
		//List<String> geneNames=getNames(geneExpressionOnChromatin);
		geneNames.add("Input");
		List<String> genomePosition=getPositions(rtrn);
		genomePosition.add("Input");
		genomePosition.add("InputChromatin");
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(geneNames, genomePosition);
		
		mwh.set("Input", "Input", totalCount);
		mwh.set("Input", "InputChromatin", totalChromatinCount);
		
		for(SingleInterval region: rtrn.keySet()){
			double val=rtrn.get(region);
			//if(val<percentile){val=percentile;}
			mwh.set("Input", region.toUCSC(), val);
		}
		
		for(String rna: geneExpression.keySet()){
			mwh.set(rna, "Input", geneExpression.get(rna));
			if(geneExpressionOnChromatin.containsKey(rna)){
				mwh.set(rna, "InputChromatin", geneExpressionOnChromatin.get(rna));
			}
		}
		
		return mwh;*/
	}

	private List<String> getNames(Map<String, Double> geneExpressionOnChromatin) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String gene: geneExpressionOnChromatin.keySet()){
			double count=geneExpressionOnChromatin.get(gene);
			if(count>this.minRNAFreq){
				//System.err.println(gene);
				rtrn.add(gene);
			}
		}
		
		return rtrn;
	}

	
	private SingleInterval merge(String name, RNAInterval region) {
		SingleInterval other=null;
		if(this.rnaRegions.containsKey(name)){other=this.rnaRegions.get(name);}
		if(other==null){return new SingleInterval(region.getReferenceName(), region.getReferenceStartPosition(), region.getReferenceEndPosition());}
		SingleInterval merge=update(other, region);
		return merge;
	}
	
	/*private void updateCoordinates(Collection<RNAInterval> regions) {
		for(SingleInterval rna: regions){
			String name=rna.getName();
			SingleInterval updated=rna;
			if(this.rnaRegions.containsKey(name)){
				updated=update(rnaRegions.get(name), rna);
			}
			this.rnaRegions.put(name, updated);
		}
	}*/

	private SingleInterval update(SingleInterval rna1, SingleInterval rna2) {
		String chr=rna1.getReferenceName();
		int start=Math.min(rna1.getReferenceStartPosition(), rna2.getReferenceStartPosition());
		int end=Math.max(rna1.getReferenceEndPosition(), rna2.getReferenceEndPosition());
		SingleInterval rtrn=new SingleInterval(chr, start, end);
		rtrn.setName(rna1.getName());
		return rtrn;
	}

	private void write(String save, Map<SingleInterval, Double> counts,  Map<SingleInterval, Double> input) throws IOException {
		FileWriter writer=new FileWriter(save+".sample.bedgraph");
		
		for(SingleInterval region: counts.keySet()){
			double count=counts.get(region);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+count+"\n");
		}
		
		writer.close();
		
		writer=new FileWriter(save+".input.bedgraph");
		
		for(SingleInterval region: input.keySet()){
			double count=input.get(region);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+count+"\n");
		}
		
		writer.close();
		
		writer=new FileWriter(save+".normalized.bedgraph");
		
		Map<SingleInterval, Double> normalized=normalize(counts, input);
		
		for(SingleInterval region: normalized.keySet()){
			double ratio=normalized.get(region);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+ratio+"\n");
		}
		
		writer.close();
	}
	
	private Map<SingleInterval, Double> normalize(Map<SingleInterval, Double> counts, Map<SingleInterval, Double> input) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		double percentile=percentile(input, 0.05);
		//double sum=sum(counts, input, percentile);
		double average=average(counts, input, percentile);
		System.err.println("percentile "+percentile+" average "+average);
		
		
		for(SingleInterval region: input.keySet()){
			double sampleCount=0.0;
			if(counts.containsKey(region)){
				sampleCount=counts.get(region);
			}
			double inputCount=input.get(region);
			//double ratio=100.0*((sampleCount/inputCount)/sum);
			double ratio=((sampleCount/inputCount)/average);
			//TODO Get 5% for cutoff
			//Normalize to total
			if(inputCount>percentile){
				rtrn.put(region, ratio);
			}
		}
		return rtrn;
	}

	private double percentile(Map<SingleInterval, Double> input, double pct) {
		Collection<Double> vals=input.values();
		List<Double> sorted=new ArrayList<Double>();
		sorted.addAll(vals);
		Collections.sort(sorted);
		return Statistics.quantile(sorted, pct);
	}

	private double sum(Map<SingleInterval, Double> counts, Map<SingleInterval, Double> input, double percentile) {
		double sum=0.0;
		
		for(SingleInterval region: input.keySet()){
			double sampleCount=0.0;
			if(counts.containsKey(region)){
				sampleCount=counts.get(region);
			}
			double inputCount=input.get(region);
			double ratio=sampleCount/inputCount;
			if(inputCount>percentile){
				sum+=ratio;
			}
		}
		
		return sum;
	}
	
	private double average(Map<SingleInterval, Double> counts, Map<SingleInterval, Double> input, double percentile) {
		double sum=0.0;
		double counter=0.0;
		
		for(SingleInterval region: input.keySet()){
			double sampleCount=0.0;
			if(counts.containsKey(region)){
				sampleCount=counts.get(region);
				counter+=1.0;
			}
			double inputCount=input.get(region);
			double ratio=sampleCount/inputCount;
			if(inputCount>percentile){
				sum+=ratio;
			}
		}
		
		return sum/counter;
	}

	/*private void write(String save, Map<SingleInterval, Double> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: counts.keySet()){
			double count=counts.get(region);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+count+"\n");
		}
		
		writer.close();
	}*/

	private Map<SingleInterval, Double> convertToCounts(Collection<Cluster> clusters, int binResolution) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		for(Cluster c: clusters){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				double count=0;
				if(rtrn.containsKey(region)){
					count=rtrn.get(region);
				}
				//count++;
				count+=(2.0/(double)c.getAllDNAIntervals().size());
				rtrn.put(region, count);
			}
		}
		return rtrn;
	}
	
	private static Collection<String> parse(String string) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			rtrn.add(line);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		problemBins.add(problemBin1);
		problemBins.add(problemBin2);
		problemBins.add(problemBin3);
		problemBins.add(problemBin4);
		problemBins.add(problemBin5);
		
		if(args.length>4){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String save=args[1];
			int resolution=new Integer(args[2]);
			Collection<String> list=parse(args[3]);
			double enrichment=new Double(args[4]);
			new NormalizedCoverage(data, save, resolution, list, enrichment);
		}
		else{System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=cluster file \n args[1]=save \n args[2]=resolution \n args[3]=gene list \n args[4]=enrichment";
}

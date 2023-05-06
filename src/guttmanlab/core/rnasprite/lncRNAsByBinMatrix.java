package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;

public class lncRNAsByBinMatrix {

	int minCount=3;
	int numPerm=100;
	int percentile=20;
	
	
	
	//TODO sample from all genomic alignments and count max overlaps, only retain those exceeding 5% of random
	//TODO Iterate over bin sizes?
	//TODO permute cluster sizes
	
	//TODO Use genes and set min bin size. If gene less than bin --> extend, if more than bin, break into bins
	
	public lncRNAsByBinMatrix(BarcodingDataStreaming data, int binResolution, String save, List<String> rnas, int minCount, Map<String, IntervalTree<String>> genes, int percentile) throws IOException {
		this.minCount=minCount;
		this.percentile=percentile;
		Map<String, Collection<SingleInterval>> coverageByRNA=new TreeMap<String, Collection<SingleInterval>>();
		
		
		Collection<SingleInterval> allRegions=new TreeSet<SingleInterval>();
		
		MatrixWithHeaders clusterCounts=new MatrixWithHeaders(rnas, CoordinateSpace.HG38.getBins(binResolution));
		
		//iterate over all clusters
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Cluster binned=c.bin(binResolution);
			Collection<String> rnaNames=c.getRNANames();
			allRegions.addAll(c.getAllDNAIntervals());
			
			for(String rna: rnaNames) {
				if(rnas.contains(rna)) {
					Collection<SingleInterval> set=getList(coverageByRNA, rna);
					Collection<SingleInterval> regions=c.getAllDNAIntervals();
					set.addAll(regions);
					Collection<SingleInterval> bins=binned.getAllDNAIntervals();
					for(SingleInterval bin: bins) {
						//System.err.println(rna+" "+bin.toUCSC());
						clusterCounts.incrementCount(rna, bin.toUCSC());
					}
				}
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		data.close();
		
		//write(coverageByRNA, save+".rnacoverage.list"); //RNA counts
		
		//clusterCounts.write(save+".clusterCounts.matrix");
		
		MatrixWithHeaders enrichment=normalizeByPermutations(allRegions, coverageByRNA, binResolution, numPerm, save, clusterCounts);
		
		
		List<String> rnasToUse=write(enrichment, clusterCounts, save+".list", genes);
		
		enrichment= enrichment.submatrixByRowNames(rnasToUse);
		enrichment.write(save+".enrichment.matrix");
	}
	
	
	


	private MatrixWithHeaders normalizeByPermutations(Collection<SingleInterval> allRegions, Map<String, Collection<SingleInterval>> coverageByRNA, int binResolution, int numPerm, String save, MatrixWithHeaders clusterCounts) throws IOException {
		List<SingleInterval> listRegions=new ArrayList<SingleInterval>();
		listRegions.addAll(allRegions);
		
		MatrixWithHeaders observed=makeMatrix(coverageByRNA, binResolution);
		//observed.write(save+".counts.matrix");
		
		MatrixWithHeaders maxRandom=new MatrixWithHeaders(observed.getRowNames(), observed.getColumnNames());
		for(int i=0; i<numPerm; i++) {
			System.err.println("perm "+i);
			MatrixWithHeaders random=sample(listRegions, coverageByRNA, binResolution);
			max(random, maxRandom);
		}
		
		//maxRandom.write(save+".maxRandom.matrix");
		
		MatrixWithHeaders enrichment=new MatrixWithHeaders(observed.getRowNames(), observed.getColumnNames());
		
		for(String row: observed.getRowNames()) {
			double max=Statistics.max(maxRandom.getRow(row));
			//System.out.println(row+"\t"+max);
			for(String column: observed.getColumnNames()) {
				double o=observed.get(row, column);
				double clusterCount=clusterCounts.get(row, column);
				if(clusterCount<minCount) {o=0.0;}
				enrichment.set(row, column, o/max);
			}
		}
		return enrichment;
	}





	private MatrixWithHeaders sample(List<SingleInterval> allRegions, Map<String, Collection<SingleInterval>> coverageByRNA, int binResolution) {
		Map<String, Collection<SingleInterval>> random=new TreeMap<String, Collection<SingleInterval>>();
		
		for(String rna: coverageByRNA.keySet()) {
			int numberOfRegions=coverageByRNA.get(rna).size();
			Collection<SingleInterval> randomList=sample(allRegions, numberOfRegions);
			random.put(rna, randomList);
		}
		
		MatrixWithHeaders rtrn=makeMatrix(random, binResolution);
		return rtrn;
	}

	private static Collection<SingleInterval> sample(List<SingleInterval> allRegions, int numberOfRegions) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		for(int i=0; i<numberOfRegions; i++) {
			int index=new Double(Math.random()*(allRegions.size()-1)).intValue();
			rtrn.add(allRegions.get(index));
		}
		return rtrn;
	}

	private void max(MatrixWithHeaders random, MatrixWithHeaders maxRandom) {
		for(String row: maxRandom.getRowNames()) {
			for(String col: maxRandom.getColumnNames()) {
				double val1=0;
				if(random.containsRow(row) && random.containsColumn(col)) {
					val1=random.get(row, col);
				}
				double val2=maxRandom.get(row, col);
				double max=Math.max(val1,  val2);
				maxRandom.set(row, col, max);
			}
		}
	}

	private void write(Map<String, Collection<SingleInterval>> coverageByRNA, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String rna: coverageByRNA.keySet()) {
			int size=coverageByRNA.get(rna).size();
			writer.write(rna+"\t"+size+"\n");
		}
		
		writer.close();
	}



	
	private List<String> write(MatrixWithHeaders norm, MatrixWithHeaders clusterCounts, String save, Map<String, IntervalTree<String>> genes) throws IOException {
		List<String> rtrn=new ArrayList<String>();
		FileWriter writer=new FileWriter(save);
		
		for(String rna: norm.getRowNames()) {
			Map<String, Double> scores=new TreeMap<String, Double>();
			for(String bin: norm.getColumnNames()) {
				double score=norm.get(rna,bin);
				double clusterCount=clusterCounts.get(rna, bin);
				if(score>=1 && clusterCount>=this.minCount) {scores.put(bin, score);}
			}
			
			if(scores.size()>0) {
				Collection<String> regions20=new TreeSet<String>();
				Collection<String> regions10=new TreeSet<String>();
				Collection<String> regions0=new TreeSet<String>();
				rtrn.add(rna);
				double max=max(scores);
				writer.write(rna);
				for(String bin: scores.keySet()) {
					double ratio=(scores.get(bin)/max)*100;
					if(ratio>=percentile) {regions20.add(bin);}
					//if(ratio>=10) {regions10.add(bin);}
					if(ratio>=0) {regions0.add(bin);}
				}
				
				Collection<String> overlappingGenes20=getOverlappingGenes(regions20, genes);
				//Collection<String> overlappingGenes10=getOverlappingGenes(regions10, genes);
				Collection<String> overlappingGenes0=getOverlappingGenes(regions0, genes);
				writer.write("\t"+regions20.size()+"\t"+overlappingGenes20.size()+"\t"+overlappingGenes0.size());
				for(String gene: overlappingGenes20) {
					writer.write("\t"+gene);
				}
				for(String region: regions20) {
					writer.write("\t"+region);
				}
				writer.write("\n");
			}
		}
		
		writer.close();
		return rtrn;
	}


	private Collection<String> getOverlappingGenes(Collection<String> regions, Map<String, IntervalTree<String>> genes) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String region: regions) {
			SingleInterval bin=new SingleInterval(region);
			Collection<String> set=getOverlappingGenes(bin, genes);
			rtrn.addAll(set);
		}
		
		return rtrn;
	}





	private Collection<String> getOverlappingGenes(SingleInterval bin, Map<String, IntervalTree<String>> genes) {
		Collection<String> rtrn=new TreeSet<String>();
		
		if(genes.containsKey(bin.getReferenceName())) {
			IntervalTree<String> tree=genes.get(bin.getReferenceName());
			Iterator<String> iter=tree.overlappingValueIterator(bin.getReferenceStartPosition(), bin.getReferenceEndPosition());
			while(iter.hasNext()) {
				String val=iter.next();
				rtrn.add(val);
			}
		}
		
		return rtrn;
	}





	private List<String> write(MatrixWithHeaders norm, MatrixWithHeaders clusterCounts, String save) throws IOException {
		List<String> rtrn=new ArrayList<String>();
		FileWriter writer=new FileWriter(save);
		
		for(String rna: norm.getRowNames()) {
			Map<String, Double> scores=new TreeMap<String, Double>();
			//writer.write(rna);
			for(String bin: norm.getColumnNames()) {
				double score=norm.get(rna,bin);
				double clusterCount=clusterCounts.get(rna, bin);
				if(score>=1 && clusterCount>=this.minCount) {scores.put(bin, score);}
			}
			
			if(scores.size()>0) {
				rtrn.add(rna);
				double max=max(scores);
				writer.write(rna+"\t"+scores.size());
				for(String bin: scores.keySet()) {
					double ratio=(scores.get(bin)/max)*100;
					if(ratio>=20) {
						writer.write("\t"+bin);
					}
				}
				writer.write("\n");
			}
		}
		
		writer.close();
		return rtrn;
	}





	private double max(Map<String, Double> scores) {
		double max=0;
		
		for(String s: scores.keySet()) {
			double val=scores.get(s);
			max=Math.max(max, val);
		}
		
		return max;
	}





	private MatrixWithHeaders normByMax(MatrixWithHeaders mwh) {
		Collection<String> list=new TreeSet<String>();
		for(String row: mwh.getRowNames()) {
			double max=Statistics.max(mwh.getRow(row));
			for(String col: mwh.getColumnNames()) {
				double count=mwh.get(row, col);
				if(count<minCount) {count=0;}
				else {list.add(row);}
				double norm=(count/max)*1000;
				mwh.set(row, col, norm);
			}
		}
		
		return mwh.submatrixByRowNames(list);
		
	}





	private MatrixWithHeaders normalizeByRNA(MatrixWithHeaders mwh, Map<String, Collection<SingleInterval>> coverageByRNA) {
		//double totalCounts=sum(inputCounts);
		
		Collection<String> rnaList=new TreeSet<String>();
		
		for(String row: mwh.getRowNames()) {
			double rnaTotal=coverageByRNA.get(row).size();
			//double fraction=rnaTotal/totalCounts;
			for(String col: mwh.getColumnNames()) {
				double score=mwh.get(row, col);
				double norm=(score/rnaTotal)*1000;
				//double norm=score/expected;
				if(score<minCount) {norm=0;}
				if(score>minCount) {rnaList.add(row);}
				mwh.set(row, col, norm);
			}
		}
		
		System.err.println(mwh.getNumberRows()+" "+rnaList.size());
		return mwh.submatrixByRowNames(rnaList);
		
	}





	private MatrixWithHeaders makeMatrix(Map<String, Collection<SingleInterval>> coverageByRNA, int binResolution) {
		List<String> rows=getRows(coverageByRNA.keySet());
		List<String> columns=getColumns(coverageByRNA, binResolution);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(String rna: coverageByRNA.keySet()) {
			Collection<SingleInterval> regions=coverageByRNA.get(rna);
			for(SingleInterval r: regions) {
				String genomeBin=r.bin(binResolution).toUCSC();
				rtrn.incrementCount(rna, genomeBin);
			}
		}
		
		return rtrn;
	}





	private List<String> getRows(Set<String> keySet) {
		List<String> rtrn=new ArrayList<String>();
		
		rtrn.addAll(keySet);
		
		return rtrn;
	}





	private List<String> getColumns(Map<String, Collection<SingleInterval>> coverageByRNA, int binResolution) {
		Collection<String> set=new TreeSet<String>();
		
		for(String rna: coverageByRNA.keySet()) {
			for(SingleInterval r: coverageByRNA.get(rna)) {
				String name=r.bin(binResolution).toUCSC();
				set.add(name);
			}
		}
		
		List<String> rtrn=new ArrayList<String>();
		for(String r: set) {
			rtrn.add(r);
			//System.out.println(r);
		}
		
		
		//rtrn.addAll(set);
		return rtrn;
	}





	private Collection<SingleInterval> getList(Map<String, Collection<SingleInterval>> coverageByRNA, String rna) {
		if(!coverageByRNA.containsKey(rna)) {coverageByRNA.put(rna, new TreeSet<SingleInterval>());}
		return coverageByRNA.get(rna);
	}





	private MatrixWithHeaders scoreHighResolution(MatrixWithHeaders mwh, BarcodingDataStreaming data, int binResolution) {
		List<String> possibleBins=CoordinateSpace.HG38.getBins(binResolution/10);
		
		List<String> list=new ArrayList<String>();
		for(String bin: possibleBins) {
			String binString=new SingleInterval(bin).bin(binResolution).toUCSC();
			if(mwh.containsColumn(binString)) {
				list.add(bin);
			}
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), list);
		
		while(data.hasNext()) {
			Cluster c=data.next();
			Cluster binned=c.bin(binResolution/10);
			for(String rna: c.getRNANames()) {
				for(SingleInterval bin: binned.getAllDNAIntervals()) {
					String col=bin.toUCSC();
					rtrn.incrementCount(rna, col);
				}
			}
		}
		
		data.close();
		
		return rtrn;
	}





	private void rowNorm(MatrixWithHeaders mwh) {
		for(String row: mwh.getRowNames()) {
			double[] vals=mwh.getRow(row);
			double sum=Statistics.sum(vals);
			for(String col: mwh.getColumnNames()) {
				double score=mwh.get(row, col);
				double norm=(score/sum)*1000;
				mwh.set(row, col, norm);
			}
		}
		
	}


	private Map<String, Integer> count(BarcodingDataStreaming data, Collection<String> lncRNAList) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Collection<String> rnaNames=c.getRNANames();
			for(String name: rnaNames) {
				if(lncRNAList.contains(name)) {
					int count=0;
					if(rtrn.containsKey(name)) {count=rtrn.get(name);}
					count++;
					rtrn.put(name, count);
				}
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		data.close();
		return rtrn;
	}


	private MatrixWithHeaders trimRows(MatrixWithHeaders mwh, Map<String, Integer> rnaCounts, int minCount) {
		List<String> list=new ArrayList<String>();
		for(String row: mwh.getRowNames()) {
			if(rnaCounts.containsKey(row)) {
				int count=rnaCounts.get(row);
				if(count>minCount) {list.add(row);}
			}
		}
		return mwh.submatrixByRowNames(list);
	}


	private MatrixWithHeaders normToMax(MatrixWithHeaders mwh) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		for(String row: mwh.getRowNames()) {
			double[] rowVals=mwh.getRow(row);
			double max=Statistics.max(rowVals);
			for(String col: mwh.getColumnNames()) {
				double val=mwh.get(row, col);
				double norm=(val/max)*1000;
				rtrn.set(row, col, norm);
			}
			
		}
		return rtrn;
	}


	private void removeSelfContacts(MatrixWithHeaders mwh) {
		for(String row: mwh.getRowNames()) {
			if(mwh.containsColumn(row)) {
				mwh.set(row, row, 0.0);
			}
		}
	}


	private void writeRNAMeta(String save, Map<String, SingleInterval> rnaPositions) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("RNA_Name\tchromosome\tStart\tEnd\n");
		
		for(String rna: rnaPositions.keySet()) {
			SingleInterval pos=rnaPositions.get(rna);
			writer.write(rna+"\t"+pos.getReferenceName()+"\t"+pos.getReferenceStartPosition()+"\t"+pos.getReferenceEndPosition()+"\n");
		}
		
		writer.close();
	}


	private void updatePositions(Map<String, SingleInterval> rnaPositions, Collection<RNAInterval> rnaIntervals,MatrixWithHeaders mwh) {
		for(RNAInterval interval: rnaIntervals) {
			Collection<String> allRNANames=interval.getAllRNANames(); //TODO This was attempt to fix exon/intron. Did it work?
			for(String name: allRNANames) {
				if(mwh.containsRow(name)){
					SingleInterval pos=interval.getSingleInterval();
					if(rnaPositions.containsKey(name)) {
						pos=merge(rnaPositions.get(name), pos);
					}
					rnaPositions.put(name, pos);
				}
			}
		}
		
	}


	private SingleInterval merge(SingleInterval singleInterval, SingleInterval pos) {
		SingleInterval rtrn=new SingleInterval(singleInterval.getReferenceName(), Math.min(singleInterval.getReferenceStartPosition(), pos.getReferenceStartPosition()), Math.max(singleInterval.getReferenceEndPosition(), pos.getReferenceEndPosition()));
		return rtrn;
	}


	private MatrixWithHeaders floorCounts(MatrixWithHeaders mwh, int minCount2) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()) {
			for(String col : mwh.getColumnNames()) {
				double val=mwh.get(row, col);
				if(val<minCount2) {val=0;}
				rtrn.set(row, col, val);
			}
		}
		
		return rtrn;
	}


	private MatrixWithHeaders normalizeByRNATotal(MatrixWithHeaders mwh, Map<String, Integer> rnaCounts, Map<String, Integer> geneLengths, double totalLength) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()) {
			for(String col: mwh.getColumnNames()) {
				int geneLength=geneLengths.get(col);
				
				double val=mwh.get(row, col);
				
				double num=val/(double)geneLength;
				
				double denom=get(rnaCounts,row)/totalLength;
				double norm=num/denom;
				rtrn.set(row, col, norm);
			}
		}
		
		return rtrn;
	}


	private double get(Map<String, Integer> rnaCounts, String row) {
		double rtrn=0;
		if(rnaCounts.containsKey(row)) {
			rtrn=rnaCounts.get(row);
		}
		
		return rtrn;
	}


	private double totalLength(Map<String, Integer> geneLengths) {
		double sum=0;
		
		for(String gene: geneLengths.keySet()) {
			sum+=geneLengths.get(gene);
		}
		
		return sum;
	}


	private MatrixWithHeaders normalizeByGeneLength(MatrixWithHeaders mwh, Map<String, Integer> geneLengths) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()) {
			for(String col: mwh.getColumnNames()) {
				double val=mwh.get(row, col);
				int length=geneLengths.get(col);
				double norm=1000000*(val/(double)length);
				rtrn.set(row, col, norm);
			}
		}
		
		return rtrn;
	}


	private Map<String, Integer> getGeneLengths(Map<String, IntervalTree<String>> genes) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(String chr: genes.keySet()) {
			IntervalTree<String> tree=genes.get(chr);
			Iterator<Node<String>> iter=tree.iterator();
			while(iter.hasNext()) {
				Node<String> val=iter.next();
				String name=val.getValue();
				int length=val.getLength();
				rtrn.put(name, length);
			}
		}
		
		
		return rtrn;
	}


	private List<String> getGeneNames(Map<String, Integer> geneLengths) {
		List<String> rtrn=new ArrayList<String>();
		
		rtrn.addAll(geneLengths.keySet());
		
		return rtrn;
	}


	private List<String> getRNAs(Map<String, Integer> rnaCounts) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String name: rnaCounts.keySet()) {
			int count=rnaCounts.get(name);
			if(count>minCount) {rtrn.add(name);}
		}
		
		return rtrn;
	}


	private Map<String, Integer> getRNACounts(BarcodingDataStreaming data) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Collection<String> names=c.getRNANames(); //TODO need to make sure this handles exons/introns well
			for(String name: names) {
				int count=0;
				if(rtrn.containsKey(name)) {count=rtrn.get(name);}
				count++;
				rtrn.put(name, count);
			}
			
			
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		data.close();
		
		return rtrn;
	}


	private void writeMeta(String save, Map<String, IntervalTree<String>> genes) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Name\tchromosome\tstart\tend\tlength\n");
		
		for(String chr: genes.keySet()) {
			Iterator<Node<String>> iter=genes.get(chr).iterator();
			while(iter.hasNext()) {
				Node<String> gene=iter.next();
				writer.write(gene.getValue()+"\t"+chr+"\t"+gene.getStart()+"\t"+gene.getEnd()+"\t"+gene.getLength()+"\n");
			}
		}
		
		writer.close();
	}


	private MatrixWithHeaders trimMatrix(MatrixWithHeaders mwh) {
		mwh=trimRows(mwh);
		mwh=trimColumns(mwh);
		return mwh;
	}


	private MatrixWithHeaders trimColumns(MatrixWithHeaders mwh) {
		List<String> list=new ArrayList<String>();
		for(String col: mwh.getColumnNames()) {
			double[] vals=mwh.getColumn(col);
			if(Statistics.max(vals)>0) {list.add(col);}
		}
		
		System.err.println(mwh.getColumnNames().size()+" "+list.size());
		
		mwh=mwh.submatrixByColumnNames(list);
		
		return mwh;
	}


	private MatrixWithHeaders trimRows(MatrixWithHeaders mwh) {
		List<String> list=new ArrayList<String>();
		for(String row: mwh.getRowNames()) {
			double[] vals=mwh.getRow(row);
			if(Statistics.max(vals)>0) {list.add(row);}
		}
		
		System.err.println(mwh.getRowNames().size()+" "+list.size());
		
		mwh=mwh.submatrixByRowNames(list);
		
		return mwh;
	}


	private Collection<String> overlappingGenes(Cluster c, Map<String, IntervalTree<String>> genes) {
		Collection<String> rtrn=new TreeSet<String>();
		
		Collection<SingleInterval> regions=c.getAllDNAIntervals();
		
		for(SingleInterval r: regions) {
			if(genes.containsKey(r.getReferenceName())) {
				Iterator<String> iter=genes.get(r.getReferenceName()).overlappingValueIterator(r.getReferenceStartPosition(), r.getReferenceEndPosition());
				while(iter.hasNext()) {
					rtrn.add(iter.next());
				}
			}
		}
		
		return rtrn;
	}


	private void iterateMatrix(Collection<String> overlappingGenes, Collection<String> rnaNames, MatrixWithHeaders mwh, Map<String, Integer> rnaCounts) {
		for(String geneName: overlappingGenes) {
			for(String rnaName: rnaNames) {
				if(mwh.hasColumn(geneName)&& mwh.hasRow(rnaName)) {
					mwh.incrementCount(rnaName, geneName);
				}
				
				int count=0;
				if(rnaCounts.containsKey(rnaName)) {
					count=rnaCounts.get(rnaName);
				}
				count++;
				rnaCounts.put(rnaName, count);
				
			}
		}
	}


	private static List<String> parseList(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		
		List<String> rtrn=new ArrayList<String>();
		for(String line: lines) {
			rtrn.add(line.split("\t")[0]);
		}
		
		return rtrn;
	}
	
	
	private static Map<String, Double> convert(TreeMap<SingleInterval, Double> loadbedgraph) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(SingleInterval r: loadbedgraph.keySet()) {
			double score=loadbedgraph.get(r);
			rtrn.put(r.toUCSC(), score);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>6) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			int binResolution=Integer.parseInt(args[1]);
			String save=args[2];
			List<String> lncRNAList=parseList(args[3]);
			//Map<String, Double> map=convert(BEDFileIO.loadbedgraph(new File(args[4])));
			System.err.println(lncRNAList.size());
			int minCount=Integer.parseInt(args[4]);
			
			Map<String, IntervalTree<String>> genes=BEDFileIO.loadRefFlatTree(args[5]);
			int percentile=Integer.parseInt(args[6]);
			
			new lncRNAsByBinMatrix(data, binResolution, save, lncRNAList, minCount, genes, percentile);
		}
		else {System.err.println(usage);}
	}
	
	

	



	static String usage=" args[0]=barcoding data \n args[1]=bin resolution \n args[2]=save \n args[3]=lncRNA list \n args[4]=min number of clusters \n args[5]=ref flat file \n args[6]=percentile to filter";
}

package guttmanlab.core.rnasprite;

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
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class lncRNAsByGeneMatrix {

	int minCount=2;
	int geneCount=10;
	
	public lncRNAsByGeneMatrix(BarcodingDataStreaming data, Map<String, IntervalTree<String>> genes, String save, Collection<String> lncRNAList) throws IOException {
		
		//TODO Mask self interactions
		//TODO split by exon/Intron
		
		
		
		Map<String, Integer> geneLengths=getGeneLengths(genes);
		List<String> geneNames=getGeneNames(geneLengths);
		
		double totalLength=totalLength(geneLengths);
		
		Map<String, Integer> lncRNACounts=count(data, lncRNAList);
		System.err.println(lncRNACounts.size());
		List<String> rnas=new ArrayList<String>();
		for(String name: lncRNACounts.keySet()) {
			int count=lncRNACounts.get(name);
			System.out.println(name+"\t"+count);
			if(count>geneCount) {rnas.add(name);}
		}
		System.err.println(rnas.size());
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rnas, geneNames);
		
		Map<String, Integer> rnaCounts=new TreeMap<String, Integer>();
		Map<String, SingleInterval> rnaPositions=new TreeMap<String, SingleInterval>();
		
		//iterate over all clusters
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Collection<String> overlappingGenes=overlappingGenes(c, genes);
			Collection<String> rnaNames=c.getRNANames();
			
			Collection<RNAInterval> rnaIntervals=c.getAllRNARegions();
			
			iterateMatrix(overlappingGenes, rnaNames, mwh, rnaCounts);
			updatePositions(rnaPositions, rnaIntervals, mwh);
			
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
			
		}
		data.close();
		
		System.err.println(mwh.getRowNames().size()+" "+mwh.getColumnNames().size());
		//mwh.write(save);
		
		System.err.println("Flooring counts");
		mwh=floorCounts(mwh, minCount);
		
		System.err.println("Trimming matrix");
		mwh=trimMatrix(mwh);
		System.err.println(mwh.getRowNames().size()+" "+mwh.getColumnNames().size());
		mwh.write(save+".counts.matrix");
		
		
		
		System.err.println("Writing and normalizing");
		
		// Normalize genes by length
		MatrixWithHeaders lengthNorm=normalizeByGeneLength(mwh, geneLengths);
		lengthNorm.write(save+".geneLengthNorm.matrix");
		System.err.println("Length normalized");
		
		MatrixWithHeaders normToMax=normToMax(lengthNorm);
		normToMax.write(save+".normToMax.matrix");
		System.err.println("max normalized");
		
		//TODO Normalize to total in the row
		MatrixWithHeaders rnaNorm=normalizeByRNATotal(mwh, rnaCounts, geneLengths, totalLength);
		rnaNorm.write(save+".rnaTotalNorm.matrix");
		System.err.println("rna total normalized");
		
		writeMeta(save+".gene.metadata", genes);
		writeRNAMeta(save+".RNA.metadata", rnaPositions);
		
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


	private static Collection<String> parseList(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		
		Collection<String> rtrn=new TreeSet<String>();
		for(String line: lines) {
			rtrn.add(line.split("\t")[0]);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>3) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Map<String, IntervalTree<String>> genes=BEDFileIO.loadRefFlatTree(args[1]);
			String save=args[2];
			Collection<String> lncRNAList=parseList(args[3]);
			System.err.println(lncRNAList.size());
			new lncRNAsByGeneMatrix(data, genes, save, lncRNAList);
		}
		else {System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=barcoding data \n args[1]=genes (BED file) \n args[2]=save \n args[3]=lncRNA list";
}

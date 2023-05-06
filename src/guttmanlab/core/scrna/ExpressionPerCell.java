package guttmanlab.core.scrna;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class ExpressionPerCell {

	public ExpressionPerCell(File[] bams, String save, Map<String, IntervalTree<String>> geneNames) throws IOException {
		
		List<String> colList=new ArrayList<String>();
		
		List<String> rowList=getRows(bams);
		rowList.add("total");
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rowList, tolist(bams));
		
		
		for(int i=0; i<bams.length; i++) {
			System.err.println(bams[i].getName());
			int count=count(bams[i]);
			if(count>1000) {
				Map<SingleInterval, Double> counts=quantify(bams[i]);
				add(counts, mwh, bams[i].getName());
				colList.add(bams[i].getName());
				mwh.set("total", bams[i].getName(), count);
			}
		}
		
		mwh=compress(mwh);
		mwh=mwh.submatrixByColumnNames(colList);
		
		Map<String, String> pidToName=getNames(mwh, geneNames);
		
		mwh.setPIDToName(pidToName);
		mwh.write(save);
		
		MatrixWithHeaders MWHbyGene=collapseByGene(mwh);
		MWHbyGene.write(save+".gene.matrix");
	}

	private MatrixWithHeaders collapseByGene(MatrixWithHeaders mwh) {
		List<String> rows=new ArrayList<String>();
		for(String row: mwh.getRowNames()) {
			if(mwh.getPIDToName().containsKey(row)) {
				String name=mwh.getPIDToName().get(row);
				if(!rows.contains(name)) {rows.add(name);}
			}
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()) {
			if(mwh.getPIDToName().containsKey(row)) {
				String name=mwh.getPIDToName().get(row);
				for(String col: mwh.getColumnNames()) {
					double current=rtrn.get(name, col);
					double add=mwh.get(row, col);
					double sum=current+add;
					rtrn.set(name, col, sum);
				}
			}
		}
		
		for(String row: rtrn.getRowNames()) {
			for(String col: rtrn.getColumnNames()) {
				double total=mwh.get("total", col);
				double sum=rtrn.get(row, col);
				double TPM=(sum/total)*1000000;
				rtrn.set(row, col, TPM);
			}
		}
		
		
		return rtrn;
	}

	private Map<String, String> getNames(MatrixWithHeaders mwh, Map<String, IntervalTree<String>> geneNames) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(String row: mwh.getRowNames()) {
			if(!row.equals("total")) {
				SingleInterval interval=new SingleInterval(row);
				if(geneNames.containsKey(interval.getReferenceName())) {
					Iterator<String> iter=geneNames.get(interval.getReferenceName()).overlappingValueIterator(interval.getReferenceStartPosition(), interval.getReferenceEndPosition());
					if(iter.hasNext()) {
						String name=iter.next();
						rtrn.put(row, name);
					}
				}
			}
		}
		return rtrn;
	}

	private List<String> getRows(File[] bams) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(int i=0; i<bams.length; i++) {
			Map<SingleInterval, Double> counts=quantify(bams[i]);
			rtrn.addAll(counts.keySet());
		}
		
		ArrayList<String> list=new ArrayList<String>();
		for(SingleInterval g: rtrn) {
			list.add(g.getSingleInterval().toUCSC(g.getOrientation()));
		}
		
		return list;
	}

	private MatrixWithHeaders compress(MatrixWithHeaders mwh) {
		List<String> list=new ArrayList<String>();
		for(String row: mwh.getRowNames()) {
			double sum=Statistics.sum(mwh.getRow(row));
			if(sum>0) {list.add(row);}
		}
		return mwh.submatrixByRowNames(list);
	}

	private List<String> tolist(File[] bams) {
		List<String> rtrn=new ArrayList<String>();
	
		for(File bam: bams) {
			if(count(bam)>1000) {
				rtrn.add(bam.getName());
			}
		}
		
		return rtrn;
	}

	
	private void add(Map<SingleInterval, Double> counts, MatrixWithHeaders mwh, String name) {
		for(SingleInterval g: counts.keySet()) {
			String row=g.getSingleInterval().toUCSC(g.getOrientation());
			mwh.set(row, name, counts.get(g));
		}
		
	}

	
	private int  count(File file) {
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			counter++;
		}
		
		reader.close();
		reads.close();
		
		return counter;
	}
	
	private Map<SingleInterval, Double> quantify(File file) {
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment interval=new SAMFragment(record);
			
			/*SingleInterval interval=new SingleInterval(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentStart()+1);
			interval.setOrientation(record.getreads);*/
			
			double count=0;
			if(rtrn.containsKey(interval.getSingleInterval())) {count=rtrn.get(interval.getSingleInterval());}
			count++;
			rtrn.put(interval.getSingleInterval(), count);
				
			
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		System.err.println(file.getName()+" "+counter);
		reader.close();
		reads.close();
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		File[] bams=new File(args[0]).listFiles();
		String save=args[1];
		Map<String, IntervalTree<String>> geneNames=BEDFileIO.loadGeneNamesFromRefFlat(args[2]);
		new ExpressionPerCell(bams, save, geneNames);
	}
}

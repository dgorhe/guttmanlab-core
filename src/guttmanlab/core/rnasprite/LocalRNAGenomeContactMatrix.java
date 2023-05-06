package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import net.sf.samtools.util.CloseableIterator;

public class LocalRNAGenomeContactMatrix {

	int windowSize=10;
	
	public LocalRNAGenomeContactMatrix(BarcodingDataStreaming data, String gene1, String save) throws IOException{
		
		List<String> positions=getAllPositions(data);
		Collection<Cluster> clusters=data.getRNAClusters(gene1);
		
		System.err.println(clusters.size());
		
		Map<String, SingleInterval> rnaRegions=getRegions(clusters);
		
		
		AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows1=getRegions(rnaRegions, gene1, windowSize);
		
		List<String> rowNames=getNames(windows1);
		List<String> columnNames=positions;
		columnNames.add("Input");
		
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rowNames, columnNames);
		
		for(Cluster c: clusters){
			update(windows1, rnaRegions.get(gene1), c, mwh);
		
		}
		
		mwh.write(save);
	}
	
	

	private List<String> getAllPositions(BarcodingDataStreaming data) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		while(data.hasNext()){
			Cluster binned=data.next().bin(1000000);
			for(SingleInterval region: binned.getAllDNAIntervals()){rtrn.add(region);}
		}
		data.close();
		
		ArrayList<String> list=new ArrayList<String>();
		
		for(SingleInterval region: rtrn){
			list.add(region.toUCSC());
		}
		return list;
	}



	private void update(AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows1, SingleInterval region, Cluster c, MatrixWithHeaders mwh) {
		//CloseableIterator<DerivedAnnotation<? extends Annotation>> iter=windows1.sortedIterator(region, false);
		Cluster binned=c.bin(1000000);
		
		for(SingleInterval rnaRegion: c.getAllRNARegions()){
			CloseableIterator<DerivedAnnotation<? extends Annotation>> iter=windows1.sortedIterator(rnaRegion, false);
			while(iter.hasNext()){
				String row=iter.next().toUCSC();
				double total=mwh.get(row, "Input");
				total++;
				mwh.set(row, "Input", total);
				for(SingleInterval genome: binned.getAllDNAIntervals()){
					String column=genome.toUCSC();
					double score=mwh.get(row, column);
					score++;
					mwh.set(row, column, score);
				}
			}
		}
		
		
		
	}



	private Map<String, SingleInterval> getRegions(Collection<Cluster> clusters) {
		Map<String, SingleInterval> rtrn=new TreeMap<String, SingleInterval>();
		
		for(Cluster c: clusters){
			for(SingleInterval rna: c.getAllRNARegions()){
				String name=rna.getName();
				SingleInterval other=null;
				if(rtrn.containsKey(name)){
					other=rtrn.get(name);
				}
				other=merge(other, rna);
				rtrn.put(name, other);
			}
		}
		
		return rtrn;
	}



	private SingleInterval merge(SingleInterval other, SingleInterval rna) {
		if(other==null){return rna;}
		return update(other, rna);
	}
	
	private SingleInterval update(SingleInterval rna1, SingleInterval rna2) {
		String chr=rna1.getReferenceName();
		int start=Math.min(rna1.getReferenceStartPosition(), rna2.getReferenceStartPosition());
		int end=Math.max(rna1.getReferenceEndPosition(), rna2.getReferenceEndPosition());
		SingleInterval rtrn=new SingleInterval(chr, start, end);
		rtrn.setName(rna1.getName());
		return rtrn;
	}



	private void update(MatrixWithHeaders mwh, AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows1,
			AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows2, SingleInterval region1,
			SingleInterval region2, int clusterSize) {
		
		CloseableIterator<DerivedAnnotation<? extends Annotation>> iter1=windows1.sortedIterator(region1, false);
		
		while(iter1.hasNext()){
			String row=iter1.next().toUCSC();
			CloseableIterator<DerivedAnnotation<? extends Annotation>> iter2=windows2.sortedIterator(region2, false);
			while(iter2.hasNext()){
				String column=iter2.next().toUCSC();
				double score=mwh.get(row, column);
				//score+=(2.0/(double)clusterSize);
				score++;
				mwh.set(row, column, score);
			}
			iter2.close();
		}
		
		iter1.close();
		
	}



	



	private List<String> getNames(AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows) {
		List<String> rtrn=new ArrayList<String>();
		
		CloseableIterator<DerivedAnnotation<? extends Annotation>> iter=windows.sortedIterator();
		while(iter.hasNext()){
			String line=iter.next().toUCSC();
			rtrn.add(line);
		}
		
		return rtrn;
	}



	private AnnotationCollection<DerivedAnnotation<? extends Annotation>> getRegions(Map<String, SingleInterval> rnaRegions, String gene1, int windowSize) {
		SingleInterval region=rnaRegions.get(gene1);
		System.err.println(gene1+" "+region.toUCSC());
		AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows=region.getWindows(windowSize, windowSize);
		return windows;
	}

	public static void main(String[] args) throws IOException{
		if(args.length>2){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String gene1=args[1];
			String save=args[2];
			new LocalRNAGenomeContactMatrix(data, gene1, save);
		} 
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=gene 1 \n args[2]=save";
	
}

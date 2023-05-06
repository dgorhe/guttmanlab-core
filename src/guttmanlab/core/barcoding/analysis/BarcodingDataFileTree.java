package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;

public class BarcodingDataFileTree implements SPRITEData{

	Map<String, IntervalTree<File>> barcodeFileTree;
	Collection<File> files;
	Iterator<File> fileIterator;
	BarcodingDataStreaming currentData;
	Collection<String> visited;
	
	public BarcodingDataFileTree(File[] files) throws IOException{
		this.barcodeFileTree=makeFileTree(files);
		this.files=new ArrayList<File>();
		this.visited=new TreeSet<String>();
		for(int i=0; i<files.length; i++){this.files.add(files[i]);}
		reset();
	}
	
	private void reset() throws IOException {
		//Reset file iterator
		this.fileIterator=files.iterator();
		
		File nextFile=fileIterator.next();
		//Make current data
		this.currentData=new BarcodingDataStreaming(nextFile);
		
		//reset visited
		this.visited=new TreeSet<String>();
		
	}

	private SingleInterval parse(String name) {
		String chr=name.split("_")[0];
		int start=new Integer(name.split("_")[1]);
		int end=new Integer(name.split("_")[2]);
		return new SingleInterval(chr, start, end);
	}
	
	private Map<String, IntervalTree<File>> makeFileTree(File[] files) {
		Map<String, IntervalTree<File>> rtrn=new TreeMap<String, IntervalTree<File>>();
		for(int i=0; i<files.length; i++){
			String name=files[i].getName();
			SingleInterval region=parse(name);
			IntervalTree<File> tree=new IntervalTree<File>();
			if(rtrn.containsKey(region.getReferenceName())){
				tree=rtrn.get(region.getReferenceName());
			}
			tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), files[i]);
			rtrn.put(region.getReferenceName(), tree);
		}
		return rtrn;
	}

	@Override
	public void close() {
		currentData.close();
		try {
			reset();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public boolean hasNext() {
		try{
			if(!currentData.hasNext()){
				if(this.fileIterator.hasNext()){
					currentData.close();
					File nextFile=fileIterator.next();
					//System.err.println(nextFile);
					currentData=new BarcodingDataStreaming(nextFile);
					return hasNext();
				}
				else{return false;}
			}
			return true;
		} catch (IOException e) {
			throw new IllegalStateException(e.toString());
		}
	}

	@Override
	public Cluster next() {
		return currentData.next();
	}

	@Override
	public Collection<Cluster> getClustersOverlappingRegion(Annotation region) throws IOException {
		return getClustersOverlappingRegion(region, 1, 1000);
	}
	
	public Collection<Cluster> getClustersOverlappingRegion(Collection<Cluster> clusters, Annotation region) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(Cluster c: clusters){
			if(c.overlapsInterval(region)){rtrn.add(c);}
		}
		
		return rtrn;
	}

	public Collection<Cluster> getClustersOverlappingRegion(Annotation region, int min, int max) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		Iterator<File> barcodingData=this.barcodeFileTree.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
		
		while(barcodingData.hasNext()){
			File file=barcodingData.next();
			System.err.println("Max cluster size "+max);
			Collection<Cluster> clusters=new BarcodingDataStreaming(file).getClustersOverlappingRegion(region, min, max);
			rtrn.addAll(clusters);
		}
		
		return rtrn;
	}

	public Collection<Cluster> getClustersOverlappingMultipleRegions(Collection<? extends Annotation> regions) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		boolean started=false;
		
		for(Annotation region: regions){
			if(!started){
				rtrn=getClustersOverlappingRegion(region);
				started=true;
			}
			else{
				rtrn=getClustersOverlappingRegion(rtrn,region);
			}
			
		}
		return rtrn;
	}

	
	
}

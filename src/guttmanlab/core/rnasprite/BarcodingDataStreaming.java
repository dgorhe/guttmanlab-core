package guttmanlab.core.rnasprite;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;
import htsjdk.samtools.util.CloseableIterator;

public class BarcodingDataStreaming implements CloseableIterator<Cluster>{

	private File barcodeFile;
	private BufferedReader reader;
	private String nextLine;
	private Cluster nextCluster;
	Collection<ClusterFilter> filters;
	
	public BarcodingDataStreaming(File barcodeFile) throws IOException{
		this.barcodeFile=barcodeFile;
		this.filters=new ArrayList<ClusterFilter>();
		this.reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
	}
	

	
	@Override
	public boolean hasNext() {
		
		try {
			nextLine = reader.readLine();
				
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			
		}
		return nextLine!=null;
	}


	@Override
	public Cluster next() {
		return Cluster.parseCluster(nextLine);
	}

	@Override
	public void close() {
		try {
			reader.close();
			this.reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}
	
	public void reset(){
		close();
	}

	

	

	public Map<SingleInterval, Double> getDNAContactFrequency(SingleInterval position, int binResolution) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		int counter=0;
		int retained=0;
		while(hasNext()){
			Cluster c=next();
			Cluster binned=c.bin(binResolution);
			if(binned.containsOverlappingDNA(position)){
				retained++;
				Collection<SingleInterval> dnaRegions=binned.getAllDNAIntervals();
				for(SingleInterval region: dnaRegions){
					if(!region.overlaps(position)){
						double count=0;
						if(rtrn.containsKey(region)){
							count=rtrn.get(region);
						}
						count+=((2.0/(double)c.getAllDNAIntervals().size()));
						rtrn.put(region, count);
					}
				}
			}
			counter++;
			if(counter%100000==0){System.err.println(counter+" "+retained);}
		}
		
		reset();
		return rtrn;
	}
	
	public Pair<MatrixWithHeaders> getRNADNAContactMatrix(int binResolution){
		return getRNADNAContactMatrix(0, binResolution);
	}
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, SingleInterval region, String gene, boolean weight){
		List<String> regions=getGenomePositions(binResolution, region);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, gene, weight);
		return counts;
	}
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, SingleInterval region, Kmer genes, boolean weight){
		List<String> regions=getGenomePositions(binResolution, region);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, genes, weight);
		return counts;
	}
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, String chr, Kmer genes, boolean weight){
		List<String> regions=getGenomePositions(binResolution, chr);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, genes, weight);
		return counts;
	}
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, Kmer genes, boolean weight){
		List<String> regions=getGenomePositions(binResolution);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, genes, weight);
		return counts;
	}
	
	
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, boolean weight, Collection<Cluster> clusters){
		List<String> regions=getGenomePositions(binResolution);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, weight, clusters);
		return counts;
	}
	
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, boolean weight, Collection<Cluster> clusters, SingleInterval region){
		List<String> regions=getGenomePositions(binResolution, region);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, weight, clusters);
		return counts;
	}
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, boolean weight){
		List<String> regions=getGenomePositions(binResolution);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, weight);
		return counts;
	}
	
	public MatrixWithHeaders getDNADNAContactMatrix(SingleInterval region, int binResolution, boolean weight){
		List<String> regions=getGenomePositions(binResolution, region);
		
		System.err.println(region.toUCSC()+" "+regions.size());
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, weight);
		return counts;
	}
	
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, boolean weight, Map<String, Integer> chrSizes){
		List<String> regions=getGenomePositions(chrSizes, binResolution);
		System.err.println(regions.size());
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, weight);
		return counts;
	}
	
	public MatrixWithHeaders getDNADNAContactMatrix(String chr, int binResolution, boolean weight, int size){
		List<String> regions=getGenomePositions(chr, size, binResolution);
		System.err.println(regions.size());
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, weight);
		return counts;
	}
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, boolean weight, String chr, int size){
		List<String> regions=getGenomePositions(chr,size, binResolution);
		System.err.println(regions.size());
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, weight);
		return counts;
	}
	
	
	public MatrixWithHeaders getDNADNAContactMatrixByChromosome(boolean weight, List<String> chromosomes){
		MatrixWithHeaders counts=new MatrixWithHeaders(chromosomes, chromosomes);
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			c=c.bin(1000000);
			for(SingleInterval region1: c.getAllDNAIntervals()){
				for(SingleInterval region2: c.getAllDNAIntervals()){
					if(!region1.equals(region2)) {
						String row=region1.getReferenceName();
						String column=region2.getReferenceName();
						if(counts.containsColumn(column) && counts.containsRow(row)) {
							double count=counts.get(row, column);
							double score=2.0/c.getAllDNAIntervals().size();
							if(weight){
								count+=score;
							}
							else{
								count++;
							}
							counts.set(row, column, count);
						}
					}
					}
				}
			counter++;
			if(counter%10000==0) {System.err.println(counter);}
			
			}
		close();
		return counts;
	}
	
	
	public MatrixWithHeaders getRNARNAContactMatrix(int binResolution, boolean weight, CoordinateSpace space){
		List<String> regions=space.getBins(binResolution);
		
		System.err.println(regions.size());
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreRNARNA(counts, binResolution, weight);
		
		
		
		return counts;
	}
	
	
	public MatrixWithHeaders getRNARNAContactMatrix(List<String> list, boolean weight){
		Map<String, Integer> scores=getAllRNAList();
		List<String> regions=getList(list, scores, 500);
		
		System.err.println(regions.size());
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreRNARNA(counts, weight);
		
		counts=normalize(counts, scores, 10);
		
		return counts;
	}
	

	private MatrixWithHeaders normalize(MatrixWithHeaders counts, Map<String, Integer> scores, int minNumInstances) {
		double total=scores.get("total");
		double minFold=5.0;
		Collection<String> include=new TreeSet<String>();
		
		System.err.println(counts.rowDimension()+" "+counts.columnDimension());
		
		
		for(String row: counts.getRowNames()) {
			double eRow=scores.get(row)/total;
			for(String column: counts.getColumnNames()) {
				double eColumn=scores.get(column)/total;
				double o=counts.get(row, column);
				double e=eRow*eColumn*total;
				double norm=o/e;
				if(o<minNumInstances) {norm=0;}
				if(norm>minFold) {include.add(row);}
				counts.set(row, column, norm);
			}
		}
		
		counts=counts.submatrixByColumnNames(include);
		counts=counts.submatrixByRowNames(include);
		
		System.err.println(counts.rowDimension()+" "+counts.columnDimension());
		return counts;
	}


	private MatrixWithHeaders makeExpected(MatrixWithHeaders counts, Map<String, Integer> scores) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(counts.getRowNames(), counts.getColumnNames());
		double total=scores.get("total");
		
		for(String row: counts.getRowNames()) {
			double eRow=scores.get(row)/total;
			for(String column: counts.getColumnNames()) {
				double eColumn=scores.get(column)/total;
				double expected=eRow*eColumn*total;
				rtrn.set(row, column, expected);
			}
		}
		
		return rtrn;
	}



	private List<String> getList(List<String> list, Map<String, Integer> map, int minNum) {
		List<String> rtrn=new ArrayList<String>();
		for(String name: list) {
			if(map.containsKey(name)) {
				int count=map.get(name);
				System.out.println(name+"\t"+count);
				if(count>minNum) {rtrn.add(name);}
			}
		}
		
		
		return rtrn;
	}



	private void scoreRNARNA(MatrixWithHeaders counts, boolean weight) {
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			Collection<String> names=new TreeSet<String>();
			for(RNAInterval r: c.getAllRNARegions()){
				names.add(r.getName());
			}
			
			for(String name1: names){
				for(String name2: names){
					if(counts.containsRow(name1) && counts.containsColumn(name2) && !name1.equals(name2)){
						counts.incrementCount(name1, name2);
					}
				}
			}
			counter++;
			if(counter %1000000 ==0){System.err.println(counter);}
		}
		close();
		
	}
	
	
	private void scoreRNARNA(MatrixWithHeaders counts, int binResolution, boolean weight) {
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			Collection<String> intronList=new TreeSet<String>();
			
			for(RNAInterval r: c.getAllRNARegions()){
				if(r.isIntron()){
					SingleInterval region=r.bin(binResolution);
					intronList.add(region.toUCSC());
				}
			}
			
			for(String name1: intronList){
				for(String name2: intronList){
					if(counts.containsRow(name1) && counts.containsColumn(name2) && !name1.equals(name2)){
						counts.incrementCount(name1, name2);
					}
				}
			}
			counter++;
			if(counter %1000000 ==0){System.err.println(counter);}
		}
		close();
		
	}



	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, Kmer genes, boolean weight, Map<String, Kmer> collapseSets){
		List<String> regions=getGenomePositions(binResolution);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, genes, weight, collapseSets);
		return counts;
	}
	
	
	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, SingleInterval region, boolean weight){
		List<String> regions=getGenomePositions(binResolution, region);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, weight);
		return counts;
	}
	
	
	private List<String> getGenomePositions(int binResolution, SingleInterval region) {
		ArrayList<String> list=new ArrayList<String>();
		
		SingleInterval start=region.bin(binResolution);
		SingleInterval end=new SingleInterval(region.getReferenceName(), region.getReferenceEndPosition()-1, region.getReferenceEndPosition()).bin(binResolution);
		
		for(int i=start.getReferenceStartPosition(); i<end.getReferenceEndPosition(); i+=binResolution){
			SingleInterval bin=new SingleInterval(start.getReferenceName(), i, i+binResolution);
			//System.err.println(bin.toUCSC());
			list.add(bin.toUCSC());
		}
		
		return list;
	}



	public MatrixWithHeaders getDNADNAContactMatrix(int binResolution, String chr, String gene){
		List<String> regions=getGenomePositions(binResolution, chr);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		scoreDNADNA(counts, binResolution, gene, false);
		return counts;
	}
	
	private List<String> getGenomePositions(int binResolution, String chr) {
		Collection<SingleInterval> temp=new TreeSet<SingleInterval>();
		while(hasNext()){
			Cluster c=next();
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region:binned.getAllDNAIntervals()){
				if(region.getReferenceName().equals(chr)){
					temp.add(region);
				}
			}
		}
		close();
		
		List<String> rtrn=new ArrayList<String>();
		for(SingleInterval region: temp){rtrn.add(region.toUCSC());}
		return rtrn;
	}
	
	public List<String> getGenomePositions(int binResolution) {
		Collection<SingleInterval> temp=new TreeSet<SingleInterval>();
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()<1000) {
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region:binned.getAllDNAIntervals()){
				temp.add(region);
			}
			}
		}
		close();
		
		List<String> rtrn=new ArrayList<String>();
		for(SingleInterval region: temp){rtrn.add(region.toUCSC());}
		return rtrn;
	}
	
	
	public List<String> getGenomePositions(Map<String, Integer> sizes, int binResolution) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String chr: sizes.keySet()) {
			for(int i=0; i<sizes.get(chr); i+=binResolution) {
				String name=chr+":"+i+"-"+(i+binResolution);
				rtrn.add(name);
			}
		}
		
		return rtrn;
	}
	
	public List<String> getGenomePositions(String chr, int size, int binResolution) {
		List<String> rtrn=new ArrayList<String>();
		
		
			for(int i=0; i<size; i+=binResolution) {
				String name=chr+":"+i+"-"+(i+binResolution);
				rtrn.add(name);
			}
		
		
		return rtrn;
	}
	
	
	
	public Collection<Cluster> getClusters(Kmer genes, Map<String, Kmer> collapseSets, int n) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(genes, collapseSets, n)){
				rtrn.add(c);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter+" "+rtrn.size());}
		}
		close();
		return rtrn;
	}
	
	public Collection<Cluster> getClusters(Kmer genes, Map<String, Kmer> collapseSets) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(genes, true, collapseSets)){
				rtrn.add(c);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter+" "+rtrn.size());}
		}
		close();
		return rtrn;
	}
	
	
	
	public static Collection<Cluster> getClusters(Collection<Cluster> clusters, Kmer genes, boolean requireAll) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		int counter=0;
		for(Cluster c: clusters){
			if(c.containsRNA(genes, requireAll)){
				rtrn.add(c);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter+" "+rtrn.size());}
		}
		return rtrn;
	}
	
	public Collection<Cluster> getClusters(Kmer genes, boolean requireAll) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(genes, requireAll)){
				rtrn.add(c);
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter+" "+rtrn.size());}
		}
		close();
		return rtrn;
	}
	
	
	public Collection<Cluster> getClusters(Kmer genes) {
		return getClusters(genes, true);
	}
	
	public Collection<Cluster> getDNAClusters(Kmer body, int binResolution) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			Cluster binned=c.bin(binResolution);
			if(containsDNA(binned, body)){
				rtrn.add(binned);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter+" "+rtrn.size());}
		}
		close();
		return rtrn;
	}
	
	
	public void getOriginalDNAClusters(Kmer body, int binResolution, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		int counter=0;
		int retained=0;
		while(hasNext()){
			Cluster c=next();
			Cluster binned=c.bin(binResolution);
			if(containsDNA(binned, body)){
				writer.write(c.toString()+"\n");
				retained++;
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter+" "+retained);}
		}
		close();
		writer.close();
	}
	
	private boolean containsDNA(Cluster c, Kmer body) {
		for(SingleInterval region: body.getIntervals()) {
			if(c.getAllDNAIntervals().contains(region)) {return true;}
		}
		return false;
	}



	/*private void scoreDNADNA(MatrixWithHeaders counts, int binResolution, boolean weight) {
		while(hasNext()){
			Cluster c=next();
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					String row=region1.toUCSC();
					String column=region2.toUCSC();
					if(!row.equals(column) && counts.containsColumn(column) && counts.containsRow(row)){
						double count=counts.get(row, column);
						double score=2.0/c.getAllDNAIntervals().size();
						//count++;
						if(weight){count+=score;}
						else{count++;}
						//System.err.println(count);
						counts.set(row, column, count);
					}
				}
			
		}
		}
		close();
	}*/
	
	private void scoreDNADNA(MatrixWithHeaders counts, int binResolution, boolean weight, Collection<Cluster> clusters) {
		for(Cluster c: clusters){
			
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					String row=region1.toUCSC();
					String column=region2.toUCSC();
					if(!row.equals(column) && counts.containsColumn(column) && counts.containsRow(row)){
						double count=counts.get(row, column);
						double score=2.0/c.getAllDNAIntervals().size();
						//count++;
						if(weight){count+=score;}
						else{count++;}
						//System.err.println(count);
						counts.set(row, column, count);
					}
				}
			
		}
		}
	}

	private void scoreDNADNA(MatrixWithHeaders counts, int binResolution, Kmer genes, boolean weight, Map<String, Kmer> collapseSets) {
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(genes, true, collapseSets)){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					String row=region1.toUCSC();
					String column=region2.toUCSC();
					if(!row.equals(column) && counts.containsColumn(column) && counts.containsRow(row)){
						double count=counts.get(row, column);
						double score=2.0/c.getAllDNAIntervals().size();
						//count++;
						if(weight){count+=score;}
						else{count++;}
						//System.err.println(count);
						counts.set(row, column, count);
					}
				}
			}
		}
		}
		close();
	}
	
	
	private void scoreDNADNA(MatrixWithHeaders counts, int binResolution, Kmer genes, boolean weight) {
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(genes, true)){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					String row=region1.toUCSC();
					String column=region2.toUCSC();
					if(!row.equals(column) && counts.containsColumn(column) && counts.containsRow(row)){
						double count=counts.get(row, column);
						double score=2.0/c.getAllDNAIntervals().size();
						//count++;
						if(weight){count+=score;}
						else{count++;}
						//System.err.println(count);
						counts.set(row, column, count);
					}
				}
			}
		}
		}
		close();
	}
	
	private void scoreDNADNA(MatrixWithHeaders counts, int binResolution, String gene, boolean weight) {
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(gene)){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					String row=region1.toUCSC();
					String column=region2.toUCSC();
					if(!row.equals(column) && counts.containsColumn(column) && counts.containsRow(row)){
						double count=counts.get(row, column);
						double score=2.0/c.getAllDNAIntervals().size();
						//count++;
						if(weight){count+=score;}
						else{count++;}
						//System.err.println(count);
						counts.set(row, column, count);
					}
				}
			}
		}
		}
		close();
	}
	
	private void scoreDNADNA(MatrixWithHeaders counts, int binResolution, boolean weight) {
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					String row=region1.toUCSC();
					String column=region2.toUCSC();
					if(counts.containsColumn(column) && counts.containsRow(row)){
						//System.err.println(row+" "+column);
						double count=counts.get(row, column);
						double score=2.0/c.getAllDNAIntervals().size();
						//count++;
						
						if(weight){
							count+=score;
						}
						else{
							count++;
						}
						//System.err.println(count);
						counts.set(row, column, count);
					}
				}
			}
			
		counter++;
		if(counter%100000==0) {System.err.println(counter);}
		}
		close();
		
		
		
	}
	


	/*private void scoreDNADNA(MatrixWithHeaders counts, int binResolution, String chr, String gene) {
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(gene)){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					String row=region1.toUCSC();
					String column=region2.toUCSC();
					if(!row.equals(column) && region1.getReferenceName().equals(chr) && region2.getReferenceName().equals(chr)){
						double count=counts.get(row, column);
						count++;
						//System.err.println(count);
						counts.set(row, column, count);
					}
				}
			}
		}
		}
		close();
	}*/



	public Pair<MatrixWithHeaders> getRNADNAContactMatrix(double minRNAFreq, int binResolution){
		MatrixWithHeaders scores= computeInput(binResolution, minRNAFreq);
		MatrixWithHeaders counts=new MatrixWithHeaders(scores.getRowNames(), scores.getColumnNames());
		
		Pair<MatrixWithHeaders> rtrn=score(scores, counts, binResolution);
		return rtrn;
	}
	
	
	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, String gene) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		//Collection<Cluster> clusters=this.getRNAClusters(gene);
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getRNANames().contains(gene)){
				Cluster binned=c.bin(binResolution);
				for(SingleInterval region: binned.getAllDNAIntervals()){
					double score=0;
					if(rtrn.containsKey(region)){
						score=rtrn.get(region);
					}
					score+=(2.0/binned.getAllDNAIntervals().size());
					rtrn.put(region, score);
				}
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		close();
		return rtrn;
	}
	
	public Map<SingleInterval, Double> getRNADNAClusterCounts(int binResolution, String gene) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		//Collection<Cluster> clusters=this.getRNAClusters(gene);
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getRNANames().contains(gene)){
				Cluster binned=c.bin(binResolution);
				for(SingleInterval region: binned.getAllDNAIntervals()){
					double score=0;
					if(rtrn.containsKey(region)){
						score=rtrn.get(region);
					}
					score++;
					rtrn.put(region, score);
				}
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		close();
		return rtrn;
	}
	
	
	public Map<String, Map<SingleInterval, Double>> getRNADNAContactsByList(int binResolution, Collection<String> genes, boolean weight) {
		Map<String, Map<SingleInterval, Double>> map=new TreeMap<String, Map<SingleInterval, Double>>();
		//Collection<Cluster> clusters=this.getRNAClusters(gene);
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			Collection<String> names=getRNAs(c, genes);
			if(names.size()>0){
				Cluster binned=c.bin(binResolution);
				for(String name: names) {
					if(!map.containsKey(name)) {map.put(name, new TreeMap<SingleInterval, Double>());}
					Map<SingleInterval, Double> rtrn=map.get(name);
					
					for(SingleInterval region: binned.getAllDNAIntervals()){
						double score=0;
						if(rtrn.containsKey(region)){
							score=rtrn.get(region);
						}
						double weightFactor=1;
						if(weight) {weightFactor=(2.0/binned.getAllDNAIntervals().size());}
						score+=weightFactor;
						rtrn.put(region, score);
					
					}
				}
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		close();
		return map;
	}
	

	
	
	private Collection<String> getRNAs(Cluster c, Collection<String> genes) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String g: genes) {
			if(c.getRNANames().contains(g)){rtrn.add(g);}
		}
		return rtrn;
	}



	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, String gene, String type, int min, int max) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()>=min && c.getClusterSize()<=max && c.containsRNA(gene, type)){
				Cluster binned=c.bin(binResolution);
				for(SingleInterval region: binned.getAllDNAIntervals()){
					double score=0;
					if(rtrn.containsKey(region)){
						score=rtrn.get(region);
					}
					//score+=(2.0/binned.getAllDNAIntervals().size());
					score++;
					rtrn.put(region, score);
				}
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		close();
		return rtrn;
	}
	
	
	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, String gene, int min, int max) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()>=min && c.getClusterSize()<=max && c.containsRNA(gene)){
				Cluster binned=c.bin(binResolution);
				for(SingleInterval region: binned.getAllDNAIntervals()){
					double score=0;
					if(rtrn.containsKey(region)){
						score=rtrn.get(region);
					}
					//score+=(2.0/binned.getAllDNAIntervals().size());
					score++;
					rtrn.put(region, score);
				}
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		close();
		return rtrn;
	}
	
	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, SingleInterval geneRegion, int min, int max) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()>=min && c.getClusterSize()<=max && c.containsRNA(geneRegion)){
				Cluster binned=c.bin(binResolution);
				for(SingleInterval region: binned.getAllDNAIntervals()){
					double score=0;
					if(rtrn.containsKey(region)){
						score=rtrn.get(region);
					}
					//score+=(2.0/binned.getAllDNAIntervals().size());
					score++;
					rtrn.put(region, score);
				}
				System.out.println(c.toString());
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		close();
		return rtrn;
	}
	
	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, Kmer kmer, boolean weight) {
		return getRNADNAContacts(binResolution, kmer, weight, true);
	}
	
	
	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, Collection<SingleInterval> rnaRegions, boolean weight) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		
		int counter=0;
		int specific=0;
		while(hasNext()){
			Cluster c=next();
			if(overlaps(c.getAllRNARegions(), rnaRegions)){
				Cluster binned=c.bin(binResolution);
				for(SingleInterval region: binned.getAllDNAIntervals()){
					double score=0;
					if(rtrn.containsKey(region)){
						score=rtrn.get(region);
					}
					if(weight){score+=(2.0/binned.getAllDNAIntervals().size());}
					else{score++;}
					rtrn.put(region, score);
				}
				specific++;
			}
			counter++;
			if(counter%10000==0){System.err.println(counter+" "+specific);}
		}
		close();
		
		
		
		
		return rtrn;
	}
	
	private boolean overlaps(Collection<RNAInterval> allRNARegions, Collection<SingleInterval> rnaRegions) {
		for(RNAInterval r1: allRNARegions){
			for(SingleInterval r2: rnaRegions){
				if(r1.overlaps(r2)){
					//System.err.println(r1.toUCSC()+" "+r1.getName()+" "+r2.toUCSC());
					return true;
				}
			}
		}
		return false;
	}


	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, Kmer kmer, boolean weight, boolean requireAll) {
		//Map<SingleInterval, Double> input=new TreeMap<SingleInterval, Double>();
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Collection<Cluster> clusters=this.getRNAClusters(kmer, requireAll);
		
		for(Cluster c: clusters){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				double score=0;
				if(rtrn.containsKey(region)){
					score=rtrn.get(region);
				}
				if(weight){score+=(2.0/binned.getAllDNAIntervals().size());}
				else{score++;}
				rtrn.put(region, score);
			}
		}
		return rtrn;
	}
	
	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, Kmer kmer, boolean weight, boolean requireAll, int min, int max) {
		//Map<SingleInterval, Double> input=new TreeMap<SingleInterval, Double>();
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Collection<Cluster> clusters=this.getRNAClusters(kmer, requireAll, min, max);
		
		System.err.println("Num clusters: "+clusters.size());
		
		for(Cluster c: clusters){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				double score=0;
				if(rtrn.containsKey(region)){
					score=rtrn.get(region);
				}
				if(weight){score+=(2.0/binned.getAllDNAIntervals().size());}
				else{score++;}
				rtrn.put(region, score);
			}
		}
		return rtrn;
	}
	
	
	
	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, SingleInterval geneCoordinates, boolean weight, int min, int max) {
		//Map<SingleInterval, Double> input=new TreeMap<SingleInterval, Double>();
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Collection<Cluster> clusters=this.getRNAClusters(geneCoordinates, min, max);
		
		System.err.println("Num clusters: "+clusters.size());
		
		for(Cluster c: clusters){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				double score=0;
				if(rtrn.containsKey(region)){
					score=rtrn.get(region);
				}
				if(weight){score+=(2.0/binned.getAllDNAIntervals().size());}
				else{score++;}
				rtrn.put(region, score);
			}
		}
		return rtrn;
	}

	private Pair<MatrixWithHeaders> score(MatrixWithHeaders scores, MatrixWithHeaders counts, int binResolution) {
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()<1000){
				Collection<String> rnas=c.getRNANames();
				score(scores, counts, rnas, c, binResolution);
				counter++;
			}
			if(counter%1000000 ==0){System.err.println(counter);}
			
		}
		close();
		Pair<MatrixWithHeaders> rtrn=new Pair<MatrixWithHeaders>(scores, counts);
		return rtrn;
	}
	
	private void score(MatrixWithHeaders scores, MatrixWithHeaders counts, Collection<String> rnas, Cluster c, int binResolution) {
		Cluster binned=c.bin(binResolution);
		double score=(2.0/(double)c.getAllDNAIntervals().size());
		
		for(SingleInterval region: binned.getAllDNAIntervals()){
			String column=region.toUCSC();
			for(String rna: rnas){
				if(scores.containsRow(rna)){
					double updatedScore=scores.get(rna, column)+score;
					scores.set(rna, region.toUCSC(), updatedScore);
					
					double count=counts.get(rna, column)+1;
					counts.set(rna, region.toUCSC(), count);
				}
			}
		}
		
	}
	
	
	
	
	
	

	private MatrixWithHeaders computeInput(int binResolution, double minRNAFreq) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Map<String, SingleInterval> rnaRegionMap=new TreeMap<String, SingleInterval>();
		Map<String, Double> geneExpression=new TreeMap<String, Double>();
		
		Collection<String> genes=new TreeSet<String>();
		Collection<String> genesToUse=new TreeSet<String>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			Collection<RNAInterval> rnaRegions=c.getAllRNARegions();
			updateCoordinates(rnaRegions, rnaRegionMap);
			genes.addAll(c.getRNANames());
			Collection<String> rnas=c.getRNANames();
			for(String rna: rnas){
				double countTotal=0;
				if(geneExpression.containsKey(rna)){countTotal=geneExpression.get(rna);}
				countTotal++;
				geneExpression.put(rna, countTotal);
				if(countTotal>=minRNAFreq){genesToUse.add(rna);}
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
		close();
		//double percentile=percentile(rtrn, 0.05);
		
		List<String> geneNames=getNames(genesToUse);
		List<String> genomePosition=getPositions(rtrn);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(geneNames, genomePosition);
		Map<String, String> rowDescriptions=convert(rnaRegionMap);
		mwh.setPIDToName(rowDescriptions);
		
		
		return mwh;
	}
	
	private Map<String, String> convert(Map<String, SingleInterval> rnaRegionMap) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String gene: rnaRegionMap.keySet()){
			SingleInterval pos=rnaRegionMap.get(gene);
			rtrn.put(gene, pos.toUCSC());
		}
		
		return rtrn;
	}

	private void updateCoordinates(Collection<RNAInterval> regions, Map<String, SingleInterval> rnaRegions) {
		for(SingleInterval rna: regions){
			String name=rna.getName();
			SingleInterval updated=rna;
			if(rnaRegions.containsKey(name)){
				updated=update(rnaRegions.get(name), rna);
			}
			rnaRegions.put(name, updated);
		}
	}
	
	private SingleInterval update(SingleInterval rna1, SingleInterval rna2) {
		String chr=rna1.getReferenceName();
		int start=Math.min(rna1.getReferenceStartPosition(), rna2.getReferenceStartPosition());
		int end=Math.max(rna1.getReferenceEndPosition(), rna2.getReferenceEndPosition());
		SingleInterval rtrn=new SingleInterval(chr, start, end);
		rtrn.setName(rna1.getName());
		return rtrn;
	}
	
	private List<String> getPositions(Map<SingleInterval, Double> inputCounts) {
		List<String> rtrn=new ArrayList<String>();
		for(SingleInterval region: inputCounts.keySet()){
			rtrn.add(region.toUCSC());
		}
		return rtrn;
	}
	
	private List<String> getNames(Collection<String> rnas) {
		List<String> rtrn=new ArrayList<String>();
		rtrn.addAll(rnas);
		return rtrn;
	}

	/*public Map<String, Integer> getRNACountsPerGene() {
		if(hasCountsPerGene){return this.countsPerGene;}
		
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
			
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			Collection<String> geneNames=c.getRNANames();
			for(String gene: geneNames){
				int count=0;
				if(rtrn.containsKey(gene)){count=rtrn.get(gene);}
				count++;
				rtrn.put(gene, count);
			}
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
			
		close();
		countsPerGene=rtrn;
		hasCountsPerGene=true;
		return rtrn;
	}*/
	
	public ArrayList<SingleInterval> getAllRNARegions() {
		ArrayList<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			Collection<RNAInterval> genes=c.getAllRNARegions();
			rtrn.addAll(genes);
			
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
			
		close();
		return rtrn;
		
	}

	public ArrayList<Integer> getRNAClusterSizes() {
		ArrayList<Integer> rtrn=new ArrayList<Integer>();
		
		while(hasNext()){
			Cluster c=next();
			rtrn.add(c.getAllRNARegions().size());
		}
		
		close();
		return rtrn;
	}

	/**
	 * Get clusters containing gene 1 and gene 2
	 * @param gene1
	 * @param gene2
	 * @return
	 */
	public Collection<Cluster> getRNAClusters(String gene1, String gene2) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(gene1) && c.containsRNA(gene2)){rtrn.add(c);}
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
		
		close();
		return rtrn;
	}

	public Collection<Cluster> getRNAClusters(String gene1) {
		return getRNAClusters(gene1, gene1);
	}
	
	
	public Collection<Cluster> getRNAClusters(Collection<Gene> genes) {
		Kmer kmer=new Kmer();
		for(Gene gene: genes){kmer.addRegion(gene.getName());}
		return getRNAClusters(kmer, false);
	}
	
	public Collection<Cluster> getRNAClusters(Kmer kmer) {
		return getRNAClusters(kmer, true);
	}
	
	
	public Collection<Cluster> getRNAClusters(Kmer kmer, boolean requireAll) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		int specific=0;
		while(hasNext()){
			Cluster c=next();
			boolean useCluster=useCluster(c, kmer, requireAll);
			if(useCluster){
				rtrn.add(c);
				specific++;
			}
			counter++;
			if(counter%1000000==0){System.err.println(counter+" "+specific);}
		}
		
		close();
		return rtrn;
	}
	
	
	public Collection<Cluster> getRNAClusters(SingleInterval geneCoordinates, int min, int max) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		int specific=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()>min && c.getClusterSize()<max){
				//System.err.println(c.getClusterSize());
				boolean useCluster=useCluster(c, geneCoordinates);
				if(useCluster){
					rtrn.add(c);
					specific++;
				}
			}
			counter++;
			if(counter%10000==0){System.err.println(counter+" "+specific);}
		}
		
		close();
		return rtrn;
	}
	
	public Collection<Cluster> getRNAClusters(Kmer kmer, boolean requireAll, int min, int max) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		int specific=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()>min && c.getClusterSize()<max){
				//System.err.println(c.getClusterSize());
				boolean useCluster=useCluster(c, kmer, requireAll);
				if(useCluster){
					rtrn.add(c);
					specific++;
				}
			}
			counter++;
			if(counter%10000==0){System.err.println(counter+" "+specific);}
		}
		
		close();
		return rtrn;
	}
	
	
	private boolean useCluster(Cluster c, SingleInterval geneCoordinates) {
		for(RNAInterval rna: c.getAllRNARegions()) {
			if(rna.overlaps(geneCoordinates)) {return true;}
		}
		return false;
	}
	
	private boolean useCluster(Cluster c, Kmer kmer, boolean requireAll) {
		if(kmer.isInput){return true;}
		
		//System.err.println(c.getClusterSize()+" "+c.getRNANames().size()+" "+kmer.regions.size());
		
		if(requireAll && c.getRNANames().containsAll(kmer.regions)){return true;}
		if(!requireAll){
		//if(kmer.getSize()>=c.getRNANames().size()){
				for(String rna: c.getRNANames()){
					//System.err.println(rna);
					if(kmer.getRegions().contains(rna)){return true;}
				}
			/*}
			else{
				for(String region: kmer.regions){
					if(c.containsRNA(region)){return true;}
				}
			}*/
		}
		
		
		return false;
		
	}

	

	public Collection<Cluster> getClustersOverlappingRegion(SingleInterval promoterRegion) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.containsOverlappingDNA(promoterRegion)){rtrn.add(c);}
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
		
		close();
		return rtrn;
	}
	
	public Collection<Cluster> getClustersOverlappingRegion(SingleInterval promoterRegion, int min, int max) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()>min && c.getClusterSize()<max){
				if(c.containsOverlappingDNA(promoterRegion)){rtrn.add(c);}
				counter++;
				//if(counter%100000==0){System.err.println(counter);}
			}
		}
		
		close();
		return rtrn;
	}
	
	public Collection<Cluster> getClustersOnChr(String chr) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.containsChromosome(chr)){rtrn.add(c);}
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
		
		close();
		return rtrn;
	}

	
	public Collection<Cluster> getClustersOnChr(String chr, int min, int max) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getClusterSize()>min && c.getClusterSize()<max){
				if(c.containsChromosome(chr)){rtrn.add(c);}
				counter++;
				if(counter%100000==0){System.err.println(counter);}
			}
		}
		
		close();
		return rtrn;
	}

	public Collection<String> getAllRNAs() {
		Collection<String> rtrn=new TreeSet<String>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			rtrn.addAll(c.getRNANames());
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		close();
		
		return rtrn;
	}
	
	public Map<String, Integer> getAllRNAList() {
		Map<String, Integer> map=new TreeMap<String, Integer>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			Collection<String> names=new TreeSet<String>();
			for(RNAInterval r: c.getAllRNARegions()) {
				if(r.hasExon()) {
					names.addAll(r.getExonNames());
				}
			}
			for(String name: names) {
				int count=0;
				if(map.containsKey(name)) {count=map.get(name);}
				count++;
				map.put(name, count);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		map.put("total", counter);
		
		close();
		return map;
	}

	public Collection<Cluster> getClustersOverlappingRNA(String rna) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.getRNANames().contains(rna)){rtrn.add(c);}
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
		
		close();
		return rtrn;
	}
	
	/*public Collection<Cluster> getClustersOverlappingRNA(String rna) throws IOException {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		List<long[]> list=indexKeys.get(rna);
		RandomAccessFile reader=new RandomAccessFile(this.barcodeFile, "r");
			
		for(long[] array: list){
			reader.seek(array[0]);
			byte[] b=new byte[(int) array[1]-1];
			reader.read(b);
			String s=new String(b);
			Cluster c=new Cluster(s);
			rtrn.add(c);
		}
			
		reader.close();
		return rtrn;
	}*/
		
	
	
	
	public BarcodingDataStreaming getSubsetOverlappingRNA(String rna) throws IOException {
		String save=this.barcodeFile.getAbsolutePath()+"."+rna+".clusters";
		
		if(new File(save).exists()){
			return new BarcodingDataStreaming(new File(save));
		}
		
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		int retained=0;
		while(hasNext()){
			Cluster c=next();
			//if(c.getClusterSize()>min && c.getClusterSize()<max){
				if(c.getRNANames().contains(rna)){
					writer.write(c.toString()+"\n");
					retained++;
				//}
			}
			counter++;
			if(counter%100000==0){System.err.println(counter+" "+retained);}
		}
		
		close();
		writer.close();
		
		System.err.println("Retained "+retained);
		
		return new BarcodingDataStreaming(new File(save));
	}

	public int numberOfClusters() {
		int counter=0;
		while(hasNext()){
			next();
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		close();
		return counter;
	}

	public Map<String, Integer> getCountsPerGene() {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		int total=0;
		while(hasNext()){
			Cluster c=next();
			for(String rna: c.getRNANames()){
				int count=0;
				if(rtrn.containsKey(rna)){count=rtrn.get(rna);}
				count++;
				rtrn.put(rna, count);
			}
			total++;
			if(total%100000 ==0){System.err.println(total);}
		}
		
		rtrn.put("Total", total);
		
		close();
		return rtrn;
	}

	public File getBarcodeFile() {
		return this.barcodeFile;
	}

	public MatrixWithHeaders getClusterRNAInteractionMatrix(Kmer kmer) {
		Collection<Cluster> clusters=getRNAClusters(kmer);
		MatrixWithHeaders mwh=makeClusterMatrix(clusters, 0.1);
		return mwh;
	}

	private MatrixWithHeaders makeClusterMatrix(Collection<Cluster> clusters, double freqCutoff) {
		List<String> rows=new ArrayList<String>();
		List<String> columns=new ArrayList<String>();
		
		Map<String, Integer> geneCount=new TreeMap<String, Integer>();
		
		for(Cluster c: clusters){
			rows.add(c.getBarcode());
			for(String gene:c.getRNANames()){
				int count=0;
				if(geneCount.containsKey(gene)){
					count=geneCount.get(gene);
				}
				count++;
				geneCount.put(gene, count);
			}
		}
		
		double cutoff=freqCutoff*clusters.size();
		for(String gene: geneCount.keySet()){
			int count=geneCount.get(gene);
			if(!gene.equals("Unassigned") && count>cutoff){
				//double ratio=(double)count/cutoff;
				//if(ratio>5){System.err.println(gene+" "+cutoff+" "+count+" "+ratio);}
				columns.add(gene);
			}
		}
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		for(Cluster c: clusters){
			for(String gene: columns){
				if(c.containsRNA(gene)){mwh.set(c.getBarcode(), gene, 1.0);}
			}
		}
		
		
		return mwh;
	}
	
	private MatrixWithHeaders makeClusterMatrix(Collection<Cluster> clusters, List<String> genesToUse, double freqCutoff) {
		List<String> rows=new ArrayList<String>();
		List<String> columns=genesToUse;
		
		for(Cluster c: clusters){
			rows.add(c.getBarcode());
		}
		
		
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		for(Cluster c: clusters){
			for(String gene: columns){
				if(c.containsRNA(gene)){mwh.set(c.getBarcode(), gene, 1.0);}
			}
		}
		
		Collection<String> list=new TreeSet<String>();
		for(String row: mwh.getRowNames()){
			double[] val=mwh.getRow(row);
			int sum=sum(val);
			if(sum>1){list.add(row);}
		}
		
		mwh=mwh.submatrixByRowNames(list);
		
		return mwh;
	}

	private int sum(double[] val) {
		int rtrn=0;
		for(int i=0; i<val.length; i++){rtrn+=val[i];}
		return rtrn;
	}

	public MatrixWithHeaders getClusterRNAInteractionMatrix(Collection<Cluster> clusters) {
		MatrixWithHeaders mwh=makeClusterMatrix(clusters, 0.01);
		return mwh;
	}
	
	
	public MatrixWithHeaders getClusterRNAInteractionMatrix(Collection<Cluster> clusters, double percentFilter) {
		MatrixWithHeaders mwh=makeClusterMatrix(clusters, percentFilter);
		return mwh;
	}
	
	public MatrixWithHeaders getClusterRNAInteractionMatrix(Collection<Cluster> clusters, List<String> genesToUse) {
		MatrixWithHeaders mwh=makeClusterMatrix(clusters, genesToUse, 0.01);
		return mwh;
	}

	public void addFilter(ClusterFilter filter){
		this.filters.add(filter);
	}

	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, Kmer kmer, double minPercentile, double maxPercentile) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Collection<Cluster> clusters=this.getRNAClusters(kmer);
		
		clusters=getClustersWithinPercentileRange(clusters, minPercentile, maxPercentile);
				
		for(Cluster c: clusters){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				double score=0;
				if(rtrn.containsKey(region)){
					score=rtrn.get(region);
				}
				score++;
				rtrn.put(region, score);
			}
		}
		return rtrn;
	}

	private Collection<Cluster> getClustersWithinPercentileRange(Collection<Cluster> clusters, double minPercentile, double maxPercentile) {
		double[] vals=new double[clusters.size()];
		int i=0;
		for(Cluster c: clusters){
			vals[i]=c.getClusterSize();
			i++;
		}
		double minClusterSize=Statistics.quantile(vals, minPercentile);
		double maxClusterSize=Statistics.quantile(vals, maxPercentile);
		
		System.err.println(minPercentile+" "+minClusterSize+" "+maxPercentile+" "+maxClusterSize);
		
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		for(Cluster c: clusters){
			if(c.getClusterSize()>minClusterSize && c.getClusterSize()<maxClusterSize){rtrn.add(c);}
		}
		return rtrn;
	}

	public Map<SingleInterval, Double> getRNADNAContacts(int binResolution, Kmer kmer, double minSize, double maxSize,	boolean weight) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
			
		int counter=0;
		int specific=0;
		while(hasNext()){
			Cluster c=next();
			if(kmer.isInput() || c.getRNANames().containsAll(kmer.regions)){
				if(c.getClusterSize()>minSize && c.getClusterSize()<maxSize){
					Cluster binned=c.bin(binResolution);
					for(SingleInterval region: binned.getAllDNAIntervals()){
						double score=0;
						if(rtrn.containsKey(region)){
							score=rtrn.get(region);
						}
						if(weight){score+=(2.0/binned.getAllDNAIntervals().size());}
						else{score++;}
						rtrn.put(region, score);
						}
					specific++;
					}
			}
			counter++;
			if(counter%10000==0){System.err.println(counter+" "+specific);}
		}
		
		close();
		
		
		return rtrn;
	}
	
	public interface ClusterFilter{
		public boolean evaluate(Cluster c);
	}

	public MatrixWithHeaders getRNADNAContacts(int binResolution, Collection<Gene> rnas) {
		Collection<Cluster> clusters=getRNAClusters(rnas);
		
		List<String> rows=new ArrayList<String>();
		for(Gene g: rnas){
			if(!rows.contains(g.getName())){rows.add(g.getName());}
		}
		
		List<String> columns=getGenomePositions(binResolution, clusters);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		int counter=0;
		for(Cluster c: clusters){
			//Cluster c=next();
			Cluster binned=c.bin(binResolution);
			for(String row: rtrn.getRowNames()){
				if(c.containsRNA(row)){
					for(SingleInterval region: binned.getAllDNAIntervals()){
						double val=rtrn.get(row, region.toUCSC());
						val++;
						rtrn.set(row,  region.toUCSC(), val);
					}
				}
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		close();
		
		return rtrn;
	}


	private List<String> getGenomePositions(int binResolution, Collection<Cluster> clusters) {
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		
		int counter=0;
		for(Cluster c: clusters){
			regions.addAll(c.bin(binResolution).getAllDNAIntervals());
			counter++;
			if(counter%1000==0){System.err.println(counter);}
		}
		
		
		List<String> rtrn=new ArrayList<String>();
		for(SingleInterval region: regions){
			rtrn.add(region.toUCSC());
		}
		
		return rtrn;
		
	}


	public Map<SingleInterval, Double> getRNACounts(String rna, boolean weighted) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		Map<String, SingleInterval> regions=new TreeMap<String, SingleInterval>();
		
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(rna.equalsIgnoreCase("Input")|| c.getRNANames().contains(rna)){
				for(RNAInterval region: c.getAllRNARegions()){
					String name=region.getName();
					double score=0;
					SingleInterval r=region;
					if(regions.containsKey(name)){
						r=merge(regions.get(name), region);
					}
					regions.put(name, r);
					if(rtrn.containsKey(name)){score=rtrn.get(name);}
					if(weighted){score+=(2.0/c.getClusterSize());}
					else{score++;}
					rtrn.put(name, score);
				}
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		close();
		
		Map<SingleInterval, Double> map=new TreeMap<SingleInterval, Double>();
		
		for(String name: rtrn.keySet()){
			double val=rtrn.get(name);
			SingleInterval region= regions.get(name);
			map.put(region,  val);
		}
		
		return map;
	}


	private SingleInterval merge(SingleInterval singleInterval, RNAInterval region) {
		String chr=singleInterval.getReferenceName();
		String name=region.getName();
		int start=Math.min(singleInterval.getReferenceStartPosition(), region.getReferenceStartPosition());
		int end=Math.max(singleInterval.getReferenceEndPosition(), region.getReferenceEndPosition());
		SingleInterval rtrn=new SingleInterval(chr, start, end);
		rtrn.setName(name);
		return rtrn;
	}


	public BarcodingDataStreaming shuffleRNA() throws IOException {
		File rtrnFile=new File(this.getBarcodeFile().getAbsolutePath()+".shuffle.clusters");
		return shuffleRNA(rtrnFile);
	}
	

	
	public IntervalTree<String> getProbabilityTree() throws IOException {
		
		//go through each cluster and get RNA reads
		Map<String, Integer> countsPerRNA=new TreeMap<String, Integer>();
		//List<RNAInterval> list=new ArrayList<RNAInterval>();
		
		int counter=0;
		while(hasNext()) {
			Cluster c=next();
			Collection<String> rnas=c.getRNANames();
			
			for(String rna: rnas) {
				increment(countsPerRNA, rna);
			}
			counter++;
			if(counter%1000000 ==0) {System.err.println(counter);}
		}
		
		close();
		
		IntervalTree<String> tree=makeProbabilityTree(countsPerRNA);
		
		return tree;
	}
	
	public IntervalTree<String> getProbabilityTree(Map<String, String> allRNA) throws IOException {
		
		//go through each cluster and get RNA reads
		Map<String, Integer> countsPerRNA=new TreeMap<String, Integer>();
		
		int counter=0;
		while(hasNext()) {
			Cluster c=next();
		
			Cluster r=c.renameRNA(allRNA);
			Collection<String> rnas=r.getRNANames();
			
			for(String rna: rnas) {
				increment(countsPerRNA, rna);
			}
			counter++;
			if(counter%1000000 ==0) {System.err.println(counter);}
		}
		
		close();
		
		IntervalTree<String> tree=makeProbabilityTree(countsPerRNA);
		
		return tree;
	}
	
	public BarcodingDataStreaming shuffleRNA(File rtrnFile) throws IOException {
		
		FileWriter writer=new FileWriter(rtrnFile);
		
		
		//go through each cluster and get RNA reads
		Map<String, Integer> countsPerRNA=new TreeMap<String, Integer>();
		//List<RNAInterval> list=new ArrayList<RNAInterval>();
		
		int totalSize=0;
		int counter=0;
		while(hasNext()) {
			Cluster c=next();
			Collection<String> rnas=c.getRNANames();
			
			for(String rna: rnas) {
				increment(countsPerRNA, rna);
				totalSize++;
			}
			counter++;
			if(counter%1000000 ==0) {System.err.println(counter);}
		}
		
		close();
		
		IntervalTree<String> tree=makeProbabilityTree(countsPerRNA);
		
		
		//rebuild clusters by randomly sampling from universe of reads
		// go cluster by cluster and repopulate
		
		counter=0;
		while(hasNext()) {
			Cluster c=next();
			int size=c.getAllRNARegions().size();
			Cluster newCluster=new Cluster(c.getBarcode());
			//newCluster.addDNAReads(c.getAllDNAIntervals());
			Collection<RNAInterval> rnas=getRandom(tree, size, totalSize);
			newCluster.addRNAReads(rnas);
			writer.write(newCluster.toString()+"\n");
			counter++;
			if(counter%1000000 ==0) {System.err.println(counter); writer.flush();}
		}
	
		writer.close();
		
		return new BarcodingDataStreaming(rtrnFile);
	}
	
	private IntervalTree<String> makeProbabilityTree(Map<String, Integer> countsPerRNA) {
		IntervalTree<String> rtrn=new IntervalTree<String>();
		int currentPos=0;
		
		for(String RNA: countsPerRNA.keySet()) {
			int length=countsPerRNA.get(RNA);
			int endPos=currentPos+length;
			rtrn.put(currentPos, endPos, RNA);
			//System.err.println(currentPos+" "+endPos+" "+RNA);
			currentPos=endPos;
		}
		return rtrn;
	}



	private void increment(Map<String, Integer> countsPerRNA, String rna) {
		int count=0;
		if(countsPerRNA.containsKey(rna)) {count=countsPerRNA.get(rna);}
		count++;
		countsPerRNA.put(rna, count);
	}



	/*public BarcodingDataStreaming shuffleRNA(File rtrnFile) throws IOException {
		
		FileWriter writer=new FileWriter(rtrnFile);
		
		
		//go through each cluster and get RNA reads
		List<RNAInterval> list=new ArrayList<RNAInterval>();
		
		int counter=0;
		while(hasNext()) {
			Cluster c=next();
			Collection<RNAInterval> intervals=c.getAllRNARegions();
			list.addAll(intervals);
			counter++;
			if(counter%1000000 ==0) {System.err.println(counter);}
		}
		
		close();
		
		//rebuild clusters by randomly sampling from universe of reads
		// go cluster by cluster and repopulate
		
		counter=0;
		while(hasNext()) {
			Cluster c=next();
			int size=c.getAllRNARegions().size();
			Cluster newCluster=new Cluster(c.getBarcode());
			//newCluster.addDNAReads(c.getAllDNAIntervals());
			Collection<RNAInterval> rnas=getRandom(list, size);
			newCluster.addRNAReads(rnas);
			writer.write(newCluster.toString()+"\n");
			counter++;
			if(counter%1000000 ==0) {System.err.println(counter); writer.flush();}
		}
	
		writer.close();
		
		return new BarcodingDataStreaming(rtrnFile);
	}*/


	private Collection<RNAInterval> getRandom(IntervalTree<String> tree, int size, int totalSize) {
		Collection<RNAInterval> rtrn=new TreeSet<RNAInterval>();
		for(int i=0; i<size; i++) {
			double random=Math.random();
			int pos=(int)(totalSize*random);
			String rna=tree.overlappingValueIterator(pos, pos+1).next();
			RNAInterval temp=new RNAInterval(new SingleInterval(rna, 0, 1));
			temp.setName(rna);
			rtrn.add(temp);
		}
		return rtrn;
	}



	public void writeClustersOverlappingRegion(SingleInterval locus, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		int counter=0;
		int retained=0;
		while(hasNext()){
			Cluster c=next();
			if(c.containsOverlappingDNA(locus)){
				retained++;
				writer.write(c.toString()+"\n");
			}
			counter++;
			if(counter%100000==0){System.err.println(counter+" "+retained);}
		}
		
		reset();
		writer.close();
	}



	public void getClusters(Collection<String> geneNames, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		int counter=0;
		int retained=0;
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(geneNames)){
				writer.write(c.toString()+"\n");
				retained++;
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter+" "+retained);}
		}
		close();
		writer.close();
		
	}



	public Collection<SingleInterval> getRegions(String gene) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		int counter=0;
		while(hasNext()){
			Cluster c=next();
			if(c.containsRNA(gene)){
				rtrn.addAll(c.getAllDNAIntervals());
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter+" "+rtrn.size());}
		}
		close();
		return rtrn;
	}
	

	/*private Collection<RNAInterval> getRandom(List<RNAInterval> list, int size) {
		Collection<RNAInterval> rtrn=new TreeSet<RNAInterval>();
		for(int i=0; i<size; i++) {
			double random=Math.random();
			int pos=(int)(list.size()*random);
			rtrn.add(list.get(pos));
		}
		return rtrn;
	}*/
	
}

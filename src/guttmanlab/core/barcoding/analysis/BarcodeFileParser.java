package guttmanlab.core.barcoding.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import net.sf.samtools.util.CloseableIterator;

/**
 * 
 * @author mitch
 * Static methods to parse barcoding files
 */
public class BarcodeFileParser {
	
	
	public static Map<String, IntervalTree<String>> parseBarcodeFile(File barcodeFile) throws NumberFormatException, IOException{
		Collection<File> barcodeFiles=new ArrayList<File>();
		barcodeFiles.add(barcodeFile);
		return parseBarcodeFile(barcodeFiles.iterator(), null);
	}
	
	public static Map<String, IntervalTree<String>> parseBarcodeFile(Iterator<File> barcodeFiles) throws NumberFormatException, IOException{
		return parseBarcodeFile(barcodeFiles, null);
	}
	
	public static Map<String, IntervalTree<String>> parseBarcodeFile(Iterator<File> barcodeFile, Annotation region) throws NumberFormatException, IOException {
		return parseBarcodeFile(barcodeFile, region, -99, -99);
	}

	public static Map<String, IntervalTree<String>> parseBarcodeFile(Iterator<File> barcodeFile, Annotation region, int minClusterSize, int maxClusterSize) throws NumberFormatException, IOException{
		Collection<Annotation> regions=new ArrayList<Annotation>();
		if(region!=null){regions.add(region);}
		return(parseBarcodeFile(barcodeFile, regions, minClusterSize, maxClusterSize));
	}
		
	
	/**
	 * Retain clusters that have ALL regions
	 * @param barcodeFile
	 * @param regions
	 * @param minClusterSize
	 * @param maxClusterSize
	 * @return
	 * @throws NumberFormatException
	 * @throws IOException
	 */
	public static Map<String, IntervalTree<String>> parseBarcodeFile(Iterator<File> barcodeFiles, Collection<Annotation> regions, int minClusterSize, int maxClusterSize) throws NumberFormatException, IOException{
		Map<String, IntervalTree<String>> rtrn=new TreeMap<String, IntervalTree<String>>();
		
		while(barcodeFiles.hasNext()){
			File barcodeFile=barcodeFiles.next();
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
			String nextLine;
			while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
				String[] tokens=nextLine.split("\t");
				//System.err.println(tokens.length+" "+nextLine);
				
				String barcode=tokens[0];
				boolean containsAll=false;
				if(regions.isEmpty()){containsAll=true;}
				
				
				boolean overlaps=false;
				boolean passClusterSize=false;
				Collection<SingleInterval> temp=new TreeSet<SingleInterval>();
				
				if(minClusterSize== -99 || maxClusterSize == -99){passClusterSize=true;}
				else if(tokens.length-1>=minClusterSize && tokens.length-1<=maxClusterSize){passClusterSize=true;}
				
				for(int i=1; i<tokens.length; i++){
					
					String chr=tokens[i].split(":")[0];
					int start=new Integer(tokens[i].split(":")[1]);
					int end=start+1;
					SingleInterval interval=new SingleInterval(chr, start, end); //This is the region in the current cluster
					temp.add(interval);
				}
				
				boolean hasAllRegions=true;
				Map<String, IntervalTree<Annotation>> localTree=makeTree(temp);
				//iterate through regions
				for(Annotation region: regions){
					//check if they overlap tree
					if (localTree.containsKey(region.getReferenceName())){
						IntervalTree<Annotation> tree=localTree.get(region.getReferenceName());
						Iterator<Node<Annotation>> overlappers=tree.overlappers(region.getReferenceStartPosition(), region.getReferenceEndPosition());
						if(!overlappers.hasNext()){hasAllRegions=false;}
					}
				}
				
				
				if(hasAllRegions && passClusterSize){
					for(SingleInterval interval: temp){
						IntervalTree<String> tree=new IntervalTree<String>();
						if(rtrn.containsKey(interval.getReferenceName())){
							tree=rtrn.get(interval.getReferenceName());
						}
						tree.put(interval.getReferenceStartPosition(), interval.getReferenceEndPosition(), barcode);
						rtrn.put(interval.getReferenceName(), tree);
					}
				}
			}
			
			reader.close();
		}
		return rtrn;
	}
	

	private static Map<String, IntervalTree<Annotation>> makeTree(Collection<SingleInterval> temp) {
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		for(Annotation interval: temp){
			String chr=interval.getReferenceName();
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			tree.put(interval.getReferenceStartPosition(), interval.getReferenceEndPosition(), interval);
			rtrn.put(chr, tree);
		}
		return rtrn;
	}

	public static List<SingleInterval> getAllGenomicPositions(File barcodeFile) throws IOException {
		List<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		String nextLine;
		while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
			String[] tokens=nextLine.split(" ");
			for(int i=1; i<tokens.length; i++){
				String chr=tokens[i].split(":")[0];
				int start=new Integer(tokens[i].split(":")[1]);
				int end=start+1;
				SingleInterval interval=new SingleInterval(chr, start, end);
				rtrn.add(interval);
			}
		}
		
		reader.close();
		return rtrn;
	}

	public static Map<String, Collection<SingleInterval>> getBarcodeToPositions(File barcodeFile) throws IOException {
		return getBarcodeToPositions(barcodeFile, "all");
	}
	
	public static Map<String, Collection<SingleInterval>> getBarcodeToPositions(Iterator<File> barcodeFiles) throws IOException {
		Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
		while(barcodeFiles.hasNext()){
			File barcodeFile=barcodeFiles.next();
			Map<String, Collection<SingleInterval>> map=getBarcodeToPositions(barcodeFile, "all");
			rtrn.putAll(map);
		}
		
		return rtrn;
	}


	public static Map<String, Collection<SingleInterval>> getBarcodeToPositions(File barcodeFile, String chrToUse) throws IOException {
		Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		String nextLine;
		int counter=0;
		while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
			//System.err.println(nextLine);
			String[] tokens=nextLine.split("\t");
			boolean hasChr=false;
			
			String barcode=tokens[0];
			List<SingleInterval> list=new ArrayList<SingleInterval>();
			
			for(int i=1; i<tokens.length; i++){
				if(tokens[i].contains(":")){
					String chr=tokens[i].split(":")[0];
					hasChr=hasChr || chr.equalsIgnoreCase(chrToUse);
					String position=tokens[i].split(":")[1];
					int start;
					int end;
					if(position.contains("-")){
						start=new Integer(position.split("-")[0]);
						end=new Integer(position.split("-")[1]);
					}
					else{
						start=new Integer(position);
						end=start+1;
					}
					
					SingleInterval interval=new SingleInterval(chr, start, end);
					list.add(interval);
				}
			}
			
			if(hasChr || chrToUse.equalsIgnoreCase("all")){
				rtrn.put(barcode, list);
				//if(list.size()>1 && list.size()<10){System.out.println(nextLine);}
			}
			
			counter++;
			if(counter %10000000 ==0 ){System.err.println(counter);}
		}
		//System.err.println("Read in "+counter +" clusters"+" "+rtrn.keySet().size());
		
		reader.close();
		return rtrn;
	}
	
	public static Map<String, Collection<SingleInterval>> getBarcodeToPositions(File barcodeFile, SingleInterval region) throws IOException {
		Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		String nextLine;
		int counter=0;
		while((nextLine=reader.readLine())!=null && (nextLine.trim().length()>0)){
			String[] tokens=nextLine.split("\t");
			boolean overlapsBin=false;
			
			String barcode=tokens[0];
			List<SingleInterval> list=new ArrayList<SingleInterval>();
			
			for(int i=1; i<tokens.length; i++){
				String chr=tokens[i].split(":")[0];
				int start=new Integer(tokens[i].split(":")[1]);
				int end=start+1;
				SingleInterval interval=new SingleInterval(chr, start, end);
				overlapsBin=overlapsBin || interval.overlaps(region);
				list.add(interval);
			}
			
			if(overlapsBin){
				rtrn.put(barcode, list);
				//if(list.size()>1 && list.size()<10){System.out.println(nextLine);}
			}
			counter++;
			if(counter %1000000 ==0 ){System.err.println(counter);}
		}
		System.err.println("Read in "+counter +" clusters"+" "+rtrn.keySet().size());
		
		reader.close();
		return rtrn;
	}

	public static Collection<String> getBarcodesOverlappingRegion(File barcodeFile, Annotation gene, int maxClusterSize) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		Map<String, Collection<SingleInterval>> barcodes=getBarcodeToPositions(barcodeFile);
		for(String barcode: barcodes.keySet()){
			Collection<SingleInterval> intervals=barcodes.get(barcode);
			if(intervals.size()<maxClusterSize){
				boolean overlaps=false;
				for(SingleInterval interval: intervals){
					if(interval.overlaps(gene)){overlaps=true;}
				}
				if(overlaps){rtrn.add(barcode);}
			}
		}
		return rtrn;
	}
	
	public static void preprocessFile(File barcodeFile, String saveDir, Map<String, Integer> chrSizes, int minSize, int maxSize) throws IOException{
		new File(saveDir).mkdir();
		/*for(String chr: chrSizes.keySet()){
			BarcodingData data=new BarcodingData(barcodeFile, chr);
			System.err.println("loaded data for "+chr);
			Collection<SingleInterval> bins=getBins(chr, chrSizes.get(chr), 1000000);
			for(SingleInterval region: bins){
				String save=saveDir+"/"+region.getFileName();
				//get all clusters overlapping this region
				Map<String, Cluster> clusters=data.getClustersOverlappingRegionByBarcodeName(region);
				System.err.println("Got clusters "+clusters.keySet().size());
				write(save, clusters);
				System.err.println("Wrote "+region.toUCSC()+" to "+save);
			}
		}*/
		
		BarcodingDataStreaming data=new BarcodingDataStreaming(barcodeFile);
		
		//TODO Keep all file handles open and then close at once
		Map<String, FileWriter> writers=new TreeMap<String, FileWriter>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()>=minSize && c.getClusterSize()<maxSize){
				Collection<String> files=new TreeSet<String>();
				for(SingleInterval region: c.getAllIntervals()){
					files.add(makeFile(region, 1000000, saveDir));
				}
				writers=write(files, writers, c);
			}
			//else{System.err.println("Skipped "+c.getBarcode()+" "+c.getClusterSize());}
			counter++;
			if(counter%100000 ==0){System.err.println(counter+" "+c.getBarcode()+" "+c.getClusterSize());}
		}
		
		data.close();
		close(writers);
	}

	private static void close(Map<String, FileWriter> writers) throws IOException {
		for(String file: writers.keySet()){
			writers.get(file).close();
		}
	}

	private static void write(String save, Map<String, Cluster> clusters) throws IOException {
		FileWriter writer= new FileWriter(save);
		
		for(String barcode: clusters.keySet()){
			Cluster c=clusters.get(barcode);
			if(c.getAllIntervals().size()<1000){
				writer.write(c+"\n");
			}
			//else{System.err.println("Skipped "+c.getBarcode()+" "+c.getAllIntervals().size());}
		}
		
		writer.close();
	}

	private static Collection<SingleInterval> getBins(String currentChr, Integer size, int resolution){
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		for(int i=0; i<size; i+=resolution){
			SingleInterval interval=new SingleInterval(currentChr, i, i+resolution);
			rtrn.add(interval);
		}
		return rtrn;
	}

	private static String makeFile(SingleInterval region, int resolution, String saveDir) throws IOException {
		int startIndex=region.getReferenceStartPosition()/resolution;
		int newStart=startIndex*resolution;
		int newEnd=(startIndex+1)*resolution;
		SingleInterval binned=new SingleInterval(region.getReferenceName(), newStart, newEnd);
		return (saveDir+"/"+binned.getFileName());
	}

	private static Map<String, FileWriter> write(Collection<String> files, Map<String, FileWriter> writers, Cluster c) throws IOException {
		for(String file: files){
			FileWriter writer;
			if(writers.containsKey(file)){
				writer=writers.get(file);
			}
			else{
				writer=new FileWriter(file);
				writers.put(file, writer);
			}
			//writer.write(c.toSPRITEFormat()+"\n");
			writer.write(c.getReadString()+"\n");
			//System.out.println(c.getReadString());
		}
		return writers;
	}
	
	public static Map<String, IntervalTree<String>> parseBarcodeFile(File data, Annotation region, int minClusterSize, int maxClusterSize) throws NumberFormatException, IOException {
		
		Collection<File> barcodingFiles=new ArrayList<File>();
		barcodingFiles.add(data);
		
		return parseBarcodeFile(barcodingFiles.iterator(), region, minClusterSize, maxClusterSize);
	}

	public static Map<String, IntervalTree<String>> parseBarcodeFile(File data, AnnotationCollection<Gene> regions, int minClusterSize, int maxClusterSize) throws NumberFormatException, IOException {
		Collection<Annotation> annotations=new TreeSet<Annotation>();
		CloseableIterator<Gene> iter=regions.sortedIterator();
		while(iter.hasNext()){
			Gene gene=iter.next();
			annotations.add(gene);
		}
		Collection<File> barcodingFiles=new ArrayList<File>();
		barcodingFiles.add(data);
		
		return parseBarcodeFile(barcodingFiles.iterator(), annotations, minClusterSize, maxClusterSize);
	}

	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File barcodeData=new File(args[0]);
			String saveDir=args[1];
			File sizesFile=new File(args[2]);
			int minSize=2;
			int maxSize=1000;
			if(args.length>3){
				minSize=new Integer(args[3]);
				maxSize=new Integer(args[4]);
			}
			
			
			Map<String, Integer> sizes=CoordinateSpace.getRefSeqLengthsFromTable(sizesFile.getAbsolutePath());
			
			BarcodeFileParser.preprocessFile(barcodeData, saveDir, sizes, minSize, maxSize);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Barcoding File \n args[1]=save directory \n args[2]=sizes file \n args[3]=min cluster size (optional, default=2) \n args[4]=max cluster size (optional, default=1000)";
}

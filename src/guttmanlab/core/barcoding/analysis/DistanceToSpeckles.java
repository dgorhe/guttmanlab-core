package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;
import net.sf.samtools.util.CloseableIterator;

public class DistanceToSpeckles {

	private int minClusterSize=1;
	private int maxClusterSize=1000;
	private int numPerm=10;
	
	public DistanceToSpeckles(Map<String, IntervalTree<File>> barcodeFileTree, Collection<SingleInterval> activeHubRegions, AnnotationCollection<Gene> genes) throws IOException{
		//BarcodingData data=getActiveHubRegions(barcodeFileTree, activeHubRegions);
		
		System.err.println("loaded data");
		
		Iterator<Gene> geneIter=genes.sortedIterator();
		while(geneIter.hasNext()){
			Gene gene=geneIter.next();
			Annotation geneSI=new SingleInterval(gene.getReferenceName(), gene.getReferenceStartPosition(), gene.getReferenceEndPosition());
			Iterator<File> iter=barcodeFileTree.get(gene.getReferenceName()).overlappingValueIterator(gene.getReferenceStartPosition(), gene.getReferenceEndPosition());
			BarcodingData data=new BarcodingData(iter);
			System.err.println(gene.getName());
			//For each gene calculate distance to each AH region
			//double[] vals=new double[activeHubRegions.size()];
			//int i=0;
			
			CloseableIterator<DerivedAnnotation<? extends Annotation>> geneWindows=geneSI.getWindows(10000, 10000).sortedIterator();
			
			while(geneWindows.hasNext()){
				DerivedAnnotation geneWindow=geneWindows.next();
			
			
				Collection<Double> vals=new ArrayList<Double>();
				
				for(Annotation ah: activeHubRegions){
					CloseableIterator<DerivedAnnotation<? extends Annotation>> windows=ah.getWindows(1000000, 1000000).sortedIterator();
					while(windows.hasNext()){
						DerivedAnnotation window=windows.next();
						double count=data.quantify(window,  geneWindow);
						vals.add(count);
					}
				}
				double average=Statistics.mean(vals);
				System.out.println(geneWindow.getReferenceName()+"\t"+geneWindow.getReferenceStartPosition()+"\t"+geneWindow.getReferenceEndPosition()+"\t"+average);
			}
		}
		
		
	}
	
	
	public DistanceToSpeckles(Map<String, IntervalTree<File>> barcodeFileTree, Collection<SingleInterval> activeHubRegions, int resolution, CoordinateSpace space, String chr) throws IOException{
		//BarcodingData data=getActiveHubRegions(barcodeFileTree, activeHubRegions);
		
		System.err.println("loaded data V2");
		
		SingleInterval chrRegion=new SingleInterval(chr, 0, space.getRefSizes().get(chr));
			
		CloseableIterator<DerivedAnnotation<? extends Annotation>> geneWindows=chrRegion.getWindows(resolution, resolution).sortedIterator();
		
		int total=chrRegion.size()/resolution;
		int counter=0;
		
		while(geneWindows.hasNext()){
			DerivedAnnotation window=geneWindows.next();
			Iterator<File> iter=barcodeFileTree.get(window.getReferenceName()).overlappingValueIterator(window.getReferenceStartPosition(), window.getReferenceEndPosition());
			BarcodingData data=new BarcodingData(iter);
		
			Collection<Double> vals=new ArrayList<Double>();
			for(Annotation ah: activeHubRegions){
				double count=data.quantify(window,  ah);
				vals.add(count);
			}
			double average=Statistics.mean(vals);
			System.out.println(window.getReferenceName()+"\t"+window.getReferenceStartPosition()+"\t"+window.getReferenceEndPosition()+"\t"+average);
			System.err.println(counter+" "+total+" "+average);
			counter++;
		}
		
		
	}
	
	public DistanceToSpeckles(BarcodingData data, Collection<SingleInterval> activeHubRegions, Collection<SingleInterval> regions, int resolution) throws IOException{
		//BarcodingData data=getActiveHubRegions(barcodeFileTree, activeHubRegions);
		
		System.err.println("loaded data");
		
		for(SingleInterval geneSI: regions){
			int total=geneSI.getLength()/resolution;
			
			CloseableIterator<DerivedAnnotation<? extends Annotation>> geneWindows=geneSI.getWindows(resolution, resolution).sortedIterator();
			
			int counter=0;
			while(geneWindows.hasNext()){
				DerivedAnnotation geneWindow=geneWindows.next();
				Collection<Double> vals=new ArrayList<Double>();
				Collection<Cluster> clusters=data.getClustersOverlappingRegion(geneWindow);
				//System.err.println(geneWindow.toUCSC()+" "+clusters.size());
				
				for(Annotation ah: activeHubRegions){
					double count=data.quantify(clusters, ah);
					//System.err.println(geneWindow.toUCSC()+" "+count);
					vals.add(count);
				}
				double average=Statistics.mean(vals);
				/*System.out.print(geneWindow.getReferenceName()+"\t"+geneWindow.getReferenceStartPosition()+"\t"+geneWindow.getReferenceEndPosition());
				for(Double val: vals){
					System.out.print("\t"+val);
				}
				System.out.print("\n");
				System.out.flush();*/
				System.out.println(geneWindow.getReferenceName()+"\t"+geneWindow.getReferenceStartPosition()+"\t"+geneWindow.getReferenceEndPosition()+"\t"+average);
				
				System.err.println(counter+" "+total+" "+average);
				counter++;
			}
		}
		
		
	}
	
	
	public DistanceToSpeckles(BarcodingDataStreaming data, Cluster hub, int resolution, CoordinateSpace space, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		
		
		for(String chr: space.getRefSizes().keySet()){
			runByChr(writer, chr, space, resolution, data, hub);
		}
		writer.close();
	}
	
	public DistanceToSpeckles(BarcodingDataStreaming data, Cluster hub, int resolution, CoordinateSpace space, String save, String chr) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		//TODO Pull clusters overlapping speckle and nucleolus directly
		runByChr(writer, chr, space, resolution, data, hub);
		
		writer.close();
	}
	

	private void runByChr(FileWriter writer, String chr, CoordinateSpace space, int resolution, BarcodingDataStreaming data, Cluster hub) throws IOException {
		SingleInterval chrRegion=new SingleInterval(chr, 0, space.getRefSizes().get(chr));
		System.err.println(chr);
		Collection<Cluster> randomHubs=getPermutations(hub, numPerm, space);
		
		Collection<Cluster> chrClusters=data.getClustersOverlappingRegion(new SingleInterval(chr, 0, space.getRefSizes().get(chr)));
		System.err.println("chr "+chrClusters.size());
		
		Collection<Annotation> allRegions=new TreeSet<Annotation>();
		allRegions.addAll(hub.getAllIntervals());
		for(Cluster randomHub: randomHubs){
			allRegions.addAll(randomHub.getAllIntervals());
		}
		
		chrClusters=BarcodingData.getClustersOverlappingRegion(chrClusters, allRegions);
		System.err.println("bodies "+chrClusters.size());
		
		
		int total=chrRegion.getLength()/resolution;
			
		CloseableIterator<DerivedAnnotation<? extends Annotation>> geneWindows=chrRegion.getWindows(resolution, resolution).sortedIterator();
			
		int counter=0;
		while(geneWindows.hasNext()){
			DerivedAnnotation geneWindow=geneWindows.next();
			writer.write(geneWindow.getReferenceName()+"\t"+geneWindow.getReferenceStartPosition());
			
			
			//Collection<Cluster> clusters=data.getClustersOverlappingRegion(geneWindow, minClusterSize, maxClusterSize);
			
			//Collection<Cluster> sc=data.getClustersOverlappingRegion(geneWindow);
			Collection<Cluster> sc=BarcodingData.getClustersOverlappingRegion(chrClusters, geneWindow);
			
			System.err.println("window "+sc.size());
			
			
			double hubDistance=hub.computeInterchromosomalDistance(geneWindow, sc);
			
			System.err.println(chr+" "+counter+" "+total+" "+hubDistance);
			writer.write("\t"+hubDistance);
			
			
			for(Cluster randomHub: randomHubs){
				double randomScore=randomHub.computeInterchromosomalDistance(geneWindow, sc);
				System.err.println(randomScore);
				writer.write("\t"+randomScore);
			}
			writer.write("\n");
			
			counter++;
		}
		
	}


	private Collection<Cluster> getPermutations(Cluster hub, int numPerm, CoordinateSpace space) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		for(int i=0;i<numPerm; i++){
			rtrn.add(hub.getPermutedCluster(space.getRefSizes()));
		}
		return rtrn;
	}


	private BarcodingData getActiveHubRegions(Map<String, IntervalTree<File>> barcodeFileTree, Collection<SingleInterval> activeHubRegions) throws IOException {
		Collection<File> files=new TreeSet<File>();
		
		for(Annotation region: activeHubRegions){
			Iterator<File> iter=barcodeFileTree.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			while(iter.hasNext()){
				File file=iter.next();
				files.add(file);
			}
		}
		return new BarcodingData(files.iterator());
	}
	
	
	public DistanceToSpeckles(BarcodingDataStreaming data, Collection<SingleInterval> activeHubRegions, Collection<SingleInterval> nucleolarRegions){
		Cluster speckleHub=new Cluster("Active", activeHubRegions);
		Cluster nucleolarHub=new Cluster("Inactive", nucleolarRegions);
		
		Map<SingleInterval, Double> speckleDistance=data.distance(speckleHub);
		Map<SingleInterval, Double> nucleolarDistance=data.distance(nucleolarHub);
		
		System.err.println(speckleDistance+" "+nucleolarDistance);
	}
	
	
	private static Map<String, IntervalTree<File>> makeFileTree(File[] files) {
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
	
	private static SingleInterval parse(String name) {
		String chr=name.split("_")[0];
		int start=new Integer(name.split("_")[1]);
		int end=new Integer(name.split("_")[2]);
		return new SingleInterval(chr, start, end);
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			
			Collection<SingleInterval> activeHubRegions=BEDFileIO.loadSingleIntervalFromFile(args[1]);
			
			Cluster activeHub=new Cluster("Active", activeHubRegions);
			Collection<SingleInterval> nucleolarRegions=BEDFileIO.loadSingleIntervalFromFile(args[2]);
			String save=args[3];
			int resolution=new Integer(args[4]);
			
			String genome=args[5];
			
			String chr=args[6];
			
			CoordinateSpace coordinateSpace=null;
			if(genome.equalsIgnoreCase("hg19")){coordinateSpace=CoordinateSpace.HG19;}
			if(genome.equalsIgnoreCase("mm10")){coordinateSpace=CoordinateSpace.MM10;}
			if(genome.equalsIgnoreCase("mm9")){coordinateSpace=CoordinateSpace.MM9;}
			
			if(coordinateSpace==null){throw new IllegalArgumentException("No genome");}
			
			
			new DistanceToSpeckles(data, activeHub, resolution, coordinateSpace, save, chr);
		
			
		}
		else{
			System.err.println(usage);
		}
	}

	static String usage=" args[0]=barcoding data (Files) \n args[1]=speckle regions \n args[2]=nucleolar regions \n args[3]=save \n args[4]=resolution \n args[5]=genome (mm9, mm10, hg19)";
	
}

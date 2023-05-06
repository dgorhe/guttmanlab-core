package guttmanlab.core.barcoding.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.Statistics;
import htsjdk.samtools.util.CloseableIterator;

public class CommandLineTools {
	
	private static void computeMaxKmer(String[] args) throws IOException{
		String usage=" args[0]=regions (BED) \n args[1]=cluster files \n args[2]=chromosoem sizes \n args[3]=save \n args[4]=max cluster size \n args[5]=num perms \n args[6]=kmer size \n args[7]=resolution";
		
		if(args.length>6){
			Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[0]);
			File[] files=new File(args[1]).listFiles();
			Map<String, IntervalTree<File>> barcodeFileTree=PermuteKmers.makeFileTree(files);
			String sizes=args[2];
			String save=args[3];
			int maxClusterSize=new Integer(args[4]);
			Map<String, Integer> chrSizes=CoordinateSpace.getRefSeqLengthsFromTable(sizes);
			
			int numPerm=100;
			if(args.length>5){
				numPerm=new Integer(args[5]);
			}
			
			int kmer=new Integer(args[6]);
			int resolution=new Integer(args[7]);
			
			PermuteKmers permute=new PermuteKmers(barcodeFileTree, chrSizes, maxClusterSize);			
			
			permute.getRandomizedCountsForMaxKmer(regions, kmer, numPerm, save, true);
		}
		else{System.err.println(usage);}
	}
	
	private static void quantifyKmer(String[] args) throws IOException{
		String usage=" args[0]=regions (BED) \n args[1]=cluster files \n args[2]=chromosoem sizes \n args[3]=max cluster size \n args[4]=kmer size \n args[5]=save \n args[6]=only consider interchromosomal (true, false) \n args[7]=num permutations";
		
		if(args.length>5){
			Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[0]);
			File[] files=new File(args[1]).listFiles();
			Map<String, IntervalTree<File>> barcodeFileTree=PermuteKmers.makeFileTree(files);
			String sizes=args[2];
			int maxClusterSize=new Integer(args[3]);
			Map<String, Integer> chrSizes=CoordinateSpace.getRefSeqLengthsFromTable(sizes);
			
			int kmerSize=new Integer(args[4]);
			
			FileWriter writer=new FileWriter(args[5]);
			
			boolean onlyInter=new Boolean(args[6]);
			
			int numPerm=new Integer(args[7]);
	
			PermuteKmers permute=new PermuteKmers(barcodeFileTree, chrSizes, maxClusterSize);
			
			Collection<Cluster> clusters=permute.getUniqueSPRITEClusters(regions, kmerSize, onlyInter);

			for(Cluster c: clusters){
				double observed=permute.weightedQuantify(c);
				double rawCount=permute.quantify(c).size();
				double[] random=permute.randomWeighted(c, numPerm, chrSizes, barcodeFileTree);
				double expected=Statistics.mean(random);
				double percentile=Statistics.percentLessThan(observed, random);
				writer.write(c.toStringNoName()+"\t"+rawCount+"\t"+observed+"\t"+(observed/expected)+"\t"+percentile+"\n");
			}
			
			writer.close();
			System.err.println("done");
			
			/*Cluster c=new Cluster("Observed", regions);
			
			
			double observedScore=permute.weightedQuantify(c);
			System.out.println(c.toStringNoName()+" "+permute.quantify(c).size()+" "+observedScore);

			double[] randomScores=permute.randomWeighted(c, 100, chrSizes, barcodeFileTree);
			
			double expected=Statistics.mean(randomScores);
			double max=Statistics.max(randomScores);
			System.err.println(observedScore+" "+expected+" "+max+" "+(observedScore/expected));
			*/
			
		}
		else{System.err.println(usage);}
	}
	
	
	private static void coverageInteractingWithRegion(String[] args) throws IOException{
		String usage=" args[0]=region \n args[1]=cluster files \n args[2]=chromosoem sizes \n args[3]=save \n args[4]=max cluster size \n args[5]=resolution";

		
		//TODO This should NOT need a resolution or kmer
		if(args.length>5){
			SingleInterval region=new SingleInterval(args[0]);
			File[] files=new File(args[1]).listFiles();
			Map<String, IntervalTree<File>> barcodeFileTree=PermuteKmers.makeFileTree(files);
			String sizes=args[2];
			String save=args[3];
			int maxClusterSize=new Integer(args[4]);
			Map<String, Integer> chrSizes=CoordinateSpace.getRefSeqLengthsFromTable(sizes);
			int resolution=new Integer(args[5]);
			
			PermuteKmers permute=new PermuteKmers(barcodeFileTree, chrSizes, maxClusterSize, 1);
			
			//Collection<Cluster> rawSpriteClusters=permute.getAllUniqueSPRITEClusters(regions, kmer, false);
			
			Collection<Cluster> rawSpriteClusters=permute.getAllUniqueSPRITEClusters(region);
			System.err.println("Number of unique SPRITE clusters: "+rawSpriteClusters.size());
			
			writeBEDGraph(save, rawSpriteClusters, resolution);
			
			//permute.writeBEDStack(save, rawSpriteClusters, resolution);
		}
		else{System.err.println(usage);}
	}
	

	private static void writeBEDGraph(String save, Collection<Cluster> rawSpriteClusters, int resolution) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(Cluster c: rawSpriteClusters){
			Cluster binned=c.bin(resolution);
			for(SingleInterval interval: binned.getAllIntervals()){
				int count=0;
				if(rtrn.containsKey(interval)){
					count=rtrn.get(interval);
				}
				count++;
				rtrn.put(interval, count);
			}
		}
		
		for(SingleInterval region: rtrn.keySet()){
			int count=rtrn.get(region);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+count+"\n");
		}
		
		writer.close();
	}

	private static void writeIGV(String[] args) throws IOException{
		String usage=" args[0]=kmer string \n args[1]=cluster files \n args[2]=chromosoem sizes \n args[3]=save \n args[4]=max cluster size \n args[5]=resolution";

		
		//TODO This should NOT need a resolution or kmer
		if(args.length>5){
			Collection<SingleInterval> regions=parseKmer(args[0]);
			File[] files=new File(args[1]).listFiles();
			Map<String, IntervalTree<File>> barcodeFileTree=PermuteKmers.makeFileTree(files);
			String sizes=args[2];
			String save=args[3];
			int maxClusterSize=new Integer(args[4]);
			Map<String, Integer> chrSizes=CoordinateSpace.getRefSeqLengthsFromTable(sizes);
			int resolution=new Integer(args[5]);
			
			PermuteKmers permute=new PermuteKmers(barcodeFileTree, chrSizes, maxClusterSize, 1);
			
			//Collection<Cluster> rawSpriteClusters=permute.getAllUniqueSPRITEClusters(regions, kmer, false);
			
			Collection<Cluster> rawSpriteClusters=permute.getAllUniqueSPRITEClusters(regions);
			System.err.println("Number of unique SPRITE clusters: "+rawSpriteClusters.size());
			permute.writeIGV(save+".igv", rawSpriteClusters, resolution);
			permute.writeClusters(rawSpriteClusters, save+".clusters");
			
			//permute.writeBEDStack(save, rawSpriteClusters, resolution);
		}
		else{System.err.println(usage);}
	}
	
	private static Collection<SingleInterval> parseKmer(String kmerString) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		String[] tokens=kmerString.split(";");
		for(int i=0; i<tokens.length; i++){
			if(tokens[i].startsWith("chr")){
				SingleInterval interval=new SingleInterval(tokens[i]);
				rtrn.add(interval);
			}
		}
		return rtrn;
	}

	private static void randomizeOnSpecificClusters(String[] args) throws IOException{
		String usage=" args[0]=clusters (cluster file) \n args[1]=barcoding files \n args[2]=sizes \n args[3]=save \n arg[4]=max cluster size \n args[5]=num perm (default=100) \n args[6]=write .igv files (default=false) \n args[7]=resolution";

		if(args.length>4){
			
			FileWriter writer=new FileWriter(args[0]+".temp");
			Collection<Cluster> temp=parse(new File(args[0]), 1000000);
			for(Cluster c: temp){
				writer.write(c.toKmerString()+"\n");
			}
			writer.close();
			
			BarcodingDataStreaming clusters=new BarcodingDataStreaming(new File(args[0]+".temp"));
			
			File[] files=new File(args[1]).listFiles();
			Map<String, IntervalTree<File>> barcodeFileTree=PermuteKmers.makeFileTree(files);
			String sizes=args[2];
			String save=args[3];
			int maxClusterSize=new Integer(args[4]);
			Map<String, Integer> chrSizes=CoordinateSpace.getRefSeqLengthsFromTable(sizes);
			
			int numPerm=100;
			if(args.length>5){
				numPerm=new Integer(args[5]);
			}
			
			boolean writeIGV=false;
			if(args.length>6){
				writeIGV=new Boolean(args[6]);
			}
			
			int resolution=new Integer(args[7]);
			
			PermuteKmers permute=new PermuteKmers(barcodeFileTree, chrSizes, maxClusterSize);
			
			/*for(Cluster c: temp){
				int count=permute.quantify(c).size();
				System.err.println(c.toKmerString()+" "+count);
			}*/
		
			permute.computePermutationsForSpecificRegion(clusters, numPerm, save, resolution);
		
			if(writeIGV){
				permute.writeClusters(clusters, save);
			}
		
			System.out.println("Done");
		}
		
		else{System.err.println(usage);}
	}


	private static void permutationsByDistance(String[] args) throws IOException{
		String usage=" args[0]=clusters (cluster file) \n args[1]=barcoding files \n args[2]=sizes \n args[3]=save \n arg[4]=max cluster size \n args[5]=num perm \n args[6]=min count to consider";

		if(args.length>6){
			
			BarcodingDataStreaming clusters=new BarcodingDataStreaming(new File(args[0]));
			File[] files=new File(args[1]).listFiles();
			Map<String, IntervalTree<File>> barcodeFileTree=PermuteKmers.makeFileTree(files);
			String sizes=args[2];
			String save=args[3];
			int maxClusterSize=new Integer(args[4]);
			int numPerm=new Integer(args[5]);
			int minCount=new Integer(args[6]);
			
			Map<String, Integer> chrSizes=CoordinateSpace.getRefSeqLengthsFromTable(sizes);
			
			PermuteKmers permute=new PermuteKmers(barcodeFileTree, chrSizes, maxClusterSize);
			permute.computePermutationsByDistance(clusters, numPerm, save, minCount);
			System.err.println("Done");
		}
		else{System.err.println(usage);}
	}


	
	private static Collection<Cluster> parse(File file, int resolution) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine=reader.readLine();
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split(",");
			String[] positions=tokens[0].split(" ");
			Cluster cluster=new Cluster(nextLine);
			for(int i=0; i<positions.length; i++){
				int start=new Integer(positions[i].split("-")[1])*resolution;
				int end=start+resolution;
				SingleInterval region=new SingleInterval(positions[i].split("-")[0], start, end);
				cluster.addRead(region);
			}
			rtrn.add(cluster);
			counter++;
		}
		reader.close();
		return rtrn;
		
		
		
		
	}

	public static void main(String[] args) throws IOException{
		//permutationsByDistance(args);
		quantifyKmer(args);
		//computeMaxKmer(args);
		//writeIGV(args);
		
		//coverageInteractingWithRegion(args);
		//randomizeOnSpecificClusters(args);
	}
	
}

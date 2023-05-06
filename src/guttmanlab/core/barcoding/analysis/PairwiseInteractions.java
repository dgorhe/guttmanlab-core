package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SingleInterval;

public class PairwiseInteractions {
	
	int resolution=100000;

	public PairwiseInteractions(BarcodingDataFileTree data, Collection<SingleInterval> kmer, String save) throws IOException{
		Collection<SingleInterval> excluded=extend(kmer, 0*resolution);
		Collection<Cluster> clusters=data.getClustersOverlappingMultipleRegions(kmer);
		BarcodingData subset=new BarcodingData(clusters);
		subset=subset.bin(resolution);
		subset.computePairwiseInteractions(save, kmer.iterator().next().getReferenceName());
	}
	
	/*public PairwiseInteractions(BarcodingDataFileTree data, Collection<SingleInterval> kmer, String save) throws IOException{
		
		for(SingleInterval region: kmer){
			System.err.println(region.toUCSC());
			Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
			regions.add(region);
			Collection<SingleInterval> excluded=extend(regions, 0*resolution);
			Collection<Cluster> clusters=data.getClustersOverlappingMultipleRegions(regions);
			BarcodingData subset=new BarcodingData(clusters);
			subset=subset.bin(resolution);
			subset.computePairwiseInteractions(save+"_"+region.getFileName()+".txt", region.getReferenceName());
		}
		
		writeKmer(save+".gmt", kmer);
		
	}*/
	
	private void writeKmer(String save, Collection<SingleInterval> kmer) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("Kmer\tKmer");
		for(SingleInterval region: kmer){
			Iterator<? extends Annotation> windows=region.getWindows(resolution, resolution).sortedIterator();
			while(windows.hasNext()){
				Annotation window=windows.next();
				writer.write("\t"+window.toUCSC());
			}	
			
		}
		writer.write("\n");
		
		writer.close();
	}

	public PairwiseInteractions(BarcodingDataStreaming data, String chr, String save) throws IOException{
		Collection<Cluster> clusters=new TreeSet<Cluster>();
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.containsChr(chr)){clusters.add(c);}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		BarcodingData subset=new BarcodingData(clusters);
		subset=subset.bin(resolution);
		subset.computePairwiseInteractions(save, chr);
		
	}
	
	private Collection<SingleInterval> extend(Collection<SingleInterval> regions, int i) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval region: regions){
			String chr=region.getReferenceName();
			int newStart=region.getReferenceStartPosition()-i;
			int newEnd=region.getReferenceEndPosition()+i;
			rtrn.add(new SingleInterval(chr, newStart, newEnd));
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException{
		BarcodingDataFileTree data=new BarcodingDataFileTree(new File(args[0]).listFiles());
		String save=args[1];
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		for(int i=2; i<args.length; i++){
			regions.add(new SingleInterval(args[i]));
		}
		new PairwiseInteractions(data, regions, save);
		
		/*BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String chr=args[1];
		String save=args[2];
		new PairwiseInteractions(data, chr, save);*/
		
		
	}
	
}

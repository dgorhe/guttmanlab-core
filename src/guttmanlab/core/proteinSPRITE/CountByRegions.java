package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.simulation.CoordinateSpace;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class CountByRegions {

	public CountByRegions(BarcodingDataStreaming data, Map<String, IntervalTree<SingleInterval>> regions, String save) throws IOException{
		Map<SingleInterval, Integer> map=new TreeMap<SingleInterval, Integer>();
		
		int count=0;
		while(data.hasNext()){
			Cluster c=data.next();
			count+=c.getAllDNAIntervals().size();
			Collection<SingleInterval> overlappingRegions= getOverlaps(c, regions);
			add(overlappingRegions, map);
		}
		
		double observed=sum(map);
		double expected=expected(regions, count);
		System.err.println(count+" "+observed+" "+expected);
		data.close();
		write(save, map);
	}
	
	
	public CountByRegions(File bam, Map<String, IntervalTree<SingleInterval>> regions, String save) throws IOException{
		Map<SingleInterval, Integer> map=new TreeMap<SingleInterval, Integer>();
		
		int count=0;
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				count++;
				Collection<SingleInterval> overlappingRegions= getOverlaps(read, regions);
				add(overlappingRegions, map);
			}
		}
		
		double observed=sum(map);
		double expected=expected(regions, count);
		System.err.println(count+" "+observed+" "+expected);
		reads.close();
		inputReader.close();
		write(save, map);
	}

	private double sum(Map<SingleInterval, Integer> map) {
		double rtrn=0;
		for(SingleInterval region: map.keySet()) {rtrn+=map.get(region);}
		return rtrn;
	}


	private double expected(Map<String, IntervalTree<SingleInterval>> regions, int count) {
		double peaks=0;
		double genome=0;
		for(String chr: regions.keySet()) {
			peaks+=size(regions.get(chr));
			genome+=CoordinateSpace.HG19.getRefSizes().get(chr);
		}
		
		double rtrn=((double)count/genome)*peaks;
		return rtrn;
	}

	private double size(IntervalTree<SingleInterval> intervalTree) {
		double sum=0;
		Iterator<SingleInterval> iter=intervalTree.valueIterator();
		while(iter.hasNext()) {
			sum+=iter.next().size();
		}
		return sum;
	}
	
	private Collection<SingleInterval> getOverlaps(SAMRecord c, Map<String, IntervalTree<SingleInterval>> trees) {
		return overlappers(toSI(c), trees);
	}

	private SingleInterval toSI(SAMRecord c) {
		SingleInterval rtrn=new SingleInterval(c.getReferenceName(), c.getAlignmentStart(), c.getAlignmentEnd());
		return rtrn;
	}


	private Collection<SingleInterval> getOverlaps(Cluster c, Map<String, IntervalTree<SingleInterval>> trees) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval region: c.getAllDNAIntervals()) {
			rtrn.addAll(overlappers(region, trees));
		}
		
		return rtrn;
	}

	private Collection<SingleInterval> overlappers(SingleInterval region, Map<String, IntervalTree<SingleInterval>> trees) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		if(trees.containsKey(region.getReferenceName())) {
			Iterator<SingleInterval> iter=trees.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			while(iter.hasNext()) {
				rtrn.add(iter.next());
			}
		}
		return rtrn;
	
	}

	private void add(Collection<SingleInterval> overlappingRegions, Map<SingleInterval, Integer> map) {
		for(SingleInterval region: overlappingRegions) {
			int count=0;
			if(map.containsKey(region)) {count=map.get(region);}
			count++;
			map.put(region, count);
		}
	}

	private void write(String save, Map<SingleInterval, Integer> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			writer.write(region.toBedgraph(map.get(region))+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		
		//File bam=new File(args[0]);
		Map<String, IntervalTree<SingleInterval>> trees=BEDFileIO.loadSingleIntervalTree(args[1]);
		
		String save=args[2];
		
		new CountByRegions(data, trees, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=cluster file \n args[1]=bed regions \n args[2]=save";
}

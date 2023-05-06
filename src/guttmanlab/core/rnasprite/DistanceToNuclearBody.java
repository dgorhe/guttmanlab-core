package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class DistanceToNuclearBody {

	
	/*public DistanceToNuclearBody(BarcodingDataStreaming data, Kmer nuclearBody, String save, int binResolution) throws IOException {
		this.binResolution=binResolution;
		nuclearBody=bin(nuclearBody, binResolution);
		
		Collection<Cluster> clusters=data.getDNAClusters(nuclearBody, binResolution);
		
		Map<SingleInterval, Double> counts=getCounts(nuclearBody, clusters);
		
		Map<SingleInterval, Double> allCounts=coverage(data);
		
		Map<SingleInterval, Double> norm=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval r: counts.keySet()) {
			double val=counts.get(r)/allCounts.get(r);
			norm.put(r, val);
		}
		
		BEDFileIO.writeBEDGraph(allCounts, save+".norm.bedgraph");
		BEDFileIO.writeBEDGraph(norm, save+".ratio.bedgraph");
		BEDFileIO.writeBEDGraph(counts, save+".bedgraph");
	}*/
	
	
	public DistanceToNuclearBody(BarcodingDataStreaming data, Kmer nuclearBody, String save, int binResolution) throws IOException {
		Map<SingleInterval, Double> counts=getCounts(data, nuclearBody, binResolution);
		
		Map<SingleInterval, Double> norm=normalizeToHub(counts, nuclearBody, binResolution);
		
		BEDFileIO.writeBEDGraph(counts, save+".bedgraph");
		BEDFileIO.writeBEDGraph(norm, save+".norm.bedgraph");
	}
	
	
	private Map<SingleInterval, Double> normalizeToHub(Map<SingleInterval, Double> counts, Kmer nuclearBody, int binSize) {
		nuclearBody=bin(nuclearBody, binSize);
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();

		List<Double> values=new ArrayList<Double>();
		for(SingleInterval region: nuclearBody.getIntervals()) {
			if(counts.containsKey(region)) {
				double val=counts.get(region);
				values.add(val);
			}
		}
		
		double denom=Statistics.quantile(values, 0.5);
		
		for(SingleInterval r: counts.keySet()) {
			double val=counts.get(r);
			double norm=val/denom;
			rtrn.put(r, norm);
		}
				
				
		return rtrn;
	}


	public DistanceToNuclearBody(BarcodingDataStreaming data, Kmer nuclearBody, String save, int binResolution, String chr, String gene) throws IOException {
		nuclearBody=bin(nuclearBody, binResolution);
		
		Collection<Cluster> clusters=data.getDNAClusters(nuclearBody, binResolution);
		clusters=filterClusters(clusters, gene);
		
		Map<SingleInterval, Double> counts=getCounts(nuclearBody, clusters);
		
		BEDFileIO.writeBEDGraph(counts, save+".bedgraph", chr);
	}
	
	public DistanceToNuclearBody(BarcodingDataStreaming data, Kmer nuclearBody, String save, int binResolution, String chr) throws IOException {
		nuclearBody=bin(nuclearBody, binResolution);
		
		Collection<Cluster> clusters=data.getDNAClusters(nuclearBody, binResolution);
		
		Map<SingleInterval, Double> counts=getCounts(nuclearBody, clusters);
		
		BEDFileIO.writeBEDGraph(counts, save+".bedgraph", chr);
	}
	
	
	private Collection<Cluster> filterClusters(Collection<Cluster> clusters, String gene) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Cluster c: clusters) {
			if(c.containsRNA(gene)) {rtrn.add(c);}
		}
		
		return rtrn;
	}


	

	private Kmer bin(Kmer nuclearBody, int binResolution2) {
		Kmer rtrn=new Kmer();
		
		for(SingleInterval region: nuclearBody.getIntervals()) {
			Collection<SingleInterval> set=getRegions(region, binResolution2);
			rtrn.addIntervals(set);
		}
		
		return rtrn;
	}
	
	
	
	
	private Collection<SingleInterval> getRegions(SingleInterval region, int binResolution) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(int i=region.getReferenceStartPosition(); i<region.getReferenceEndPosition(); i+=binResolution){
			SingleInterval temp=new SingleInterval(region.getReferenceName(), i, i+binResolution);
			rtrn.add(temp);
			//System.err.println(region1.toUCSC()+" "+region.toUCSC());
		}
		return rtrn;
	}
	
	
	private Map<SingleInterval, Double> getCounts(BarcodingDataStreaming data, Kmer nuclearBody, int binSize) {
		nuclearBody=bin(nuclearBody, binSize);
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		
		while(data.hasNext()) {
			Cluster c=data.next();
			c=c.bin(binSize);
			Cluster body=hitsInBody(c, nuclearBody);
			for(SingleInterval region: c.getAllDNAIntervals()) {
				//TODO ensure that this cluster contains nuclear body region besides ones on same chromosome
				boolean hasNonChr=hasNonChr(body, region.getReferenceName());
				if(hasNonChr) {
					double count=0;
					if(rtrn.containsKey(region)) {count=rtrn.get(region);}
					//count++;
					count+=(2.0/c.getClusterSize());
					rtrn.put(region, count);
				}
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		data.close();
		return rtrn;
	}
	
	private Map<SingleInterval, Double> getCounts(Kmer nuclearBody, Collection<Cluster> clusters) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(Cluster c: clusters) {
			Cluster body=hitsInBody(c, nuclearBody);
			for(SingleInterval region: c.getAllDNAIntervals()) {
				//TODO ensure that this cluster contains nuclear body region besides ones on same chromosome
				boolean hasNonChr=hasNonChr(body, region.getReferenceName());
				if(hasNonChr) {
					double count=0;
					if(rtrn.containsKey(region)) {count=rtrn.get(region);}
					//count++;
					count+=(2.0/c.getClusterSize());
					rtrn.put(region, count);
				}
			}
		}
		return rtrn;
	}
	
	private Cluster hitsInBody(Cluster c, Kmer nuclearBody) {
		Cluster rtrn=new Cluster("body");
		
		for(SingleInterval region: c.getAllDNAIntervals()) {
			if(nuclearBody.getIntervals().contains(region)) {rtrn.addDNARead(region);}
		}
		
		return rtrn;
	}
	
	private boolean hasNonChr(Cluster body, String chr) {
		for(SingleInterval region: body.getAllDNAIntervals()) {
			if(!region.getReferenceName().equals(chr)) {return true;}
		}
		return false;
	}
	
	
	public static void main(String[] args) throws IOException {
		if(args.length>3) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Kmer kmer=new Kmer();
			kmer.addIntervals(BEDFileIO.loadSingleIntervalFromFile(args[1]));
			String save=args[2];
			int binResolution=Integer.parseInt(args[3]);
			
			if(args.length>5) {
				String chr=args[4];
				String gene=args[5];
				new DistanceToNuclearBody(data, kmer, save, binResolution, chr, gene);
			}
			else if(args.length>4) {
				String chr=args[4];
				new DistanceToNuclearBody(data, kmer, save, binResolution, chr);
			}
			else {
				new DistanceToNuclearBody(data, kmer, save, binResolution);
			}
		}
		else {
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=data \n args[1]=kmer \n args[2]=save \n args[3]=bin resolution \n args[4]=chr (optional) \n args[5]=RNA (optional)";
}

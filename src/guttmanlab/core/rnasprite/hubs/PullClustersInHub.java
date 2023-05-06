package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.rnasprite.Kmer;

public class PullClustersInHub {
	int binResolution=1000000;
	
	public PullClustersInHub(BarcodingDataStreaming data, Kmer nuclearBody, String save) throws IOException {
		nuclearBody=bin(nuclearBody, binResolution);
		
		data.getOriginalDNAClusters(nuclearBody, binResolution, save);
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
	
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Kmer kmer=new Kmer();
			kmer.addIntervals(BEDFileIO.loadSingleIntervalFromFile(args[1]));
			String save=args[2];
			new PullClustersInHub(data, kmer, save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=kmer \n args[2]=save";
}

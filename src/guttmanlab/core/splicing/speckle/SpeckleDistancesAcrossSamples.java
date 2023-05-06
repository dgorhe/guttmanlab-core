package guttmanlab.core.splicing.speckle;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Kmer;

public class SpeckleDistancesAcrossSamples {
	
	private static MatrixWithHeaders makeMatrix(File[] files, Map<SingleInterval, Double>[] scores) {
		Collection<SingleInterval> allRegions=new TreeSet<SingleInterval>();
		List<String> columns=new ArrayList<String>();
		for(int i=0; i<scores.length; i++) {
			allRegions.addAll(scores[i].keySet());
			columns.add(files[i].getName());
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(getRow(allRegions), columns);
		
		for(int i=0; i<files.length; i++) {
			String col=files[i].getName();
			for(SingleInterval region: allRegions) {
				double score=get(scores[i],region);
				rtrn.set(region.toUCSC(), col, score);
			}
		}
		return rtrn;
	}

	
	private static List<String> getRow(Collection<SingleInterval> allRegions) {
		List<String> rtrn=new ArrayList<String>();
		
		for(SingleInterval r: allRegions) {rtrn.add(r.toUCSC());}
		
		return rtrn;
	}


	private static double get(Map<SingleInterval, Double> map, SingleInterval region) {
		if(map.containsKey(region)) {return map.get(region);}
		return 0;
	}


	public static void main(String[] args) throws IOException {
		if(args.length>3) {
		File[] files=new File(args[0]).listFiles();
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[1]);
		String save=args[2];
		
		Kmer kmer=new Kmer();
		kmer.addIntervals(regions);
		

		int binResolution= Integer.parseInt(args[3]);
		
		Map<SingleInterval, Double>[] scores=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i]);
			BarcodingDataStreaming data=new BarcodingDataStreaming(files[i]);
			scores[i]=DistanceToNuclearBody.scoreBins(data, kmer, binResolution)[0];
		}
		
		MatrixWithHeaders mwh=makeMatrix(files, scores);
		//MatrixWithHeaders normalized=mwh.rankNormalize();
		
		mwh.write(save);
		
		System.err.println("done");
		}
		else {System.err.println(usage);}
		
	}
	
	

	static String usage=" args[0]=data files \n args[1]=kmer \n args[2]=save \n args[3]=bin resolution";
	
}

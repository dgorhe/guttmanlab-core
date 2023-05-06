package guttmanlab.core.splicing.speckle;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class CountWindows {

	private static Map<SingleInterval, Integer> getCounts(Collection<SingleInterval> regions, int binResolution) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(SingleInterval r: regions) {
			Collection<SingleInterval> bins=r.allBins(binResolution);
			for(SingleInterval bin: bins) {
				int count=0;
				if(rtrn.containsKey(bin)) {count=rtrn.get(bin);}
				count++;
				rtrn.put(bin, count);
			}
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[0]);
		int binResolution=Integer.parseInt(args[1]);
		int threshold=Integer.parseInt(args[2]);
		Map<SingleInterval, Integer> counts=getCounts(regions, binResolution);
		
		for(SingleInterval r: counts.keySet()) {
			int count=counts.get(r);
			if(count>threshold) {System.out.println(r.toShortBED());}
		}
		
	}

	
	
}

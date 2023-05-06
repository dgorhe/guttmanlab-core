package guttmanlab.core.rnasprite.hubs;

import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class MakeBins {
	
	private static Collection<SingleInterval> getBins(SingleInterval region, int binResolution) {
		Collection<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		
		rtrn.addAll(getRegions(region, binResolution));
		
		
		return rtrn;
	}
	
	private static Collection<SingleInterval> getRegions(SingleInterval region, int binResolution) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(int i=region.getReferenceStartPosition(); i<region.getReferenceEndPosition(); i+=binResolution){
			SingleInterval temp=new SingleInterval(region.getReferenceName(), i, i+binResolution);
			rtrn.add(temp);
			//System.err.println(region1.toUCSC()+" "+region.toUCSC());
		}
		return rtrn;
	}

	public static void main(String[] args) {
		SingleInterval region=new SingleInterval(args[0]);
		int binResolution=Integer.parseInt(args[1]);
		Collection<SingleInterval> bins=getBins(region, binResolution);
		
		System.out.println("Name\tUse");
		for(SingleInterval bin: bins) {
			System.out.println(bin.toUCSC()+"\tuse");
		}
		
	}
	
}

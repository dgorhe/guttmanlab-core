package guttmanlab.core.rnasprite;

import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class ContactsAtDefinedRegions {
	
	int binResolution=1000000;

	public ContactsAtDefinedRegions(BarcodingDataStreaming data, Collection<SingleInterval> regions) {
		Collection<SingleInterval> bins=getBins(regions);
		
		
	}
	
	private Collection<SingleInterval> getBins(Collection<SingleInterval> regionsCombined) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval region: regionsCombined){
			rtrn.addAll(get1MbRegions(region));
		}
		
		return rtrn;
	}
	
	private Collection<SingleInterval> get1MbRegions(SingleInterval region) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(int i=region.getReferenceStartPosition(); i<region.getReferenceEndPosition(); i+=binResolution){
			SingleInterval temp=new SingleInterval(region.getReferenceName(), i, i+binResolution);
			rtrn.add(temp);
			//System.err.println(region1.toUCSC()+" "+region.toUCSC());
		}
		return rtrn;
	}
	
}

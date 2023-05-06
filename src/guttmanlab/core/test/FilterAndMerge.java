package guttmanlab.core.test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.barcoding.analysis.BarcodingDataStreaming;
import guttmanlab.core.barcoding.analysis.Cluster;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class FilterAndMerge {
	
	int resolution=60;

	static String usage=" args[0]=cluster file \n args[1]=resolution (60, default) \n args[2]=save";
	
	public static void main(String[] args) throws IOException{
		if(args.length<2){System.err.println(usage);}
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		int resolution=new Integer(args[1]);
		String save=args[2];
		//FileWriter writer=new FileWriter(save);
		
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()>1){
				Cluster binned=c.bin(resolution);
				if(binned.getClusterSize()>1){
					regions.addAll(binned.getAllIntervals());
					//writer.write(binned.toString()+"\n");
				}
			}
		}
		data.close();
		//writer.close();
		
		
		List<String> rows=getNames(regions);
		List<String> columns=getNames(regions);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()>1){
				Cluster binned=c.bin(resolution);
				if(binned.getClusterSize()>1){
					Collection<SingleInterval> intervals=binned.getAllIntervals();
					for(SingleInterval region1: intervals){
						for(SingleInterval region2: intervals){
							if(region1!=region2){
								double count=mwh.get(region1.toUCSC(), region2.toUCSC());
								count++;
								mwh.set(region1.toUCSC(), region2.toUCSC(), count);
							}
						}
					}
				}
			}
		}
		data.close();
		mwh.write(save);
	}

	private static List<String> getNames(Collection<SingleInterval> regions) {
		List<String> rtrn=new ArrayList<String>();
		
		for(SingleInterval region: regions){
			rtrn.add(region.toUCSC());
		}
		
		return rtrn;
	}
	
}

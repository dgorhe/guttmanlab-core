package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.simulation.CoordinateSpace;

public class HUHContactMatrix {

	private static Collection<Cluster> parse(String string) throws IOException {
		Collection<String> lines=BEDFileIO.loadLines(string);
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		for(String line: lines){
			Collection<String> huhTags=parseHUH(line);
			Cluster c=new Cluster(line.split("\t")[0]);
			Collection<SingleInterval> regions=parseCluster(line);
			c.addDNAReads(regions);
			rtrn.add(c);	
		}
		
		return rtrn;
	}
	
	
	private static Collection<SingleInterval> parseCluster(String line) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		//NYBot10_Stg.Odd2Bo51.Even2Bo19.Odd2Bo29.Even2Bo90.Odd2Bo17.Even2Bo9.Odd2Bo86	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HEK-EZH2-16-chr1:90057393	HEK-EZH2-16-chr6:78716936
		
		String[] tokens=line.split("\t");
		boolean huh=false;
		for(int i=1; i<tokens.length; i++){
			if(tokens[i].startsWith("HUHNOSPLINT")){huh=true;}
			else{
				String[] region=tokens[i].split("-");
				String pos=region[region.length-1];
				int start=new Integer(pos.split(":")[1]);
				int end=start+1;
				SingleInterval r=new SingleInterval(pos.split(":")[0], start, end);
				rtrn.add(r);
			}
		}
		//if(huh){return rtrn;}
		return rtrn;
	}

	
	private static Collection<String> parseHUH(String line) {
		Collection<String> rtrn=new TreeSet<String>();
		//NYBot10_Stg.Odd2Bo51.Even2Bo19.Odd2Bo29.Even2Bo90.Odd2Bo17.Even2Bo9.Odd2Bo86	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HUHNOSPLINT	HEK-EZH2-16-chr1:90057393	HEK-EZH2-16-chr6:78716936
		
		String[] tokens=line.split("\t");
		
		for(int i=1; i<tokens.length; i++){
			if(tokens[i].startsWith("HUHNOSPLINT")){rtrn.add(tokens[i]);}
			
		}
		//if(huh){return rtrn;}
		return rtrn;
	}

	
	
	
	
	private static MatrixWithHeaders getContacts(Collection<Cluster> data, int binResolution) {
		List<String> rows=getRegions(data, binResolution);
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, rows);
		
		for(Cluster c: data){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					if(!region1.equals(region2)){
						rtrn.incrementCount(region1.toUCSC(), region2.toUCSC());
					}
				}
			}
		}
		
		
		return rtrn;
	}
	
	
	
	private static List<String> getRegions(Collection<Cluster> data, int resolution) {
		Collection<SingleInterval> list=new TreeSet<SingleInterval>();
		
		for(Cluster c: data){
			Cluster binned=c.bin(resolution);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				list.add(region);
			}
		}
		
		List<String> rtrn=new ArrayList<String>();
		for(SingleInterval region: list){
			rtrn.add(region.toUCSC());
		}
		
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException{
		if(args.length>2){
			Collection<Cluster> data=parse(args[0]);
			int binResolution=new Integer(args[1]);
			String save=args[2];
			
			
			MatrixWithHeaders contacts=getContacts(data, binResolution);
			contacts.write(save);
	
		}
		else{System.err.println(usage);}
	}

	

	static String usage=" args[0]=clusters \n args[1]=bin resolution \n args[2]=save";


	

	
}

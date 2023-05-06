package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;

public class QuantifyChromosome {

	static Collection<String> excludedList;
	static int minCount=9;
	
	public static void main(String[] args) throws IOException{
		excludedList=new ArrayList<String>();
		excludedList.add("Gm42418");
		excludedList.add("18S");
		excludedList.add("28S");
		excludedList.add("ITS1");
		excludedList.add("5.8S");
		excludedList.add("5S");
		
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String rna=args[1];
		int minClusterSize=new Integer(args[2]);
		int maxClusterSize=new Integer(args[3]);
		Map<String, Integer> map=new TreeMap<String, Integer>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.getClusterSize()>minClusterSize && c.getClusterSize()<maxClusterSize && c.containsRNA(rna)){ //TODO Count number of instances
				if(count(c, rna)>minCount){
				//if(!hasExcludedRNA(c, excludedList)){
					/*System.err.print(c.getBarcode());
					for(SingleInterval rnaRegion: c.getAllRNARegions()){
						System.err.print(" "+rnaRegion.getName()+"_"+rnaRegion.toUCSC());
					}
					System.err.println();*/
					//System.out.println(c.getAllDNAIntervals().size()+" "+c.getAllRNARegions().size());
					for(SingleInterval region: c.getAllDNAIntervals()){
						int count=0;
						if(map.containsKey(region.getReferenceName())){count=map.get(region.getReferenceName());}
						count++;
						map.put(region.getReferenceName(), count);
					}
				}
			//}
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
			
		}
		data.close();
		
		for(String chr: map.keySet()){
			System.out.println(chr+" "+map.get(chr));
		}
	}

	private static int count(Cluster c, String rna) {
		int count=0;
		for(SingleInterval rnaRegion: c.getAllRNARegions()){
			if(rnaRegion.getName().equals(rna)){count++;}
		}
		//System.err.println(c.getBarcode()+" "+rna+" "+count);
		return count;
	}

	private static boolean hasExcludedRNA(Cluster c, Collection<String> excludedList2) {
		for(String name: excludedList2){
			if(c.containsRNA(name)){return true;}
		}
		return false;
	}
}

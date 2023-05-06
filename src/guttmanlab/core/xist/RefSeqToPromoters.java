package guttmanlab.core.xist;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class RefSeqToPromoters {

	
	public static void main(String[] args) throws IOException {
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		FileWriter writer=new FileWriter(args[1]);
		int binResolution=Integer.parseInt(args[2]);
		
		
		Collection<SingleInterval> promoters=new TreeSet<SingleInterval>();
		for(Gene g: genes) {
			SingleInterval promoter=new SingleInterval(g.getReferenceName(), g.get5PrimePosition()-100, g.get5PrimePosition()+100);
			promoter.setOrientation(g.getOrientation());
			promoters.add(promoter);
			//writer.write(g.get5Prime(1000).toBED()+"\n");
		}
		
		Map<SingleInterval, Integer> counts=quantify(promoters, binResolution);
		
		for(SingleInterval p: counts.keySet()) {
			writer.write(p.toBedgraph(counts.get(p))+"\n");
		}
		
		
		writer.close();
	}

	private static Map<SingleInterval, Integer> quantify(Collection<SingleInterval> promoters, int binResolution) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(SingleInterval p: promoters) {
			SingleInterval binned=p.bin(binResolution);
			int count=0;
			if(rtrn.containsKey(binned)) {count=rtrn.get(binned);}
			count++;
			rtrn.put(binned,  count);
		}
		
		return rtrn;
	}
	
	
	
	
}

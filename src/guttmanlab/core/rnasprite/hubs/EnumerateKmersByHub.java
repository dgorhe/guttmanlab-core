package guttmanlab.core.rnasprite.hubs;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.IntegerDistributionAbstractTest;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.Kmer;

public class EnumerateKmersByHub {

	
	private static Collection<Kmer> parseHubs(String string) throws IOException {
		List<String> lines= BEDFileIO.loadLines(string);
		
		Map<String, Collection<String>> temp=new TreeMap<String, Collection<String>>();
		for(String line: lines) {
			String hub=line.split("\t")[0];
			String className=line.split("\t")[2];
			Collection<String> set=new TreeSet<String>();
			if(temp.containsKey(hub)) {set=temp.get(hub);}
			set.add(className);
			temp.put(hub, set);
		}
		
		
		Collection<Kmer> rtrn=new ArrayList<Kmer>();
		for(String hub: temp.keySet()) {
			Kmer kmer=new Kmer();
			kmer.setName(hub);
			kmer.addRegions(temp.get(hub));
			rtrn.add(kmer);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		Collection<Kmer> hubs=parseHubs(args[0]);
		String save=args[1];
		//int k=Integer.parseInt(args[2]);
		
		FileWriter writer=new FileWriter(save);
		/*(writer.write("kmer\thub\n");
		for(Kmer hub: hubs) {
			Collection<Kmer> kmers=hub.enumerateSubK(k);
			for(Kmer kmer: kmers) {writer.write(kmer.toString()+"\t"+hub.getName()+"\n");}
		*/
		
		for(Kmer hub: hubs) {
			writer.write(hub.toString()+"\t"+hub.getName()+"\n");
		}
		
		writer.close();
	}
	
}

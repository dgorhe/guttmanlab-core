package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.Kmer;

public class BatchSubmit {

	
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
	
	private static Collection<Kmer> enumerateSub(Kmer hub) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		rtrn.add(hub);
		for(int i=3; i<hub.getSize(); i++) {
			Collection<Kmer> kmers=hub.enumerateSubK(i);
			rtrn.addAll(kmers);
		}
		
		return rtrn;
	}
	
	
	
	private static void run(File barcodes, String allRNAs, String kmer, String save, String name, String probTree) throws IOException, InterruptedException {
		String cmd="java -jar -Xmx16000m /groups/guttman/mguttman/scripts/RandomDistOfHubs.jar "+barcodes.getAbsolutePath()+" "+probTree+" "+allRNAs+" "+kmer+" "+save;
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Process p=Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
		
	}

	
	private static Collection<Kmer> parseKmers(String string) throws IOException {
		Collection<Kmer> rtrn=new ArrayList<Kmer>();
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			Kmer k=new Kmer(tokens[0]);
			k.setName(tokens[0]);
			rtrn.add(k);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length>4) {
			Collection<Kmer> kmers=parseKmers(args[0]);
			File barcodeFile=new File(args[1]);
			
			
			String allRNAs=args[2];
			String saveBase=args[3];
			String probTree=args[4];
		
			int counter=0;
			for(Kmer k: kmers) {
				String name=k.getName();
				String save=saveBase+"/"+name+".scores";
				run(barcodeFile, allRNAs, name, save, name, probTree);
				counter++;
			}	
		}
		else {System.err.println(usage);}
	}

	

	static String usage=" args[0]=kmers \n args[1]=barcodes \n args[2]=all rnas \n args[3]=saveDir \n args[4]=prob tree";

	
	
}

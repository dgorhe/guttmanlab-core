package guttmanlab.core.rnasprite.hubs;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.Kmer;

public class RandomizeKmers {
	
	int randomSize=1000;

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
	
	private static void write(Collection<Kmer> randomHubs, String string) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Kmer k: randomHubs) {
			writer.write(k.toString()+"\t"+k.getName()+"\n");
		}
		
		writer.close();
	}
	
	private static Collection<String> getRNAs(Map<String, String> allRNAs, String region) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String name: allRNAs.keySet()) {
			String set=allRNAs.get(name);
			if(set.equals(region)) {rtrn.add(name);}
		}
		
		return rtrn;
	}

	private static Map<String, String> parseRNAs(String string) throws IOException {
		List<String> lines= BEDFileIO.loadLines(string);
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String line: lines) {
			String rna=line.split("\t")[1].replaceAll("\"", "");
			String className=line.split("\t")[2];
			rtrn.put(rna, className);
		}
		
		return rtrn;
	}
	
	private static List<String> getRegions(Collection<Kmer> hubs) {
		List<String> rtrn=new ArrayList<String>();
		
		for(Kmer hub: hubs) {
			rtrn.addAll(hub.getRegions());
		}
		return rtrn;
	}
	
	
	private static Collection<Kmer> randomizeHubs(Collection<Kmer> hubs, List<String> regions) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(Kmer hub: hubs) {
			int size=hub.getSize();
			Kmer rand=sample(regions, size, "randomized_"+hub.getName());
			rtrn.add(rand);
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
	
	private static Kmer sample(List<String> regions, int size, String name) {
		Kmer rtrn=new Kmer();
		rtrn.setName(name);
		
		for(int i=0; i<size; i++) {
			int index=new Double(Math.random()*regions.size()).intValue();
			rtrn.addRegion(regions.remove(index));
		}
		return rtrn;
	}

	public static Collection<Kmer> enumerateRandomKmers(Collection<Kmer> hubs, Collection<Kmer> possibleKmers){
		int maxKmer=maxKmer(hubs);
		List<String> allRegions=getRegions(hubs);
		
		
		Collection<Kmer> randomK=new TreeSet<Kmer>();
		for(int i=3; i<=maxKmer; i++) {
			randomK.addAll(random(allRegions, i, 100, possibleKmers));
		}
		
		return randomK;
	}
	
	public static Collection<Kmer> enumerateRandomKmers(Collection<Kmer> hubs, Collection<Kmer> possibleKmers, int k){
		Kmer all=new Kmer();
		for(Kmer hub: hubs) {all.addRegions(hub.getRegions());}
		Collection<Kmer> allKmers= all.enumerateSubK(k);
		
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		for(Kmer kmer: allKmers) {
			if(contains(possibleKmers,kmer)) {kmer.setName("actual");}
			else {
				kmer.setName("random");
				rtrn.add(kmer);
			}
		}
		
		return rtrn;
	}
	
	
	public static Collection<Kmer> enumerateRandomKmers(Collection<Kmer> hubs, int k){
		Kmer all=new Kmer();
		for(Kmer hub: hubs) {all.addRegions(hub.getRegions());}
		
		
		
		Collection<Kmer> allKmers= all.enumerateSubK(k);
		
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		for(Kmer kmer: allKmers) {
			kmer.setName("all");
			rtrn.add(kmer);
			
		}
		
		
		System.err.println(all.getSize()+" "+k+" "+rtrn.size());
		return rtrn;
	}
	
	private static Collection<Kmer> random(List<String> allRegions, int k, int numPerm, Collection<Kmer> possibleKmers) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(int i=0; i<numPerm; i++) {
			Kmer rand=pickK(allRegions, k);
			rand.setName("random");
			if(rand.getSize()==k && !contains(possibleKmers,rand)) {rtrn.add(rand);}
		}
		
		return rtrn;
	}

	private static boolean contains(Collection<Kmer> possibleKmers, Kmer rand) {
		for(Kmer k: possibleKmers) {
			if(rand.equals(k)) {return true;}
		}
		return false;
	}

	private static Kmer pickK(List<String> allRegions, int k) {
		Kmer kmer=new Kmer();
		for(int i=0; i<k; i++) {
			kmer.addRegion(allRegions.get(new Double(Math.random()*allRegions.size()).intValue()));
		}
		return kmer;
	}

	private static int maxKmer(Collection<Kmer> hubs) {
		int max=0;
		for(Kmer hub: hubs) {max=Math.max(max, hub.getSize());}
		return max;
	}

	private static Collection<Kmer> enumerateAll(Collection<Kmer> hubs) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		for(Kmer hub: hubs) {
			Collection<Kmer> subK=enumerateSub(hub);
			rtrn.addAll(subK);
		}
		return rtrn;
	}
	
	private static Collection<Kmer> enumerateAll(Collection<Kmer> hubs, int k) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		for(Kmer hub: hubs) {
			Collection<Kmer> subK=hub.enumerateSubK(k);
			rtrn.addAll(subK);
		}
		return rtrn;
	}

	
	private static Collection<Kmer> subset(Collection<Kmer> randomKmer, int num) {
		if(randomKmer.size()<num) {return randomKmer;}
		double proportion=(double)num/(double)randomKmer.size();
		
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(Kmer k: randomKmer) {
			double rand=Math.random();
			if(rand<proportion) {rtrn.add(k);}
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			int k=Integer.parseInt(args[2]);
			Collection<Kmer> hubs=parseHubs(args[0]);
			
			Collection<Kmer> randomKmer=enumerateRandomKmers(hubs, k);
			
			write(randomKmer, args[1]+"."+k+".kmers");
			
			//write(subset(randomKmer, 1000), args[1]+".random.subset");
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=hubs \n args[1]=save \n args[2]=k";
	

	

	

	
}

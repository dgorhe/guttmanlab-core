package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.rnasprite.Kmer;

public class ScoreKmers {

	public ScoreKmers(BarcodingDataStreaming data, Map<String, String> allRNA, String save) throws IOException {
		Map<Kmer, Integer> scores=new TreeMap<Kmer, Integer>();
		
		
		FileWriter writer=new FileWriter(save);
		int counter=0;
		int filtered=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Kmer full=getHits(c, allRNA);
			//System.err.println(full);
			
			//Collection<Kmer> subK=getSubK(full);
			if(full.getSize()>1) {
				writer.write(c+"\n");
				filtered++;
			}
			
			if (full.getSize()>2) {
				int count=0;
				if(scores.containsKey(full)) {count=scores.get(full);}
				count++;
				scores.put(full, count);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+filtered);}
		}
		
		data.close();
		
		writer.close();
		writeScores(save+".scores", scores);
		
		//Map<Kmer, Collection<Cluster>> clusters=getClusters(scores, allRNA, data);
		
		//write(save, clusters);
		
	}
	
	private void writeScores(String string, Map<Kmer, Integer> scores) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Kmer k: scores.keySet()) {
			writer.write(k+"\t"+k.getSize()+"\t"+scores.get(k)+"\n");
		}
		
		writer.close();
	}

	private void write(String saveDir, Map<Kmer, Collection<Cluster>> clusters) throws IOException {
		
		for(Kmer k: clusters.keySet()) {
			String save=saveDir+"/"+k;
			write(save, clusters.get(k));
		}
		
	}

	private void write(String save, Collection<Cluster> list) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Cluster c: list) {
			writer.write(c.toString()+"\n");
		}
		
		writer.close();
	}

	private Map<Kmer, Collection<Cluster>> getClusters(Map<Kmer, Integer> scores, Map<String, String> allRNA, BarcodingDataStreaming data) {
		Map<Kmer, Collection<Cluster>> rtrn=new TreeMap<Kmer, Collection<Cluster>>();
		Collection<Kmer> kmers=filter(scores);
		
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Kmer full=getHits(c, allRNA);
			Collection<Kmer> matches=hitMatches(full, kmers);
			for(Kmer k: matches) {
				if(!rtrn.containsKey(k)) {rtrn.put(k, new ArrayList<Cluster>());}
				Collection<Cluster> list=rtrn.get(k);
				list.add(c);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		data.close();
		return rtrn;
	}

	private Collection<Kmer> hitMatches(Kmer full, Collection<Kmer> kmers) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		for(Kmer k: kmers) {
			if(hitMatches(full, k)) {rtrn.add(k);}
		}
		
		return rtrn;
	}

	private Map<Kmer, Collection<Kmer>> getSubKs(Collection<Kmer> kmers) {
		Map<Kmer, Collection<Kmer>> rtrn=new TreeMap<Kmer, Collection<Kmer>>();
		
		for(Kmer k: kmers) {
			Collection<Kmer> list=new TreeSet<Kmer>();
			for(String region: k.getRegions()) {
				list.add(k.remove(region));
			}
			rtrn.put(k,  list);
		}
		
		return rtrn;
	}

	private Collection<Kmer> filter(Map<Kmer, Integer> scores) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(Kmer k: scores.keySet()) {
			int score=scores.get(k);
			if(score>1) {rtrn.add(k);}
		}
		
		return rtrn;
	}

	private boolean hitMatches(Kmer full, Kmer setK) {
		int size=setK.getSize()-1;
		
		int count=0;
		for(String region: setK.getRegions()) {
			if(full.containsRegion(region)) {count++;}
		}
		
		if(count>=size) {return true;}
		
		return false;
		
		
	}
	
	/*private Collection<Kmer> hitMatches(Kmer full, Collection<Kmer> kmers) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		for(Kmer k: kmers) {
			if(match(full, k)) {rtrn.add(k);}
		}
		return rtrn;
	}*/

	private boolean match(Kmer full, Kmer k) {
		for(String region: k.getRegionList()) {
			Kmer sub=k.remove(region);
			if(matchExact(full, sub)) {return true;}
		}
		return false;
	}

	private boolean matchExact(Kmer full, Kmer sub) {
		for(String region: sub.getRegions()) {
			if(!full.containsRegion(region)) {return false;}
		}
		return true;
	}

	/*private void write(String save, Map<Kmer, Integer> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Kmer k: scores.keySet()) {
			writer.write(k.toString()+"\t"+k.getSize()+"\t"+scores.get(k)+"\n");
		}
		
		writer.close();
	}*/

	private Collection<Kmer> getSubK(Kmer full) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(int i=3; i<=full.getSize(); i++) {
			rtrn.addAll(full.enumerateSubK(i));
		}
		
		return rtrn;
	}

	private void add(Map<Kmer, Integer> scores, Collection<Kmer> subK) {
		for(Kmer k: subK) {
			int count=0;
			if(scores.containsKey(k)) {count=scores.get(k);}
			count++;
			scores.put(k,  count);
		}
		
	}
	
	private void add(Map<Kmer, Integer> scores, Kmer k) {
		int count=0;
			if(scores.containsKey(k)) {count=scores.get(k);}
			count++;
			scores.put(k,  count);
		
		
	}

	private Kmer getHits(Cluster c, Map<String, String> allRNA) {
		Kmer rtrn=new Kmer();
		for(String name: c.getRNANames()) {
			if(allRNA.containsKey(name)) {rtrn.addRegion(allRNA.get(name));}
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
	
	public static void main(String[] args) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		Map<String, String> allRNA=parseRNAs(args[1]);
		String save=args[2];
		new ScoreKmers(data, allRNA, save);
	}
	
}

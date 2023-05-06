package guttmanlab.core.nanobody;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

public class ScoreClusters {

	public ScoreClusters(File clusters, File input, String save) throws IOException {
		//merge clusters
		Collection<Cluster> allClusters=parse(clusters);
		
		//score each cluster
		Map<Cluster, String> scores=scoreEach(allClusters, input);
	
		write(save, scores);
	}
	
	
	private Map<Cluster, String> scoreEach(Collection<Cluster> allClusters, File input) throws IOException {
		Collection<String> allProteins=new TreeSet<String>();
		int counter=0;
		for(Cluster c: allClusters) {
			allProteins.addAll(c.proteins);
			counter++;
			if(counter%100000==0) {System.err.println(counter+" "+allClusters.size());}
		}
		
		Map<String, Integer> scores=score(input, allProteins);
		
		Map<Cluster, String> rtrn=new TreeMap<Cluster, String>();
		for(Cluster c: allClusters) {
			int sum=score(c, scores);
			String seq=getMostAbundantSeq(c, scores);
			c.setScore(sum);
			rtrn.put(c, seq);
		}
		
		return rtrn;
	}

	
	private String getMostAbundantSeq(Cluster c, Map<String, Integer> scores) {
		String mostAbundant=c.proteins.iterator().next();
		int maxScore=0;
		
		for(String protein: c.proteins) {
			int score=get(scores, protein);
			if(score>maxScore) {
				mostAbundant=protein;
				maxScore=score;
			}
		}
		
		return mostAbundant;
	}


	private Map<String, Integer> score(File input, Collection<String> set) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String line;
		int counter=0;
		while ((line = reader.readLine()) != null) {
			String seq=line.split("\t")[0];
			if(!seq.contains("*") && set.contains(seq)) {
				int count=0;
				if(rtrn.containsKey(seq)) {
					count=rtrn.get(seq);
				}
				count++;
				rtrn.put(seq, count);
				
			}
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		reader.close();
		return rtrn;
	}
	

	private int score(Cluster cluster, Map<String, Integer> map) {
		int sum=0;
		for(String c: cluster.proteins) {
				sum+=get(map,c);
		}
		return sum;
	}


	private int get(Map<String, Integer> map, String c) {
		int rtrn=0;
		if(map.containsKey(c)) {rtrn=map.get(c);}
		return rtrn;
	}


	private void write(String save, Map<Cluster, String> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Cluster c: scores.keySet()) {
			String name=c.getName();
			String seq=scores.get(c);
			writer.write(name+"\t"+seq+"\t"+c.getScore()+"\n");
		}
		
		writer.close();
	}


	private Collection<Cluster> parseAll(File[] clusters) throws IOException {
		Collection<Cluster> allClusters=new TreeSet<Cluster>();
		
		for(int i=0; i<clusters.length; i++) {
			System.err.println(clusters[i].getAbsolutePath());
			Collection<Cluster> list=parse(clusters[i]);
			allClusters.addAll(list);
		}
		
		return allClusters;
	}


	private Collection<Cluster> parse(File file) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String line;
		int counter=0;
		while ((line = reader.readLine()) != null) {
			String[] tokens=line.split("\t");
			Collection<String> list=new TreeSet<String>();
			for(int i=3; i<tokens.length; i++) {
				list.add(tokens[i]);
			}
			Cluster c=new Cluster(list);
			c.setName(tokens[0]);
			c.setScore(Integer.parseInt(tokens[2]));
			rtrn.add(c);
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		reader.close();
		return rtrn;
	}


	private Collection<Cluster> mergeClusters(Collection<Cluster> clusters) throws IOException {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		Map<String, Cluster> map=makeMap(clusters);
		
		for(Cluster c: clusters) {
			Cluster mergedClster=getAllClusters(c, map);
			rtrn.add(mergedClster);
		}
		
		
		
		return rtrn;
	}
	
	private Map<String, Cluster> makeMap(Collection<Cluster> clusters) {
		Map<String, Cluster> rtrn=new TreeMap<String, Cluster>();
		
		for(Cluster c: clusters) {
			for(String protein: c.proteins) {
				Collection<String> list=new TreeSet<String>();
				list.addAll(c.proteins);
				if(rtrn.containsKey(protein)) {
					list.addAll(rtrn.get(protein).proteins);
				}
				rtrn.put(protein, new Cluster(list));
			}
		}
		
		return rtrn;
	}

	private Cluster getAllClusters(Cluster c, Map<String, Cluster> map) {
		Collection<String> list=new TreeSet<String>();
		
		for(String p: c.proteins) {
			Cluster c1=map.get(p);
			list.addAll(c1.proteins);
		}
		
		
		return new Cluster(list);
	}
	
	public static void main(String[] args) throws IOException {
		File clusters=new File(args[0]);
		File input=new File(args[1]);
		String save=args[2];
		new ScoreClusters(clusters, input, save);
	}
	
}

package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.barcoding.analysis.Cluster;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;

public class InterchromosomalContacts {

	int numPerms=1000;
	Map<String, Integer> chrSizes=CoordinateSpace.MM10.getRefSizes();
	int binResolution=1000000;
	
	public InterchromosomalContacts(File[] files, Collection<SingleInterval> regions, String save) throws IOException{
		
		Map<Cluster, Collection<Cluster>> perms=new TreeMap<Cluster, Collection<Cluster>>();
		Collection<Cluster> allClusters=new TreeSet<Cluster>();
		
		FileWriter writer=new FileWriter(save);
		Collection<Cluster> kmers=generatePairs(regions);
		allClusters.addAll(kmers);
		
		for(Cluster kmer: kmers){
			System.err.println(kmer.toString());
			List<Cluster> random=getPerms(kmer, numPerms);
			perms.put(kmer, random);
			allClusters.addAll(random);
		}
		
		Map<Cluster, Double> scores=score(files, allClusters);
		
		for(Cluster kmer: perms.keySet()){
			double observed=scores.get(kmer);
			double[] randomScores=getScores(scores, perms.get(kmer));
			//double max=Statistics.max(randomScores);
			double percentile=Statistics.percentLessThan(observed, randomScores);
			double average=Statistics.mean(randomScores);
			writer.write(kmer.toString(false)+"\t"+observed+"\t"+average+"\t"+percentile+"\n");
			//for(int i=0; i<randomScores.length; i++){writer.write("\t"+randomScores[i]);}
			//writer.write("\n");
		}
		
		writer.close();
	}

	
	private double[] getScores(Map<Cluster, Double> scores, Collection<Cluster> kmers) {
		double[] rtrn=new double[kmers.size()];
		
		int i=0;
		for(Cluster c: kmers){
			double score=0;
			if(scores.containsKey(c)){
				score=scores.get(c);
			}
			rtrn[i]=score;
			i++;
		}
		
		return rtrn;
	}


	private Map<Cluster, Double> score(File[] files, Collection<Cluster> allClusters) throws IOException{
		Map<Cluster, Double> num=new TreeMap<Cluster, Double>();
		Map<Cluster, Double> denom=new TreeMap<Cluster, Double>();
		
		for(int i=0; i<files.length; i++){
			if(i%10==0){System.err.println(i+" "+files.length);}
			MatrixWithHeaders mwh=new MatrixWithHeaders(files[i]);
			
			for(Cluster c: allClusters){
				if(hasKmer(mwh, c)){increment(num, c);}
				increment(denom, c);
			}
		}
		
		System.err.println(num.size()+" "+denom.size());
		
		Map<Cluster, Double> rtrn=new TreeMap<Cluster, Double>();
		for(Cluster c: num.keySet()){
			double numerator=num.get(c);
			double denominator=denom.get(c);
			double fraction=numerator/denominator;
			rtrn.put(c, fraction);
		}
		
		return rtrn;
	}
	
	private void increment(Map<Cluster, Double> num, Cluster c) {
		double val=0;
		if(num.containsKey(c)){val=num.get(c);}
		val=val+1;
		num.put(c, val);
	}


	private double[] fractionOfCells(File[] files, List<Cluster> random, Cluster observed) throws IOException {
		double[] counter=new double[random.size()+1];
		double[] total=new double[random.size()+1];
		
		for(int i=0; i<files.length; i++){
			if(i%10==0){System.err.println(i+" "+files.length);}
			MatrixWithHeaders mwh=new MatrixWithHeaders(files[i]);
			
			if(hasKmer(mwh, observed)){counter[0]++;}
			total[0]++;
			
			for(int j=0; j<random.size(); j++){
				if(hasKmer(mwh, random.get(j))){counter[j+1]++;}
				total[j+1]++;
			}
		}
		
		double[] rtrn=new double[random.size()+1];
		for(int j=0; j<rtrn.length; j++){
			rtrn[j]=counter[j]/total[j];
		}
		return rtrn;
	}

	private double fractionOfCells(File[] files, Cluster kmer) throws IOException {
		double counter=0;
		double total=0;
		for(int i=0; i<files.length; i++){
			MatrixWithHeaders mwh=new MatrixWithHeaders(files[i]);
			if(hasKmer(mwh, kmer)){counter++;}
			total++;
		}
		
		return counter/total;
	}

	private boolean hasKmer(MatrixWithHeaders mwh, Cluster kmer) {
		Iterator<SingleInterval> iter=kmer.getAllIntervals().iterator();
		SingleInterval region1=iter.next();
		SingleInterval region2=iter.next();
		
		Collection<String> bins1=getBins(region1, this.binResolution);
		Collection<String> bins2=getBins(region2, this.binResolution);
		
		for(String bin1: bins1){
			for(String bin2: bins2){
				if(mwh.containsColumn(bin1) && mwh.containsColumn(bin2)){
					double score=mwh.get(bin1, bin2);
					if(score>0){return true;}
				}
			}
		}
		return false;
	}

	private Collection<String> getBins(SingleInterval region1, int binResolution2) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(int i=region1.getReferenceStartPosition(); i<region1.getReferenceEndPosition(); i+=binResolution2){
			SingleInterval region=new SingleInterval(region1.getReferenceName(), i, i+binResolution2);
			rtrn.add(region.toUCSC());
			//System.err.println(region1.toUCSC()+" "+region.toUCSC());
		}
		return rtrn;
	}

	private Collection<Cluster> generatePairs(Collection<SingleInterval> regions) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(SingleInterval region1: regions){
			for(SingleInterval region2: regions){
				if(!region1.getReferenceName().equals(region2.getReferenceName())){
					Cluster c=new Cluster("pairs");
					c.addRead(region1);
					c.addRead(region2);
					rtrn.add(c);
				}
			}
		}
		return rtrn;
	}

	private List<Cluster> getPerms(Cluster kmer, int numPerm) {
		List<Cluster> rtrn=new ArrayList<Cluster>();
		for(int i=0; i<numPerm; i++){
			Cluster randomCluster=kmer.getPermutedCluster(chrSizes, binResolution);
			rtrn.add(randomCluster);
		}
		return rtrn;
	}
	
	
	private static Map<String, MatrixWithHeaders> parse(File[] listFiles) throws IOException {
		Map<String, MatrixWithHeaders> rtrn=new TreeMap<String, MatrixWithHeaders>();
		
		for(int i=0; i<listFiles.length; i++){ //TODO Fix this
			if(i%10==0){System.err.println(i+" "+listFiles.length);}
			rtrn.put(listFiles[i].getName(), new MatrixWithHeaders(listFiles[i]));
		}
		System.err.println("Done loading matrices");
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		//Map<String, MatrixWithHeaders> matrices=parse(new File(args[0]).listFiles());
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		
		//Nucleolus
		/*regions.add(new SingleInterval("chr12:5000000-17000000"));
		regions.add(new SingleInterval("chr12:25000000-32000000"));
		regions.add(new SingleInterval("chr15:3000000-6000000"));
		regions.add(new SingleInterval("chr15:67000000-71000000"));
		regions.add(new SingleInterval("chr16:5000000-8000000"));
		regions.add(new SingleInterval("chr18:3000000-10000000"));
		regions.add(new SingleInterval("chr18:13000000-24000000"));
		regions.add(new SingleInterval("chr18:39000000-42000000"));
		regions.add(new SingleInterval("chr18:57000000-60000000"));
		regions.add(new SingleInterval("chr19:11000000-24000000"));
		regions.add(new SingleInterval("chr19:25000000-28000000"));
		regions.add(new SingleInterval("chr19:29000000-37000000"));
		regions.add(new SingleInterval("chr19:48000000-53000000"));
		regions.add(new SingleInterval("chr19:58000000-61000000"));*/
		
		//Speckle
		/*regions.add(new SingleInterval("chr2:164000000-174000000"));
		regions.add(new SingleInterval("chr2:177000000-181000000"));
		regions.add(new SingleInterval("chr4:128000000-142000000"));
		regions.add(new SingleInterval("chr4:147000000-155000000"));
		regions.add(new SingleInterval("chr5:112000000-126000000"));
		regions.add(new SingleInterval("chr8:123000000-127000000"));
		regions.add(new SingleInterval("chr11:95000000-103000000"));
		regions.add(new SingleInterval("chr11:115000000-121000000"));
		regions.add(new SingleInterval("chr13:55000000-58000000"));
		regions.add(new SingleInterval("chr15:76000000-79000000"));
		regions.add(new SingleInterval("chr17:25000000-30000000"));*/
		
		//Centromeres
		regions.add(new SingleInterval("chr1:3000000-13000000"));
		regions.add(new SingleInterval("chr2:3000000-13000000"));
		regions.add(new SingleInterval("chr3:3000000-13000000"));
		regions.add(new SingleInterval("chr4:3000000-13000000"));
		regions.add(new SingleInterval("chr5:3000000-13000000"));
		regions.add(new SingleInterval("chr6:3000000-13000000"));
		regions.add(new SingleInterval("chr7:3000000-13000000"));
		regions.add(new SingleInterval("chr8:3000000-13000000"));
		regions.add(new SingleInterval("chr9:3000000-13000000"));
		regions.add(new SingleInterval("chr10:3000000-13000000"));
		regions.add(new SingleInterval("chr11:3000000-13000000"));
		regions.add(new SingleInterval("chr12:3000000-13000000"));
		regions.add(new SingleInterval("chr13:3000000-13000000"));
		regions.add(new SingleInterval("chr14:3000000-13000000"));
		regions.add(new SingleInterval("chr15:3000000-13000000"));
		regions.add(new SingleInterval("chr16:3000000-13000000"));
		regions.add(new SingleInterval("chr17:3000000-13000000"));
		regions.add(new SingleInterval("chr18:3000000-13000000"));
		regions.add(new SingleInterval("chr19:3000000-13000000"));
		regions.add(new SingleInterval("chrX:3000000-13000000"));
		
		new InterchromosomalContacts(files, regions, save);
		
		
		
	}
	
}

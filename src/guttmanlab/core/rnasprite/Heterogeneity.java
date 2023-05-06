package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.datastructures.MatrixWithHeaders;

public class Heterogeneity {

	
	BarcodingDataStreaming data;
	Map<String, Integer> countsPerGene;
	int clusterSize;
	double foldCutoff=1;
	int minCount=100;
	double freqCutoff;
	
	public Heterogeneity(BarcodingDataStreaming allData, String rna, double cutoff, String save) throws IOException{
		
		this.freqCutoff=cutoff;
		getCounts(allData, rna);
		
		int minCount=new Double(cutoff*clusterSize).intValue();
		List<String> genes=getAllGenes(data, minCount);
		genes=filterPairwise(data, countsPerGene, genes, rna);
		
		System.err.println("Genes "+genes.size());
		
		Map<Kmer, Integer> kmerCounts=getKmerCounts(genes); //TODO Write the clusters directly --> MWH of clusters
		Map<Kmer, Integer> collapsed=collapse(kmerCounts);
		
		Map<Kmer, Double> enrichment=computeEnrichment(collapsed);
		
		//write(save+".unique", kmerCounts);
		//write(save+".cumulative", collapsed);
		write(save+".enrichment", enrichment, collapsed);
	}
	
	private Map<Kmer, Double> computeEnrichment(Map<Kmer, Integer> collapsed) {
			/**
		 	double geneCount=countsPerGene.get(gene);
			double expected=rnaCount*(geneCount/total);
			double pairwiseCount=pairwiseCounts.get(gene);
			double ratio=pairwiseCount/expected;
			**/
			
		Map<Kmer, Double> rtrn=new TreeMap<Kmer, Double>();
		
		double totalCount=this.countsPerGene.get("Total");
		
		for(Kmer kmer: collapsed.keySet()){
			double kmerCount=collapsed.get(kmer);
			double minEnrichment=Double.MAX_VALUE;
			//for each region in kmer, exclude 1 and make subK
			for(String region: kmer.getRegions()){
				Kmer subK=kmer.remove(region);
				double subKCount=0.0;
				if(collapsed.containsKey(subK)){
					subKCount=collapsed.get(subK);
				
					//else{subKCount=kmerCount;}
					double singleCount=this.countsPerGene.get(region);
					double expected=subKCount*(singleCount/totalCount);
					double enrichment=kmerCount/expected;
					//System.err.println(kmer.toString()+" "+region+" "+kmerCount+" "+subKCount+" "+singleCount+" "+totalCount+" "+enrichment);
					minEnrichment=Math.min(enrichment, minEnrichment);
				}
			}
			rtrn.put(kmer, minEnrichment);
		}
		
		return rtrn;
	}

	private void write(String save, Map<Kmer, ? extends Number> kmerCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		
		for(Kmer kmer: kmerCounts.keySet()){
			Number count=kmerCounts.get(kmer);
			writer.write(kmer.toString()+"\t"+count+"\n");
		}
		
		writer.close();
	}
	
	private void write(String save, Map<Kmer, ? extends Number> enrichment, Map<Kmer, ? extends Number> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		
		for(Kmer kmer: enrichment.keySet()){
			Number e=enrichment.get(kmer);
			Number count=counts.get(kmer);
			writer.write(kmer.toString()+"\t"+count+"\t"+e+"\n");
		}
		
		writer.close();
	}
	
	private List<String> filterPairwise(BarcodingDataStreaming data, Map<String, Integer> countsPerGene, List<String> genes, String rna) {
		List<String> rtrn=new ArrayList<String>();
		
		double total=countsPerGene.get("Total");
		double rnaCount=countsPerGene.get(rna);
		
		Map<String, Integer> pairwiseCounts=new TreeMap<String, Integer>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			for(String gene: c.getRNANames()){
				int count=0;
				if(pairwiseCounts.containsKey(gene)){count=pairwiseCounts.get(gene);}
				count++;
				pairwiseCounts.put(gene, count);
			}
			counter++;
			if(counter%10000==0){System.err.println(counter);}
		}
		
		data.close();
		
		for(String gene: genes){
			double geneCount=countsPerGene.get(gene);
			double expected=rnaCount*(geneCount/total);
			//double pairwiseCount=data.getClustersOverlappingRNA(gene).size();
			double pairwiseCount=pairwiseCounts.get(gene);
			double ratio=pairwiseCount/expected;
			System.err.println(gene+" "+pairwiseCount+" "+expected+" "+ratio);
			if(ratio>foldCutoff){rtrn.add(gene);}
		}
		
		return rtrn;
	}
	
	
	private List<String> getAllGenes(BarcodingDataStreaming data, int minCount) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		int counter=0;
		
		while(data.hasNext()){
			Cluster c=data.next();
			for(String gene: c.getRNANames()){
				int count=0;
				if(rtrn.containsKey(gene)){count=rtrn.get(gene);}
				count++;
				rtrn.put(gene, count);
			}
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		
		data.close();
		
		List<String> list=new ArrayList<String>();
		for(String gene: rtrn.keySet()){
			if(!gene.equals("Unassigned")){
				int count=rtrn.get(gene);
				if(count>minCount){list.add(gene);}
			}
		}
		
		return list;
	}
	
	
	private void getCounts(BarcodingDataStreaming allData, String rna) throws IOException {
		String save=allData.getBarcodeFile().getAbsolutePath()+"."+rna+".clusters";
		FileWriter writer=new FileWriter(save);
		
		this.countsPerGene=new TreeMap<String, Integer>();
		int counter=0;
		int retained=0;
		while(allData.hasNext()){
			Cluster c=allData.next();
			if(c.getRNANames().contains(rna)){
				writer.write(c.toString()+"\n"); 
				retained++;
			}
			for(String gene: c.getRNANames()){
				int count=0;
				if(countsPerGene.containsKey(gene)){count=countsPerGene.get(gene);}
				count++;
				countsPerGene.put(gene, count);
			}
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
		
		this.countsPerGene.put("Total", counter);
		writer.close();
		allData.close();
		System.err.println("Retained "+retained);
		this.clusterSize=retained;
		this.data=new BarcodingDataStreaming(new File(save));
	}
	
	/*private Map<Kmer, Integer> getKmerCounts(List<String> genes, String save) throws IOException {
		Map<Kmer, Collection<Cluster>> kmerCounts=new TreeMap<Kmer, Collection<Cluster>>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			Collection<Cluster> count=new ArrayList<Cluster>();
			Kmer kmer=getKmer(c, genes);
			if(kmerCounts.containsKey(kmer)){count=kmerCounts.get(kmer);}
			count.add(c);
			kmerCounts.put(kmer, count);
		}
		data.close();
		
		Map<Kmer, Integer> rtrn=new TreeMap<Kmer, Integer>();
		for(Kmer kmer: kmerCounts.keySet()){
			System.err.println(kmer.toString());
			Collection<Cluster> clusters=kmerCounts.get(kmer);
			if(clusters.size()>minCount){
				MatrixWithHeaders mwh=makeClusterMatrix(clusters);
				mwh.write(save+"."+kmer.toFileName());
			}
			rtrn.put(kmer, clusters.size());
		}
		
		return rtrn;
	}*/
	
	private Map<Kmer, Integer> getKmerCounts(List<String> genes) {
		Map<Kmer, Integer> kmerCounts=new TreeMap<Kmer, Integer>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			int count=0;
			Kmer kmer=getKmer(c, genes);
			if(kmerCounts.containsKey(kmer)){count=kmerCounts.get(kmer);}
			count++;
			kmerCounts.put(kmer, count);
		}
		
		
		data.close();
		return kmerCounts;
	}

	private MatrixWithHeaders makeClusterMatrix(Collection<Cluster> clusters) {
		List<String> rows=new ArrayList<String>();
		List<String> columns=new ArrayList<String>();
		
		Map<String, Integer> geneCount=new TreeMap<String, Integer>();
		
		for(Cluster c: clusters){
			rows.add(c.getBarcode());
			for(String gene:c.getRNANames()){
				int count=0;
				if(geneCount.containsKey(gene)){
					count=geneCount.get(gene);
				}
				count++;
				geneCount.put(gene, count);
			}
		}
		
		double cutoff=this.freqCutoff*clusters.size();
		for(String gene: geneCount.keySet()){
			int count=geneCount.get(gene);
			if(!gene.equals("Unassigned") && count>cutoff){
				double ratio=(double)count/cutoff;
				if(ratio>5){System.err.println(gene+" "+cutoff+" "+count+" "+ratio);}
				columns.add(gene);
			}
		}
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		for(Cluster c: clusters){
			for(String gene: columns){
				if(c.containsRNA(gene)){mwh.set(c.getBarcode(), gene, 1.0);}
			}
		}
		
		
		return mwh;
	}

	private Kmer getKmer(Cluster c, List<String> genes) {
		Kmer newCluster=new Kmer();
		for(String rna: c.getRNANames()){
			if(genes.contains(rna)){
				newCluster.addRegion(rna);
			}
		}
		return newCluster;
	}
	
	private Map<Kmer, Integer> collapse(Map<Kmer, Integer> kmerCounts) {
		Map<Kmer, Integer> rtrn=new TreeMap<Kmer, Integer>();
		
		for(Kmer kmer: kmerCounts.keySet()){
			int count=getAllSubsets(kmerCounts, kmer);
			rtrn.put(kmer, count);
		}
		
		return rtrn;
	}
	
	private int getAllSubsets(Map<Kmer, Integer> kmerCounts, Kmer kmer) {
		int count=0;
		for(Kmer kmer2: kmerCounts.keySet()){
			if(kmer.isSubset(kmer2)){
				count+=kmerCounts.get(kmer2);
				//System.out.println(kmer.toString()+" "+kmer2.toString()+" "+kmerCounts.get(kmer2)+" "+count);
			}
		}
		return count;
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String rna=args[1];
			double cutoff=new Double(args[2]);
			String save=args[3];
			new Heterogeneity(data, rna,cutoff, save);
		}
		else{
			System.err.println(usage);
		}
	}

	static String usage=" args[0]=cluster file \n args[1]=rna \n args[2]=fraction cutoff \n args[3]=save";
	
}

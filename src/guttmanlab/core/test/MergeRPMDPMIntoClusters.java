package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.rnasprite.RNAInterval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class MergeRPMDPMIntoClusters {
	int binResolution;

	public MergeRPMDPMIntoClusters(File file, int minClusterSize, int maxClusterSize, Collection<Gene> genes, int binResolution,  String save, boolean makeBedgraph) throws IOException{
		this.binResolution=binResolution;
		SAMFileReader reader=new SAMFileReader(file);
		
		Map<String, Cluster> clusters=new TreeMap<String, Cluster>();
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			String name=read.getReadName();
			parse(name, read, clusters);
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		reads.close();
		reader.close();
	
	
		int count=0;
		int dna=0;
		int rna=0;
		int both=0;
		for(String barcode: clusters.keySet()){
			Cluster c=clusters.get(barcode);
				if(c.getClusterSize()>minClusterSize){
				if(c.hasDNA()){dna++;}
				if(c.hasRNA()){rna++;}
				if(c.hasDNA() && c.hasRNA()){both++;}
				count++;
			}
		}
		
		MatrixWithHeaders mwh=getDNAContacts(clusters, minClusterSize, maxClusterSize);
		mwh.write(save+".dnacontacts");
		
		double ratio=(double)both/(double)count;
		double ratio2=(double)both/(double)rna;
		System.err.println("DNA:"+dna+" RNA:"+rna+" RNA-DNA:"+both+" Total:"+count+" Fraction RNA-DNA out of total:"+ratio+" Fraction RNA-DNA out of all RNA:"+ratio2);
		
		if(!makeBedgraph){
			FileWriter writer=new FileWriter(save);
			for(Gene gene: genes){
				int numClusters=this.numClusters(clusters, gene.getGenomicRegion(), minClusterSize, maxClusterSize);
				writer.write(gene.getName()+"\t"+numClusters+"\n");
				System.err.println(gene.getName()+" "+numClusters);
			}
			writer.close();
		}
		
		if(makeBedgraph){
			for(Gene gene: genes){
				System.err.println(gene.getName());
				String saveToUse=save+"."+gene.getName()+".bedgraph";
				Map<SingleInterval, Integer> counts=getDNA(clusters, gene.getGenomicRegion(), minClusterSize,maxClusterSize);
				write(saveToUse, counts);
			}
		}
	}
	
	
	public MergeRPMDPMIntoClusters(File file, int minClusterSize, int maxClusterSize, int binResolution,  String save) throws IOException{
		this.binResolution=binResolution;
		SAMFileReader reader=new SAMFileReader(file);
		
		Map<String, Cluster> clusters=new TreeMap<String, Cluster>();
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			String name=read.getReadName();
			parse(name, read, clusters);
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		reads.close();
		reader.close();
	
		
		MatrixWithHeaders mwh=getDNAContacts(clusters, minClusterSize, maxClusterSize);
		mwh.write(save+".dnacontacts");
	}
	
	
	
	
	private MatrixWithHeaders getDNAContacts(Map<String, Cluster> clusters, int min, int max) {
		List<String> regions=getRegions(clusters);
		MatrixWithHeaders rtrn=new MatrixWithHeaders(regions, regions);
		
		for(Cluster c: clusters.values()){
			Cluster binned=c.bin(this.binResolution);
			if(c.getClusterSize()>min && c.getClusterSize()<max){
			for(SingleInterval r1:binned.getAllDNAIntervals()){
				for(SingleInterval r2: binned.getAllDNAIntervals()){
					if(!r1.equals(r2)){
						String row=r1.toUCSC();
						String column=r2.toUCSC();
						double count=rtrn.get(row, column);
						count++;
						rtrn.set(row, column, count);
					}
				}
			}
			}
		}
		return rtrn;
	}




	private List<String> getRegions(Map<String, Cluster> clusters) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(Cluster c: clusters.values()){
			Cluster binned=c.bin(this.binResolution);
			for(SingleInterval region: binned.getAllDNAIntervals()){
				rtrn.add(region);
			}
		}
		
		List<String> list=new ArrayList<String>();
		for(SingleInterval r: rtrn){list.add(r.toUCSC());}
		
		return list;
	}




	private Map<SingleInterval, Integer> getDNA(Map<String, Cluster> clusters, SingleInterval region, int min, int max) {
		List<Cluster> list=new ArrayList<Cluster>();
		for(Cluster c: clusters.values()){
			if(overlaps(c, region)){
				if(c.getClusterSize()>min && c.getClusterSize()<max){
					list.add(c);
				}
				//System.out.println(c.toString());
			}
		}
		
		System.err.println(list.size());
		
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		for(Cluster c: list){
			Cluster binned=c.bin(this.binResolution);
			for(SingleInterval r: binned.getAllDNAIntervals()){
				int count=0;
				if(rtrn.containsKey(r)){count=rtrn.get(r);}
				count++;
				rtrn.put(r, count);
			}
		}
		return rtrn;
	}
	
	private int numClusters(Map<String, Cluster> clusters, SingleInterval region, int min, int max) {
		List<Cluster> list=new ArrayList<Cluster>();
		for(Cluster c: clusters.values()){
			if(overlaps(c, region)){
				if(c.getClusterSize()>min && c.getClusterSize()<max){
					list.add(c);
				}
			}
		}
		
		return list.size();
	}


	private boolean overlaps(Cluster c, SingleInterval region) {
		for(RNAInterval r: c.getAllRNARegions()){
			if(overlaps(r, region)){return true;}
		}
		return false;
	}


	private boolean overlaps(RNAInterval r, SingleInterval region) {
		if(!r.getReferenceName().equals(region.getReferenceName())){return false;}
		int start=Math.max(r.getReferenceStartPosition(), region.getReferenceStartPosition());
		int end=Math.min(r.getReferenceEndPosition(), region.getReferenceEndPosition());
		if(start<end){return true;}
		return false;
	}


	private void write(String save, Map<SingleInterval, Integer> counts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: counts.keySet()){
			writer.write(region.toBedgraph(counts.get(region))+"\n");
		}
		
		writer.close();
	}


	private void parse(String name, SAMRecord read, Map<String, Cluster> clusters) {
		//FS10000829:26:BPA73104-2019:1:1101:10000:1010::[BIODPM][R8-3_TermStag_bot_3][R7-1_bot_A1][R6-1_bot_B3][R5-1_bot_A9][R4-1_bot_A7][R3-1_bot_A11][R2-1_bot_B5][R1-1_bot_A10]
		String barcode=name.split("::")[1];
		
		String[] tokens=barcode.split("\\[");
		boolean isDNA=isDNA(tokens[1]);
		barcode=getBarcode(tokens);
		//System.err.println(tokens[1]+ " "+barcode+" "+isDNA);
		
		
		Cluster cluster=new Cluster(barcode);
		if(clusters.containsKey(barcode)){
			cluster=clusters.get(barcode);
		}
		
		SAMFragment frag=new SAMFragment(read);
		SingleInterval region=frag.getGenomicInterval();
		if(isDNA){
			cluster.addDNARead(region);
		}
		else{
			RNAInterval r=new RNAInterval(region);
			cluster.addRNARead(r);
		}
		
		clusters.put(barcode, cluster);
		
	}


	private boolean isDNA(String string) {
		if(string.contains("DPM")){return true;}
		return false;
	}


	private String getBarcode(String[] tokens) {
		String rtrn="";
		for(int i=2; i<tokens.length; i++){
			String e=tokens[i].replaceAll("\\]", "_");
			rtrn+=e;
		}
		return rtrn;
	}


	public static void main(String[] args) throws IOException{
		/*if(args.length>5){
		File file=new File(args[0]);
		Collection<Gene> regions=BEDFileIO.loadRegionsFromFile(args[1]);
		int resolution=new Integer(args[2]);
		String save=args[3];
		boolean makeBedgraph=new Boolean(args[4]);
		int min=new Integer(args[5]);
		int max=new Integer(args[6]);
		
		new MergeRPMDPMIntoClusters(file, min, max, regions, resolution, save, makeBedgraph);
		}*/
		
		if(args.length>4){
			File file=new File(args[0]);
			int resolution=new Integer(args[1]);
			String save=args[2];
			int min=new Integer(args[3]);
			int max=new Integer(args[4]);
			
			new MergeRPMDPMIntoClusters(file, min, max, resolution, save);
			}
		
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=bin resolution \n args[2]=save \n args[3]=min cluster size \n args[4]=max cluster size";
	
}

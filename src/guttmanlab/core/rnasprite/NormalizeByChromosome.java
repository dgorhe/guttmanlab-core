package guttmanlab.core.rnasprite;

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
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;

public class NormalizeByChromosome {

	static SingleInterval problemBin1=new SingleInterval("chr2:79490000-79500000");
	static SingleInterval problemBin2=new SingleInterval("chr11:3119270-3192250");
	static SingleInterval problemBin3=new SingleInterval("chr15:99734977-99736026");
	static SingleInterval problemBin4=new SingleInterval("chr3:5173978-5175025");
	static SingleInterval problemBin5=new SingleInterval("chr13:58176952-58178051");
	
	static Collection<SingleInterval> problemBins=new TreeSet<SingleInterval>();

	
	private static void write(String save, Map<SingleInterval, Double> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			double val=map.get(region);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+val+"\n");
		}
		
		writer.close();
	}
	
	
	
	private static Collection<SingleInterval> writeSmoothed(String save, Map<SingleInterval, Double> map, int binResolution, int n, double cutoff) throws IOException{
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			Collection<SingleInterval> regions=getNeighbors(region, binResolution, n);
			double val=average(regions, map);
			if(val>cutoff){
				rtrn.add(region);
			}
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+val+"\n");
		}
		writer.close();
		return rtrn;
	}
	
	private static Collection<SingleInterval> getNeighbors(SingleInterval region, int binResolution, int n) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		//upstream
		for(int i=0; i<n; i++){
			int start=region.getReferenceStartPosition()-((i+1)*binResolution);
			int end=start+binResolution;
			SingleInterval newRegion=new SingleInterval(region.getReferenceName(), start, end);
			rtrn.add(newRegion);
		}
		
		
		//downstream
		for(int i=0; i<n; i++){
			int start=region.getReferenceStartPosition()+((i+1)*binResolution);
			int end=start+binResolution;
			SingleInterval newRegion=new SingleInterval(region.getReferenceName(), start, end);
			rtrn.add(newRegion);
		}
		
		return rtrn;
	}



	private static double average(Collection<SingleInterval> regions, Map<SingleInterval, Double> map) {
		double sum=0;
		double count=0;
		
		for(SingleInterval region: regions){
			if(map.containsKey(region)){
				sum+=map.get(region);
				count++;
			}
		}
		
		return sum/count;
	}

	private static Map<SingleInterval, Double> removeProblemBins(Map<SingleInterval, Double> map) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval region: map.keySet()){
			if(!overlaps(region, problemBins)){
				rtrn.put(region, map.get(region));
			}
		}
		
		return rtrn;
	}


	private static boolean overlaps(SingleInterval region, Collection<SingleInterval> problemBins) {
		for(SingleInterval problemBin: problemBins){
			if(region.overlaps(problemBin)){return true;}
		}
		return false;
	}

	
	private static Map<String, String> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(String line: lines){
			rtrn.put(line.split("\t")[0].replaceAll("\"", ""), line.split("\t")[1]);
		}
		return rtrn;
	}

	
	private static Collection<String> getClasses(Map<String, String> genes) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String rna: genes.keySet()){
			String rnaClass=genes.get(rna);
			rtrn.add(rnaClass);
		}
		
		return rtrn;
	}
	
	private static Kmer getKmer(Map<String, String> genes, String rnaClass) {
		Kmer rtrn=new Kmer();
		for(String rna: genes.keySet()){
			String rnaClass2=genes.get(rna);
			if(rnaClass2.equals(rnaClass)){rtrn.addRegion(rna);}
		}
		return rtrn;
	}


	private static Map<SingleInterval, Double> normalize(Map<SingleInterval, Double> map, String chr) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		Collection<Double> vals=new ArrayList<Double>();
		for(SingleInterval region: map.keySet()){
			if(region.getReferenceName().equals(chr)){vals.add(map.get(region));}
		}
		
		double avg=Statistics.mean(vals);
		
		for(SingleInterval region: map.keySet()){
			if(region.getReferenceName().equals(chr)){
				double val=map.get(region);
				double norm=val/avg;
				rtrn.put(region, norm);
			}
		}
		
		return rtrn;
	}


	private static void writeRegions(String save, Collection <SingleInterval> regions) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: regions){
			writer.write(region.toBED()+"\n");
		}
		
		writer.close();
	}
	
	
	private static Map<SingleInterval, Double> getRowScoreMap(MatrixWithHeaders mwh, String rnaName) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(String column: mwh.getColumnNames()){
			double val=mwh.get(rnaName, column);
			SingleInterval region=new SingleInterval(column);
			rtrn.put(region, val);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		problemBins.add(problemBin1);
		problemBins.add(problemBin2);
		problemBins.add(problemBin3);
		problemBins.add(problemBin4);
		problemBins.add(problemBin5);
		
		if(args.length>4){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Collection<SingleInterval> rnas=BEDFileIO.loadSingleIntervalFromFile(args[1]);
			int binResolution=new Integer(args[2]);
			String save=args[3];
			double cutoff=new Double(args[4]);
			
			//MatrixWithHeaders mwh=data.getRNADNAContacts(binResolution, rnas);
			
			Collection<String> genes=new TreeSet<String>();
			for(SingleInterval rna: rnas){
				String rnaName=rna.getName();
				if(!genes.contains(rnaName)){
					System.err.println(rnaName+" "+rna.getReferenceName());
					String chr=rna.getReferenceName();
					Map<SingleInterval, Double> map=data.getRNADNAContacts(binResolution, rnaName);
					map=normalize(map, chr);
					map=removeProblemBins(map);
					write(save+"."+rnaName+".bedgraph", map);
					Collection<SingleInterval> regions=writeSmoothed(save+"."+rnaName+".smoothed.bedgraph", map, binResolution, 5, cutoff);
					writeRegions(save+"."+rnaName+".bed", regions);
					//data.close();
					genes.add(rnaName);
				}
			}
		
			
			
		}
		else{System.err.println(usage);}
	}


	


	static String usage=" args[0]=clusters \n args[1]=RNA list (BED file) \n args[2]=bin resolution \n args[3]=save \n args[4]=fold enrichment over chromosome";
	
}

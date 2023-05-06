package guttmanlab.core.spidr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class ValuePerTSS {
	
	private int n=100;

	public ValuePerTSS(Collection<Gene> genes, File input, int index, String save) throws IOException {
		Collection<SingleInterval> tssPositions=getTSSPositions(genes, n);
		
		Map<SingleInterval, Double> scores=new TreeMap<SingleInterval, Double>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		String header;
		String totals;
		
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			if(counter==0){
				System.err.println(tokens[index]);
				header=nextLine;
			}
			else if(counter==1) {
				totals=nextLine;
			}
			else{
				SingleInterval region=new SingleInterval(tokens[0]);
				if(tssPositions.contains(region)) {
					double val=Double.parseDouble(tokens[index]);
					scores.put(region, val);
				}
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}	
		}
		reader.close();
		
		
		Map<Gene, List<Double>> scoreByPosition=getScoreByPosition(genes, scores, n);
		write(scoreByPosition, save);
	}

	private Map<Gene, List<Double>> getScoreByPosition(Collection<Gene> genes, Map<SingleInterval, Double> scores, int n2) {
		Map<Gene, List<Double>> rtrn=new TreeMap<Gene, List<Double>>();
		
		for(Gene gene: genes) {
			int tss=gene.get5PrimePosition();
			List<Double> list=new ArrayList<Double>();
			for(int i=-n2; i<n2; i++) {
				//TODO check strand
				SingleInterval r=new SingleInterval(gene.getReferenceName(), tss+i, tss+i+1);
				if(gene.getOrientation().equals(Strand.NEGATIVE)) {
					r=new SingleInterval(gene.getReferenceName(), tss-i, (tss-i)+1);
				}
				
				r.setOrientation(Strand.antisense(gene.getOrientation())); //TODO might need to take inverse strand
				double val=get(r, scores);
				list.add(val);
			}
			rtrn.put(gene, list);
		}
		
		return rtrn;
	}

	private double get(SingleInterval r, Map<SingleInterval, Double> scores) {
		if(scores.containsKey(r)) {return scores.get(r);}
		return 0;
	}

	private Collection<SingleInterval> getTSSPositions(Collection<Gene> genes, int n2) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(Gene gene: genes) {
			int tss=gene.get5PrimePosition();
			
			for(int i=-n2; i<n2; i++) {
				SingleInterval r=new SingleInterval(gene.getReferenceName(), tss+i, tss+i+1);
				r.setOrientation(Strand.antisense(gene.getOrientation()));
				System.out.println(r.toBED());
				rtrn.add(r);
			}
			
		}
		
		return rtrn;
	}

	private void write(Map<Gene, List<Double>> scoreByPosition, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(Gene g: scoreByPosition.keySet()) {
			writer.write(g.getName()+"\t"+g.toUCSCStrand());
			for(Double val: scoreByPosition.get(g)) {
				writer.write("\t"+val);
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>3) {
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
			File input=new File(args[1]);
			int index=Integer.parseInt(args[2]);
			String save=args[3];
			new ValuePerTSS(genes, input, index, save);
		}
		else {System.err.println(usage);}
		
	}
	static String usage=" args[0]=genes \n args[1]=count matrix \n args[2]=index position \n args[3]=save";
}

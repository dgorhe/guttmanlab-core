package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class RNADNACorrelationMatrix {

	//int minRNAFreq=10;
	int binResolution=1000000;
	
	static SingleInterval problemBin1=new SingleInterval("chr2:79490000-79500000");
	static SingleInterval problemBin2=new SingleInterval("chr11:3119270-3192250");
	static SingleInterval problemBin3=new SingleInterval("chr15:99734977-99736026");
	static SingleInterval problemBin4=new SingleInterval("chr3:5173978-5175025");
	static SingleInterval problemBin5=new SingleInterval("chr13:58176952-58178051");
	
	static Collection<SingleInterval> problemBins=new TreeSet<SingleInterval>();
	
	public RNADNACorrelationMatrix(BarcodingDataStreaming data, String save, Map<String, String> rnaGroups) throws IOException{
		MatrixWithHeaders mwh=data.getRNADNAContactMatrix(binResolution).getValue1();
		MatrixWithHeaders submatrix=subset(mwh, rnaGroups);
		submatrix=excludeProblemBins(submatrix, problemBins);
		
		submatrix.write(save+".contacts");
		
		MatrixWithHeaders correlation=correlationMatrix(submatrix);
		correlation.write(save+".pearson");
		
		MatrixWithHeaders spearman=spearmanMatrix(submatrix);
		spearman.write(save+".spearman");
		
		MatrixWithHeaders euclidian=euclidianDistance(submatrix);
		euclidian.write(save+".euclidian");
	
	}
	
	

	private MatrixWithHeaders excludeProblemBins(MatrixWithHeaders mwh, Collection<SingleInterval> problemBins) {
		Collection<String> list=new ArrayList<String>();
		for(SingleInterval problemBin: problemBins){
			String region=problemBin.bin(binResolution).toUCSC();
			list.add(region);
		}
		
		return mwh.excludeByColumnNames(list);
	}

	private MatrixWithHeaders subset(MatrixWithHeaders mwh, Map<String, String> rnaGroups) {
		Map<String, Collection<String>> groups=new TreeMap<String, Collection<String>>();
		for(String rna: rnaGroups.keySet()){
			String group=rnaGroups.get(rna);
			Collection<String> temp=new TreeSet<String>();
			if(groups.containsKey(group)){temp=groups.get(group);}
			temp.add(rna);
			groups.put(group, temp);
		}
		
		List<String> rows=new ArrayList<String>();
		
		for(String group: groups.keySet()){rows.add(group);}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, mwh.getColumnNames());
		
		for(String group: groups.keySet()){
			Collection<String> rnas=groups.get(group);
			double[] row=getVector(mwh, rnas);
			if(row!=null){
				rtrn.setRow(group, row);
			}
		}
		return rtrn;
	}

	private double[] getVector(MatrixWithHeaders mwh, Collection<String> rnas) {
		double[] rtrn=null;
		
		for(String rna: rnas){
			if(mwh.containsRow(rna)){
				double[] vals=mwh.getRow(rna);
				rtrn=merge(vals, rtrn);
			}
		}
		
		return rtrn;
	}

	private double[] merge(double[] vals1, double[] vals2) {
		if(vals2==null){return vals1;}
		if(vals1==null){return vals2;}
		
		double[] rtrn=new double[vals1.length];
		for(int i=0; i<vals1.length; i++){
			rtrn[i]=vals1[i]+vals2[i];
		}
		return rtrn;
	}

	private MatrixWithHeaders correlationMatrix(MatrixWithHeaders mwh) {
		PearsonsCorrelation p=new PearsonsCorrelation();
		//SpearmanCorrelation s=new SpearmanCorrelation(vals1, vals2);
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getRowNames());
		int counter=0;
		for(String gene1: mwh.getRowNames()){
			System.err.println(gene1);
			double[] vals1=mwh.getRow(gene1);
			for(String gene2: mwh.getRowNames()){
				double[] vals2=mwh.getRow(gene2);
				double correlation=p.correlation(vals1, vals2);
				rtrn.set(gene1, gene2, correlation);
			}
			counter++;
		}
		return rtrn;
	}
	
	private MatrixWithHeaders euclidianDistance(MatrixWithHeaders mwh) {
		EuclideanDistance e=new EuclideanDistance();
		//SpearmansCorrelation s=new SpearmansCorrelation();
		//SpearmanCorrelation s=new SpearmanCorrelation(vals1, vals2);
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getRowNames());
		int counter=0;
		for(String gene1: mwh.getRowNames()){
			System.err.println(gene1);
			double[] vals1=mwh.getRow(gene1);
			for(String gene2: mwh.getRowNames()){
				double[] vals2=mwh.getRow(gene2);
				double correlation=e.compute(vals1, vals2);
				rtrn.set(gene1, gene2, correlation);
			}
			counter++;
		}
		return rtrn;
	}
	
	private MatrixWithHeaders spearmanMatrix(MatrixWithHeaders mwh) {
		SpearmansCorrelation s=new SpearmansCorrelation();
		//SpearmanCorrelation s=new SpearmanCorrelation(vals1, vals2);
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getRowNames());
		int counter=0;
		for(String gene1: mwh.getRowNames()){
			System.err.println(gene1);
			double[] vals1=mwh.getRow(gene1);
			for(String gene2: mwh.getRowNames()){
				double[] vals2=mwh.getRow(gene2);
				double correlation=s.correlation(vals1, vals2);
				rtrn.set(gene1, gene2, correlation);
			}
			counter++;
		}
		return rtrn;
	}
	
	
	private static Map<String, String> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(String line: lines){
			rtrn.put(line.split("\t")[0].replaceAll("\"", ""), line.split("\t")[1]);
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		problemBins.add(problemBin1);
		problemBins.add(problemBin2);
		problemBins.add(problemBin3);
		problemBins.add(problemBin4);
		problemBins.add(problemBin5);
		if(args.length>2){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String save=args[1];
			Map<String, String> genes=parse(args[2]);
			new RNADNACorrelationMatrix(data, save, genes);
		}
		else{System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=barcoding data \n args[1]=save \n args[2]=genes";
}

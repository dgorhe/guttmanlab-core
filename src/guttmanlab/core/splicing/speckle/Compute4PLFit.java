package guttmanlab.core.splicing.speckle;

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

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.EnumeratedDistribution;

import flanagan.physchem.ImmunoAssay;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;

public class Compute4PLFit {

	private static void writeTable(File[] order, Map<String, SingleInterval> coordinates, Map<String, Pair<Integer>>[] scores, String save, int minCutoff, Map<String, Double> ratios, Map<SingleInterval, Double> speckleScores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		//Map<String, Map<SingleInterval, Double>> c50ByGene=new TreeMap<String, Map<SingleInterval, Double>>();
		
		System.err.println(coordinates.size());
		//Collection<String> allGenes=new TreeSet<String>();
		
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		
		writer.write("Gene Name\t coordinates \t length \t min total\tmin ratio\texpression\tspeckle score\tobserved range\tfit range\tderivative\tslope\tslope [0-1]\ttop\tbottom\tC50\thill\tr2");
		for(int i=0; i<scores.length; i++) {
			writer.write("\t"+x[i]);
		}
		for(int i=0; i<scores.length; i++) {
			writer.write("\t"+x[i]);
		}
		writer.write("\n");
		
		int counter=0;
		for(String gene: coordinates.keySet()) {
			SingleInterval region=coordinates.get(gene);
			double minTotal=getMinTotal(scores, gene);
			double[] y=new double[order.length];
			if(minTotal>=minCutoff) {
				for(int i=0; i<scores.length; i++) {
					double ratio=getRatio(scores[i], gene);
					y[i]=ratio;
				}
				double[] fitParms=fit4PL(x,y);
				if(fitParms!=null) {
					//get name
					//String parentName=getParentName(gene);
					//add(c50ByGene, parentName, region, fitParms[2]);
					double speckleScore=-1;
					if(speckleScores.containsKey(region)) {speckleScore=speckleScores.get(region);}
					double ratio=ratios.get(gene);
					writer.write(gene+"\t"+region.toUCSC()+"\t"+region.getGenomicLength()+"\t"+minTotal+"\t"+ratio);
					double[] ynorm=norm(y, fitParms);
					double expression=minTotal/(double)region.getGenomicLength();
					double slope=slope(fitParms);
					double slope2=slope(fitParms, 0, 1);
					double range1=range(y);
					double range2=range(fitParms[0], fitParms[1]);
					double derivative=derivative(x,y);
					writer.write("\t"+expression+"\t"+speckleScore+"\t"+range1+"\t"+range2+"\t"+derivative+"\t"+slope+"\t"+slope2+"\t"+fitParms[0]+"\t"+fitParms[1]+"\t"+fitParms[2]+"\t"+fitParms[3]+"\t"+fitParms[4]);
					for(int i=0; i<y.length; i++) {writer.write("\t"+y[i]);}
					for(int i=0; i<y.length; i++) {writer.write("\t"+ynorm[i]);}
					//System.out.println(gene+"\t"+derivative);
					writer.write("\n");
				}
			}
			counter++;
			if(counter%1000 ==0) {System.err.println(counter);}
			//if(counter>100) {System.err.println("done"); break;}
		}
		
		
		writer.close();
		
		//write(save+".c50ByGene.scores", c50ByGene);
	}
	
	
	private static void writeTable(File[] order, Map<String, SingleInterval> coordinates, Map<String, Pair<Integer>>[] scores, String save, Map<String, Collection<String>> genesToJunctions, Map<String, Pair<Integer>>[] junctionScores, Map<String, Integer> sum) throws IOException {
		Map<String, String> speckleScores=new TreeMap<String, String>();
		writeTable(order, coordinates, scores, save, genesToJunctions, junctionScores, sum, speckleScores);
	}
	
	private static void writeTable(File[] order, Map<String, SingleInterval> coordinates, Map<String, Pair<Integer>>[] scores, String save, Map<String, Collection<String>> genesToJunctions, Map<String, Pair<Integer>>[] junctionScores, Map<String, Integer> sum, Map<String, String> linesByGene) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<String, Integer> junctionSum=sum(junctionScores);
		
		FileWriter writerBG=new FileWriter(save+".bedgraph");
		FileWriter writerJunc=new FileWriter(save+".junctions");
		
		double r2cutoff=0.9;
		
		System.err.println(coordinates.size());
		//Collection<String> allGenes=new TreeSet<String>();
		
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		
		//double[] x= {0, 11.75, 16.92, 21.50, 26.92, 31.92, 47.63, 62.57, 77.50, 92.50, 122.50, 250.92};
		
		/*writer.write("Gene Name\t coordinates \t length \t min total\tmin ratio\texpression\tspeckle score\tobserved range\tfit range\tderivative\tslope\tslope [0-1]\ttop\tbottom\tC50\thill\tr2");
		
		for(int i=0; i<scores.length; i++) {
			writer.write("\t"+x[i]);
		}
		for(int i=0; i<scores.length; i++) {
			writer.write("\t"+x[i]);
		}
		writer.write("\n");*/
		
		int counter=0;
		for(String gene: coordinates.keySet()) {
			SingleInterval coordinate=coordinates.get(gene);
		
			String line="";
			if(linesByGene.containsKey(gene)) {line=linesByGene.get(gene);}
			
			double[] y=getY(scores, gene);
			
			//System.err.println(x.length+" "+y.length);
			
			double[] fitParms=fit4PL(x,y);
			if(fitParms!=null && fitParms[4]>r2cutoff) {
				double slope2=slope(fitParms, 0, 1);
				
				writer.write(line+"\t"+gene+"\t"+coordinate.toUCSC()+"\t"+sum.get(gene)+"\t"+coordinates.get(gene).getGenomicLength()+"\t"+fitParms[2]+"\t"+slope2+"\t"+fitParms[3]);
				
				List<Double> jVals=new ArrayList<Double>();
				//jVals.add(fitParms[2]);
				for(String junction: genesToJunctions.get(gene)) {
					double[] y2=getY(junctionScores, junction);
					int sumJ=junctionSum.get(junction);
					double[] fitParms2=fit4PL(x,y2);
					if(fitParms2!=null && fitParms2[4]>r2cutoff) {
						SingleInterval jr=new SingleInterval(junction);
						SingleInterval trimmed=jr.getMidPoint();
						SingleInterval binned=trimmed.bin(1000);
						writerBG.write(binned.toBedgraph(50-fitParms2[2])+"\n");
						//double slope=slope(fitParms2, 0, 1);
						jVals.add(fitParms2[2]);
						
						writerJunc.write(line+"\t"+junction+"\t"+sumJ+"\t"+fitParms2[2]+"\t"+fitParms2[4]);
						for(int i=0; i<y2.length; i++) {writerJunc.write("\t"+y2[i]);}
						writerJunc.write("\n");
						
					}
				}
				
				if(jVals.isEmpty()) {jVals.add(fitParms[2]);}
				
				double median=Statistics.quantile(jVals, 0.5);
				writer.write("\t"+median);
				//for(Double val: jVals) {writer.write("\t"+val);}
				
				for(int i=0; i<y.length; i++) {writer.write("\t"+y[i]);}
				
				//System.out.println(gene+"\t"+derivative);
				writer.write("\n");
			}
			
			counter++;
			if(counter%1000 ==0) {System.err.println(counter);}
			//if(counter>100) {System.err.println("done"); break;}
		}
		
		
		writer.close();
		writerBG.close();	
		writerJunc.close();
		//write(save+".c50ByGene.scores", c50ByGene);
	}
	
	
	private static double minCount(Map<String, Pair<Integer>>[] scores, String gene) {
		double[] counts=new double[scores.length];
		for(int i=0; i<scores.length; i++) {
			double count=getCount(scores[i], gene);
			counts[i]=count;
		}
		return Statistics.min(counts);
	}
	
	private static double[] getY(Map<String, Pair<Integer>>[] scores, String gene) {
		double[] y=new double[scores.length];
		for(int i=0; i<scores.length; i++) {
			double ratio=getRatio(scores[i], gene);
			y[i]=ratio;
		}
		return y;
	}
	
	private static double[] getYFraction(Map<String, Pair<Double>>[] scores, String gene) {
		double[] y=new double[scores.length];
		for(int i=0; i<scores.length; i++) {
			double ratio=getRatioDouble(scores[i], gene);
			y[i]=ratio;
		}
		return y;
	}
	
	private static int[] getSum(Map<String, Pair<Integer>>[] scores, String gene) {
		int[] y=new int[scores.length];
		for(int i=0; i<scores.length; i++) {
			int ratio=getSum(scores[i], gene);
			y[i]=ratio;
		}
		return y;
	}


	private static void write(String save, Map<String, Map<SingleInterval, Double>> c50ByGene) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String gene: c50ByGene.keySet()) {
			writer.write(gene);
			for(SingleInterval r: c50ByGene.get(gene).keySet()) {
				writer.write("\t"+r.toUCSC()+"\t"+c50ByGene.get(gene).get(r));
			}
			writer.write("\n");
		}
		
		writer.close();
	}


	private static String getParentName(String gene) {
		return gene.split("_j")[0];
	}


	private static void add(Map<String, Map<SingleInterval, Double>> c50ByGene, String parentName, SingleInterval junction, double d) {
		if(!c50ByGene.containsKey(parentName)) {c50ByGene.put(parentName, new TreeMap<SingleInterval, Double>());}
		Map<SingleInterval, Double> map=c50ByGene.get(parentName);
		map.put(junction, d);
	}


	private static double range(double[] vals) {
		return Statistics.max(vals)-Statistics.min(vals);
	}
	
	private static double range(double top, double bottom) {
		return top-bottom;
	}

	private static double[] norm(double[] y, double[] fitParms) {
		///top, bottom, C50, HillSlope
		
		double[] rtrn=new double[y.length];
		
		for(int i=0; i<y.length; i++) {
			double val=y[i];
			double scaleFactor=1.0/(fitParms[0]-fitParms[1]);
			rtrn[i]=(val-fitParms[1])*scaleFactor;
		}
		
		return rtrn;
	}

	
	private static double slope(double[] fitParms, double min, double max) {
		double a = max;
		double b = fitParms[3];
		double c = fitParms[2];
		double d = min;
		
		///top, bottom, C50, HillSlope
		
		double num=((b*(a-d)/c))* Math.pow((b-1)/(b+1), (1-(1/b)));
		double denom=Math.pow(1+((b-1)/(b+1)), 2);
		double slope=Math.abs(num/denom);
		return slope;
	}
	
	private static double slope(double[] fitParms) {
		double a = fitParms[0];
		double b = fitParms[3];
		double c = fitParms[2];
		double d = fitParms[1];
		
		///top, bottom, C50, HillSlope
		
		double num=((b*(a-d)/c))* Math.pow((b-1)/(b+1), (1-(1/b)));
		
		
		//double num=Math.pow(((a-d)/c)* (b-1)/(b+1), (1-(1/b)));
		
		double denom=Math.pow(1+((b-1)/(b+1)), 2);
		
		
				
		//double num=-1*(b*Math.pow((b - 1)/(b + 1), (b - 1)*(a-d)));
		//double denom=c*Math.pow(Math.pow(Math.pow((b-1)/(b+1), 1/b), b)+1,2);
				
				//(c*((((b - 1)/(b + 1))^(1/b))^b + 1)^2);
		
		double slope=Math.abs(num/denom);
		//double slope = -(b*(((b - 1)/(b + 1))^(1/b))^(b - 1)*(a - d))/(c*((((b - 1)/(b + 1))^(1/b))^b + 1)^2);
		return slope;
	}

	private static double[] fit4PL(double[] x, double[] y) {
		try {
		ImmunoAssay assay = new ImmunoAssay(x, y);
		
		assay.suppressPrint();
		assay.suppressYYplot();
		
		try {
			assay.fourParameterLogisticFit();
		}catch(java.awt.HeadlessException ex) {}
		double[] vals=assay.getModelParameterValues(); ///top, bottom, C50, HillSlope
		
		//double ss=assay.getTotalSumOfWeightedSquares();
		double r2=assay.getCoefficientOfDetermination();
		
		double[] rtrn=new double[vals.length+1];
		for(int i=0; i<vals.length; i++) {
			rtrn[i]=vals[i];
		}
		rtrn[vals.length]=r2;
		return rtrn;
		}catch(java.lang.IllegalArgumentException ex) {
			double[] rtrn= {-1,-1,-1,-1,-1};
			return rtrn;
		}
	}

	private static double derivative(double[] x, double[] y) {
		Collection<Double> vals=new ArrayList<Double>();
		for(int i=1; i<x.length; i++) {
			double num=y[i]-y[i-1];
			double denom=x[i]-x[i-1];
			double rate=num/denom;
			vals.add(rate);
		}
		return Statistics.max(vals);
	}
	
	private static double getMinTotal(Map<String, Pair<Integer>>[] scores, String gene) {
		double[] totals=new double[scores.length];
		for(int i=0; i<scores.length; i++) {
			double total=getTotal(scores[i], gene);
			totals[i]=total;
		}
		return Statistics.min(totals);
	}

	private static double getTotal(Map<String, Pair<Integer>> map, String gene) {
		Pair<Integer> pair=new Pair<Integer>(0,0);
		if(map.containsKey(gene)) {
			pair=map.get(gene);
		}
		return pair.getValue1()+pair.getValue2();
	}

	private static int getSum(Map<String, Pair<Integer>> map, String gene) {
		int sum=0;
		if(map.containsKey(gene)) {
			Pair<Integer> pair=map.get(gene);
			sum=pair.getValue1()+pair.getValue2();
			//ratio=(double)pair.getValue1()/sum;
		}
		return sum;
	}
	
	private static double getRatioDouble(Map<String, Pair<Double>> map, String gene) {
		double ratio=-1;
		if(map.containsKey(gene)) {
			Pair<Double> pair=map.get(gene);
			double sum=pair.getValue1()+pair.getValue2();
			if(sum>0) {
				ratio=(double)pair.getValue1()/sum;
			}
		}
		return ratio;
	}
	
	private static double getRatio(Map<String, Pair<Integer>> map, String gene) {
		double ratio=-1;
		if(map.containsKey(gene)) {
			Pair<Integer> pair=map.get(gene);
			double sum=pair.getValue1()+pair.getValue2();
			ratio=(double)pair.getValue1()/sum;
		}
		return ratio;
	}
	
	private static Pair<Integer> getPair(Map<String, Pair<Integer>> map, String gene) {
		//double ratio=0;
		if(map.containsKey(gene)) {
			Pair<Integer> pair=map.get(gene);
			return pair;
		}
		return new Pair<Integer>(0,0);
	}
	
	private static double getCount(Map<String, Pair<Integer>> map, String gene) {
		double ratio=0;
		if(map.containsKey(gene)) {
			Pair<Integer> pair=map.get(gene);
			double sum=pair.getValue1()+pair.getValue2();
			ratio=sum;
		}
		return ratio;
	}
	
	
	private static Map<String, double[]>[] parsePermScores(File[] files) throws NumberFormatException, IOException {
		Map<String, double[]>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn[i]=parsePermScore(files[i]);
		}
		
		return rtrn;
	}
	
	
	private static Map<String, Double>[] parsePerm(File[] files, int index) throws NumberFormatException, IOException {
		Map<String, Double>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			rtrn[i]=new TreeMap<String, Double>();
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				String name=tokens[0]+"_"+tokens[1];
				double score= Double.parseDouble(tokens[index+10]);
				rtrn[i].put(name, score);
			}
			reader.close();	
		}
		
		
		
		return rtrn;
	}
	
	private static Map<String, int[]> parsePermScore(File file , int numPerm) throws NumberFormatException, IOException {
		Map<String, int[]> rtrn=new TreeMap<String, int[]>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			//if(junctionsToUse.contains(name)) {
				//int[] vals= new int[tokens.length-7];
				
				int[] vals=new int[Math.min(tokens.length-7, numPerm)];
				
				for(int i=0; i<vals.length; i++) {
					vals[i]=Integer.parseInt(tokens[i+7]);
				}
				
				
				//Pair<Integer> pair=new Pair<Integer>(Integer.parseInt(tokens[5]), Integer.parseInt(tokens[6]));
				rtrn.put(name, vals);
			//}
		}
		reader.close();
		
		return rtrn;
	}
	
	private static Map<String, double[]> parsePermScore(File file) throws NumberFormatException, IOException {
		Map<String, double[]> rtrn=new TreeMap<String, double[]>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			double[] vals= new double[tokens.length-11];
			
			for(int i=11; i<tokens.length; i++) {
				vals[i-11]=Double.parseDouble(tokens[i]);
			}
			
			
			//Pair<Integer> pair=new Pair<Integer>(Integer.parseInt(tokens[5]), Integer.parseInt(tokens[6]));
			rtrn.put(name, vals);
		}
		reader.close();
		
		return rtrn;
	}
	
	private static Map<String, double[]>[] parseEstimates(File[] files) throws NumberFormatException, IOException {
		Map<String, double[]>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn[i]=parseEstimate(files[i]);
		}
		
		return rtrn;
	}
	
	private static Map<String, double[]> parseEstimate(File file) throws NumberFormatException, IOException {
		Map<String, double[]> rtrn=new TreeMap<String, double[]>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			double[] vals= new double[3];
			vals[0]=Double.parseDouble(tokens[8]);
			vals[1]=Double.parseDouble(tokens[9]);
			vals[2]=Double.parseDouble(tokens[10]);
			//Pair<Integer> pair=new Pair<Integer>(Integer.parseInt(tokens[5]), Integer.parseInt(tokens[6]));
			rtrn.put(name, vals);
		}
		reader.close();
		
		return rtrn;
	}
	
	
	private static Map<String, Pair<Integer>>[] parse(File[] files) throws NumberFormatException, IOException {
		Map<String, Pair<Integer>>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn[i]=parse(files[i]);
		}
		
		return rtrn;
	}
	
	private static Map<String, Pair<Double>>[] parseFractions(File[] files) throws NumberFormatException, IOException {
		Map<String, Pair<Double>>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn[i]=parseFractions(files[i]);
		}
		
		return rtrn;
	}
	
	private static Map<String, Pair<Double>>[] parseJunctionFractions(File[] files) throws NumberFormatException, IOException {
		Map<String, Pair<Double>>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn[i]=parseJunctionFractions(files[i]);
		}
		
		return rtrn;
	}
	
	private static Map<String, Pair<Integer>>[] parseTotalCounts(File[] files) throws NumberFormatException, IOException {
		Map<String, Pair<Integer>>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn[i]=parseTotalCounts(files[i]);
		}
		
		return rtrn;
	}
	
	private static Map<String, Pair<Integer>>[] parseJunctionIntron(File[] files) throws NumberFormatException, IOException {
		Map<String, Pair<Integer>>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn[i]=parseJunctionCounts(files[i]);
		}
		
		return rtrn;
	}
	
	private static Map<String, Integer> parseJunctionCounts(File[] files) throws NumberFormatException, IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(int i=0; i<files.length; i++) {
			
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				String name=tokens[0]+"_"+tokens[1];
				int count=Integer.parseInt(tokens[3]);
				if(rtrn.containsKey(name)) {
					count+=rtrn.get(name);
				}
				rtrn.put(name, count);
			}
			reader.close();
			
			
		}
		
		return rtrn;
	}
	
	
	private static Map<String, Integer>[] parseJunctionCountArrays(File[] files) throws NumberFormatException, IOException {
		Map<String, Integer>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn[i]=new TreeMap<String, Integer>();
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				String name=tokens[0]+"_"+tokens[1];
				int count=Integer.parseInt(tokens[3]);
				rtrn[i].put(name, count);
			}
			reader.close();
		}
		
		return rtrn;
	}
	
	
	private static Map<String, Double> parseRatios(File[] files) throws IOException {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			
			
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				String name=tokens[0]+"_"+tokens[1];
				//SingleInterval region=new SingleInterval(tokens[1]);
				double ratio=Double.parseDouble(tokens[4]);
				rtrn.put(name,  ratio);
			}
			reader.close();
			
		}
		
		return rtrn;
	}
	
	
	private static Map<String, Boolean> parseUniqueness(File[] files) throws IOException {
		Map<String, Boolean> rtrn=new TreeMap<String, Boolean>();
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			
			
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				String name=tokens[0]+"_"+tokens[1];
				//SingleInterval region=new SingleInterval(tokens[1]);
				boolean b=Boolean.parseBoolean(tokens[2]);
				rtrn.put(name,  b);
			}
			reader.close();
			
		}
		
		return rtrn;
	}
	
	private static Map<String, Double> parseCounts(File[] files) throws IOException {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			
			
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				String name=tokens[0]+"_"+tokens[1];
				//SingleInterval region=new SingleInterval(tokens[1]);
				//boolean b=Boolean.parseBoolean(tokens[2]);
				double count=Double.parseDouble(tokens[3]);
				rtrn.put(name, count);
			}
			reader.close();
			
		}
		
		return rtrn;
	}
	
	private static Map<String, SingleInterval> parseCoordinates(File[] files) throws IOException {
		Map<String, SingleInterval> rtrn=new TreeMap<String, SingleInterval>();
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			
			
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				String name=tokens[0]+"_"+tokens[1];
				SingleInterval region=new SingleInterval(tokens[1]);
				rtrn.put(name,  region);
			}
			reader.close();
			
		}
		
		return rtrn;
	}

	
	private static Map<String, Pair<Integer>> parse(File file) throws NumberFormatException, IOException {
		Map<String, Pair<Integer>> rtrn=new TreeMap<String, Pair<Integer>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			Pair<Integer> pair=new Pair<Integer>(Integer.parseInt(tokens[2]), Integer.parseInt(tokens[3]));
			rtrn.put(name, pair);
		}
		reader.close();
		
		return rtrn;
	}
	
	
	private static Map<String, Pair<Integer>> parseTotalCounts(File file) throws NumberFormatException, IOException {
		Map<String, Pair<Integer>> rtrn=new TreeMap<String, Pair<Integer>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			Pair<Integer> pair=new Pair<Integer>(Integer.parseInt(tokens[3]), Integer.parseInt(tokens[4]));
			rtrn.put(name, pair);
		}
		reader.close();
		
		return rtrn;
	}
	
	private static Map<String, Pair<Integer>> parseJunctionCounts(File file) throws NumberFormatException, IOException {
		Map<String, Pair<Integer>> rtrn=new TreeMap<String, Pair<Integer>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			Pair<Integer> pair=new Pair<Integer>(Integer.parseInt(tokens[2]), Integer.parseInt(tokens[4]));
			rtrn.put(name, pair);
		}
		reader.close();
		
		return rtrn;
	}
	
	private static Map<String, Pair<Double>> parseFractions(File file) throws NumberFormatException, IOException {
		Map<String, Pair<Double>> rtrn=new TreeMap<String, Pair<Double>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			Pair<Double> pair=new Pair<Double>(Double.parseDouble(tokens[7]), Double.parseDouble(tokens[8]));
			rtrn.put(name, pair);
		}
		reader.close();
		
		return rtrn;
	}
	
	
	private static Map<String, Pair<Double>> parseJunctionFractions(File file) throws NumberFormatException, IOException {
		Map<String, Pair<Double>> rtrn=new TreeMap<String, Pair<Double>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0]+"_"+tokens[1];
			Pair<Double> pair=new Pair<Double>(Double.parseDouble(tokens[3]), (Double.parseDouble(tokens[9])*150));
			rtrn.put(name, pair);
		}
		reader.close();
		
		return rtrn;
	}
	
	private static File[] order(List<String> files) {
		File[] rtrn=new File[files.size()];
		
		int counter=0;
		for(String f: files) {
			rtrn[counter]=new File(f);
			counter++;
		}
		
		return order(rtrn);
	}
	
	private static File[] order(File[] files) {
		Map<Integer, File> ordered=new TreeMap<Integer, File>();
		
		for(int i=0; i<files.length; i++) {
			String name=files[i].getName().split("min")[0];
			Integer num=Integer.parseInt(name);
			ordered.put(num, files[i]);
		}
		
		
		File[] rtrn=new File[files.length];
		
		int index=0;
		for(Integer key: ordered.keySet()) {
			rtrn[index]=ordered.get(key);
			index++;
		}
		return rtrn;
	}
	
	/*public static void run(File[] files, String save, int minCutoff, BarcodingDataStreaming data, Kmer kmer, int binResolution) throws IOException {
		files=order(files);
		Map<String, Pair<Integer>>[] scores=parse(files);
		Map<String, SingleInterval> coordinates=parseCoordinates(files);
		Map<String, Double> ratios=parseRatios(files);
		Map<String, Integer> sum=sum(scores);
		
		Map<SingleInterval, Double> speckleScores=DistanceToNuclearBody.distance(data, kmer, binResolution, coordinates.values());
		
		
		
		
		Map<String, Collection<String>> genesToJunctions=parseJunctions(files);
		Map<String, Pair<Integer>>[] junctionScores=parseJunctionScores(files);
		
		writeTable(files, coordinates, scores, save, genesToJunctions, junctionScores, sum, speckleScores);
		
		
		//writeTable(files, coordinates, scores, save, minCutoff, ratios, speckleScores);
		
		
		
		
		System.exit(0);
	}*/
	
	
	/*public static void run(File[] files, String save, int minCutoff) throws IOException {
		files=order(files);
		Map<String, Pair<Integer>>[] scores=parse(files);
		Map<String, SingleInterval> coordinates=parseCoordinates(files);
		Map<String, Double> ratios=parseRatios(files);
		Map<String, Integer> sum=sum(scores);
		
		Map<String, Collection<String>> genesToJunctions=parseJunctions(files);
		Map<String, Pair<Integer>>[] junctionScores=parseJunctionScores(files);
		
		writeTable(files, coordinates, scores, save, genesToJunctions, junctionScores, sum);
		
		System.exit(0);
	}*/
	
	public static void run(File[] files, String save, int minCutoff, Map<String, String> linesByGene) throws IOException {
		FileWriter writer=new FileWriter(save);
		//FileWriter bg=new FileWriter(save+".actual.bedgraph");
		//FileWriter bgLow=new FileWriter(save+".low.bedgraph");
		//FileWriter bgHigh=new FileWriter(save+".high.bedgraph");
		
		files=order(files);
		Map<String, Pair<Integer>>[] scores=parse(files); //gene+coordinate
		
		Map<String, Double> junctionCounts=getJunctionCount(scores);
		
		Map<String, Double> maxVals=getMax(junctionCounts);
		
		int numPerm=10;
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		Map<String, int[]>[] perms=parsePerms(files, numPerm);
		
		int counter=0;
		Collection<String> junctionsToUse=new TreeSet<String>();
		Map<String, Collection<String>> junctionsByGene=new TreeMap<String, Collection<String>>();
		for(String gene: junctionCounts.keySet()) {
			String name=gene.split("_")[0];
			SingleInterval region=new SingleInterval(gene.split("_")[1]);
			double maxVal=maxVals.get(name);
			double[] y=getY(scores, gene);
			//int[] sum=getSum(scores, gene);
			//double max=Statistics.max(y);
			double min=Statistics.min(y);
			double count=junctionCounts.get(gene);
			double ratio=count/maxVal;
			if(count>10 && ratio>0.1 && min>=0) {
				/*System.out.print(gene+"\t"+count);
				for(int i=0; i<y.length; i++) {System.out.print("\t"+y[i]);}
				System.out.println();*/
				add(junctionsByGene, name, gene);
				junctionsToUse.add(gene);
				
				double[] fit=fit4PL(x,y);
				//bg.write(region.toBedgraph(50-fit[2])+"\n");
				
				
				List<Double> c50=getPerms(perms, scores, x, gene);
				
				double minPerm=-1;
				double maxPerm=-1;
				if(!c50.isEmpty()) {
					//writer.write("\t"+Statistics.quantile(c50, 0.5)+"\t"+Statistics.quantile(c50, 0.05)+"\t"+Statistics.quantile(c50, 0.95));
					minPerm=Statistics.min(c50);
					maxPerm=Statistics.max(c50);
				}
				writer.write(gene+"\t"+fit[2]+"\t"+fit[4]+"\t"+minPerm+"\t"+maxPerm);
				
				for(int i=0; i<y.length; i++) {
					writer.write("\t"+y[i]);
				}
				writer.write("\n");
				
			}
			counter++;
			if(counter%1000==0) {System.err.println(counter+" "+junctionCounts.size());}
		}
		
		
		
		/*int counter=0;
		for(String gene: junctionsByGene.keySet()) {
			Collection<String> junctions=junctionsByGene.get(gene);
			double sumCount=sum(junctions, junctionCounts);
			SingleInterval region=coordinates(junctions);
			
			if(sumCount>100) {
				double[] yOriginal=getYScaledFractions(scores, junctions, junctionCounts);
				double[] fit=fit4PL(x,yOriginal);
				writer.write(gene+"\t"+sumCount+"\t"+fit[2]+"\t"+fit[4]);
				bg.write(region.toBedgraph(fit[2])+"\n");
				
				
				List<Double> c50=new ArrayList<Double>();
				for(int i=0; i<numPerm; i++) {
					double[] y=getYScaledFractions(scores, junctions, junctionCounts, perms, i);
					fit=fit4PL(x,y);
					if(fit[2]!=-1) {
						c50.add(fit[2]);
					}
					//writer.write("\t"+fit[4]+"\t"+fit[2]);
				}
				
				writer.write("\t"+Statistics.quantile(c50, 0.5)+"\t"+Statistics.quantile(c50, 0.05)+"\t"+Statistics.quantile(c50, 0.95));
				bgLow.write(region.toBedgraph(Statistics.quantile(c50, 0.05))+"\n");
				bgHigh.write(region.toBedgraph(Statistics.quantile(c50, 0.95))+"\n");
				
				for(int i=0; i<yOriginal.length; i++) {
					writer.write("\t"+yOriginal[i]);
				}
				
				writer.write("\n");
				
				
			}
			counter++;
			if(counter%10==0) {System.err.println(gene+" "+counter+" "+junctionsByGene.size());}
		}*/
		
		writer.close();
		//bg.close();
		//bgLow.close();
		//bgHigh.close();
		System.exit(0);
	}
	
	
	public static void writeFraction(List<String> fileNames, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		FileWriter bg=new FileWriter(save+".bedgraph");
		
		File[] files=order(fileNames);
		Map<String, Pair<Double>>[] scores=parseFractions(files); //gene+coordinate
		Map<String, Integer> junctionCounts=parseJunctionCounts(files);
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		int counter=0;
		for(String gene: junctionCounts.keySet()) {
			int junctionCount=junctionCounts.get(gene);
			if(junctionCount>0) {
				String name=gene.split("_")[0];
				SingleInterval region=new SingleInterval(gene.split("_")[1]);
				double[] y=getYFraction(scores, gene);
				double min=Statistics.min(y);
				double[] fit=fit4PL(x,y);	
				writer.write(name+"\t"+region.toUCSC()+"\t"+junctionCount+"\t"+fit[2]+"\t"+fit[4]);
				for(int i=0; i<y.length; i++) {
					writer.write("\t"+y[i]);
				}
				writer.write("\n");
				bg.write(region.toBedgraph(fit[2])+"\n");
			}
			counter++;
			
			if(counter%1000==0) {System.err.println(counter+" "+junctionCounts.size());}
		}
		writer.close();
		bg.close();
	}
	
	public static void processPerm(List<String> fileNames, String save, int index, Map<String, Double> rpkmMap) throws NumberFormatException, IOException {
		FileWriter writer=new FileWriter(save);
		//FileWriter bg=new FileWriter(save+".bedgraph");
		
		File[] files=order(fileNames);
		//Map<String, Pair<Double>>[] scores=parseJunctionFractions(files); //gene+coordinate
		
		Map<String, Double>[] scores=parsePerm(files, index);
		Map<String, Integer> junctionCounts=parseJunctionCounts(files);
		Map<String, Integer> geneMax=getGeneMax(junctionCounts);
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		int counter=0;
		for(String gene: junctionCounts.keySet()) {
			int junctionCount=junctionCounts.get(gene);
			if(junctionCount>0) {
				String name=gene.split("_")[0];
				if(rpkmMap.containsKey(name)) {
					double rpkm=rpkmMap.get(name);
					int geneMaxScore=geneMax.get(name);
					double junctionFraction=(double)junctionCount/(double)geneMaxScore;
					SingleInterval region=new SingleInterval(gene.split("_")[1]);
					double[] y=getYPerm(scores, gene);
					double min=Statistics.min(y);
					double[] fit=fit4PL(x,y);	
					/*if(min>=0 && fit[2]>=0 && junctionFraction>0.1) {
						bg.write(region.toBedgraph(50-fit[2])+"\n");
					}*/
					if(min>=0) {
						writer.write(name+"\t"+region.toUCSC()+"\t"+rpkm+"\t"+junctionCount+"\t"+junctionFraction+"\t"+fit[2]+"\t"+fit[4]);
						for(int i=0; i<y.length; i++) {
							writer.write("\t"+y[i]);
						}
						writer.write("\n");
					}
				}
			}
			counter++;
			if(counter%1000==0) {System.err.println(counter+" "+junctionCounts.size());}
		}
		writer.close();
		//bg.close();
	}

	

	
	
	private static double[] getYPerm(Map<String, Double>[] scores, String gene) {
		double[] rtrn=new double[scores.length];
		
		for(int i=0; i<scores.length; i++) {
			if(scores[i].containsKey(gene)) {
				rtrn[i]=scores[i].get(gene);
			}
			else {rtrn[i]=-1;}
		}
		
		return rtrn;
	}


	public static void writeJunctionFraction(List<String> fileNames, String save, Map<String, Double> rpkmMap) throws IOException {
		FileWriter writer=new FileWriter(save);
		//FileWriter bg=new FileWriter(save+".bedgraph");
		
		File[] files=order(fileNames);
		Map<String, Pair<Double>>[] scores=parseJunctionFractions(files); //gene+coordinate
		Map<String, Integer> junctionCounts=parseJunctionCounts(files);
		Map<String, Integer> geneMax=getGeneMax(junctionCounts);
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		int counter=0;
		for(String gene: junctionCounts.keySet()) {
			int junctionCount=junctionCounts.get(gene);
			if(junctionCount>0) {
				String name=gene.split("_")[0];
				if(rpkmMap.containsKey(name)) {
					double rpkm=rpkmMap.get(name);
					int geneMaxScore=geneMax.get(name);
					double junctionFraction=(double)junctionCount/(double)geneMaxScore;
					SingleInterval region=new SingleInterval(gene.split("_")[1]);
					double[] y=getYFraction(scores, gene);
					double min=Statistics.min(y);
					if(min>=0) {
						double[] fit=fit4PL(x,y);
						if(fit[2]>x[x.length-1]) {fit[2]=-1; fit[3]=-1; fit[4]=-1;}
						writer.write(name+"\t"+region.toUCSC()+"\t"+rpkm+"\t"+junctionCount+"\t"+junctionFraction+"\t"+fit[2]+"\t"+fit[4]+"\t"+fit[3]);
						for(int i=0; i<y.length; i++) {
							writer.write("\t"+y[i]);
						}
						writer.write("\n");
					}
					/*if(min>=0 && fit[2]>=0 && junctionFraction>0.1) {
						bg.write(region.toBedgraph(50-fit[2])+"\n");
					}*/
				}
			}
			counter++;
			if(counter%1000==0) {System.err.println(counter+" "+junctionCounts.size());}
		}
		writer.close();
		//bg.close();
	}
	
	
	public static void writeJunctionCounts(List<String> fileNames, String save, Map<String, Double> rpkmMap) throws IOException {
		FileWriter writer=new FileWriter(save);
		//FileWriter bg=new FileWriter(save+".bedgraph");
		
		File[] files=order(fileNames);
		Map<String, Integer>[] scores=parseJunctionCountArrays(files); //gene+coordinate
		Map<String, Integer> junctionCounts=parseJunctionCounts(files);
		Map<String, Integer> geneMax=getGeneMax(junctionCounts);
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		int counter=0;
		for(String gene: junctionCounts.keySet()) {
			int junctionCount=junctionCounts.get(gene);
			if(junctionCount>0) {
				String name=gene.split("_")[0];
				if(rpkmMap.containsKey(name)) {
					double rpkm=rpkmMap.get(name);
					int geneMaxScore=geneMax.get(name);
					double junctionFraction=(double)junctionCount/(double)geneMaxScore;
					SingleInterval region=new SingleInterval(gene.split("_")[1]);
					double[] y=getYCounts(scores, gene);
					double min=Statistics.min(y);
					if(min>=0) {
						double[] fit=fit4PL(x,y);
						if(fit[2]>x[x.length-1]) {fit[2]=-1; fit[3]=-1; fit[4]=-1;}
						writer.write(name+"\t"+region.toUCSC()+"\t"+rpkm+"\t"+junctionCount+"\t"+junctionFraction+"\t"+fit[2]+"\t"+fit[4]+"\t"+fit[3]);
						for(int i=0; i<y.length; i++) {
							writer.write("\t"+y[i]);
						}
						writer.write("\n");
					}
				}
			}
			counter++;
			if(counter%1000==0) {System.err.println(counter+" "+junctionCounts.size());}
		}
		writer.close();
		//bg.close();
	}
	
	private static double[] getYCounts(Map<String, Integer>[] scores, String gene) {
		double[] rtrn=new double[scores.length];
		for(int i=0; i<rtrn.length; i++) {
			double score=0;
			if(scores[i].containsKey(gene)){
				score=scores[i].get(gene);
			}
			rtrn[i]=score;
		}
		
		return rtrn;
	}


	private static Map<String, Integer> getGeneMax(Map<String, Integer> junctionCounts) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(String junction: junctionCounts.keySet()) {
			int score=junctionCounts.get(junction);
			String name=junction.split("_")[0];
			if(rtrn.containsKey(name)) {
				score=Math.max(rtrn.get(name), score);
			}
			rtrn.put(name, score);
		}
		
		return rtrn;
	}


	public static void writeJunction(List<String> fileNames, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		File[] files=order(fileNames);
		Map<String, Pair<Integer>>[] scores=parseJunctionIntron(files); //gene+coordinate
		Map<String, Integer> junctionCounts=parseJunctionCounts(files);
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		int counter=0;
		for(String gene: junctionCounts.keySet()) {
			String name=gene.split("_")[0];
			SingleInterval region=new SingleInterval(gene.split("_")[1]);
			double[] y=getY(scores, gene);
			double min=Statistics.min(y);
			int junctionCount=junctionCounts.get(gene);
			double[] fit=fit4PL(x,y);	
			writer.write(name+"\t"+region.toUCSC()+"\t"+junctionCount+"\t"+fit[2]+"\t"+fit[4]);
			for(int i=0; i<y.length; i++) {
				writer.write("\t"+y[i]);
			}
			writer.write("\n");
			counter++;
			if(counter%1000==0) {System.err.println(counter+" "+junctionCounts.size());}
		}
		writer.close();
	}
	
	public static void writeTotalCounts(List<String> fileNames, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		File[] files=order(fileNames);
		Map<String, Pair<Integer>>[] scores=parseTotalCounts(files); //gene+coordinate
		Map<String, Integer> junctionCounts=parseJunctionCounts(files);
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		int counter=0;
		for(String gene: junctionCounts.keySet()) {
			String name=gene.split("_")[0];
			SingleInterval region=new SingleInterval(gene.split("_")[1]);
			double[] y=getY(scores, gene);
			double min=Statistics.min(y);
			int junctionCount=junctionCounts.get(gene);
			double[] fit=fit4PL(x,y);	
			writer.write(name+"\t"+region.toUCSC()+"\t"+junctionCount+"\t"+fit[2]+"\t"+fit[4]);
			for(int i=0; i<y.length; i++) {
				writer.write("\t"+y[i]);
			}
			writer.write("\n");
			counter++;
			if(counter%1000==0) {System.err.println(counter+" "+junctionCounts.size());}
		}
		writer.close();
	}
	
	
	private static Collection<String> allJunctions(Map<String, Pair<Double>>[] scores) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(int i=0; i<scores.length; i++) {rtrn.addAll(scores[i].keySet());}
		
		return rtrn;
	}


	private static List<Double> getPerms(Map<String, int[]>[] perms, Map<String, Pair<Integer>>[] scores, double[] x, String junction) {
		List<Double> rtrn=new ArrayList<Double>();
		int numPerm=10;
		for(int i=0; i<numPerm; i++) {
			double[] y=getY(scores, junction, perms, i);
			double[] fit=fit4PL(x,y);
			if(fit[2]!=-1) {
				rtrn.add(fit[2]);
			}
			//writer.write("\t"+fit[4]+"\t"+fit[2]);
		}
		return rtrn;
		
	}

	
	

	

	private static SingleInterval coordinates(Collection<String> junctions) {
		int start=0;
		int end=0;
		String chr="";
		boolean started=false;
		
		for(String j: junctions) {
			SingleInterval region=new SingleInterval(j.split("_")[1]);
			if(!started) {
				started=true;
				start=region.getReferenceStartPosition();
				end=region.getReferenceEndPosition();
				chr=region.getReferenceName();
			}
			else {
				start=Math.min(start, region.getReferenceStartPosition());
				end=Math.max(end, region.getReferenceEndPosition());
			}
		}
		return new SingleInterval(chr, start, end);
	}


	private static Map<String, int[]>[] parsePerms(File[] files, int numPerm) throws NumberFormatException, IOException {
		Map<String, int[]>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			rtrn[i]=parsePermScore(files[i], numPerm);
		}
		
		return rtrn;
	}


	private static double sum(Collection<String> junctions, Map<String, Double> junctionCounts) {
		double sumCount=0;
		for(String j: junctions) {
			sumCount+=junctionCounts.get(j);
		}
		return sumCount;
	}


	private static double[] getYScaled(Map<String, Pair<Integer>>[] scores, Collection<String> junctions, Map<String, Double> junctionCounts) {
		double sumCount=0;
		for(String j: junctions) {
			sumCount+=junctionCounts.get(j);
		}
		
		Pair<Double>[] vals=new Pair[scores.length];
		
		for(int i=0; i<scores.length; i++) {
			Pair<Double> sum=new Pair<Double>(0.0,0.0);
			for(String junc: junctions) {
				double weight=junctionCounts.get(junc)/sumCount;
				Pair<Integer> pair=get(scores[i],junc);
				sum.setValue1(sum.getValue1()+(weight*pair.getValue1()));
				sum.setValue2(sum.getValue2()+(weight*pair.getValue2()));
			}
			vals[i]=sum;
		}
		
		double[] rtrn=new double[scores.length];
		for(int i=0; i<vals.length; i++) {
			double num=vals[i].getValue1();
			double denom=vals[i].getValue1()+vals[i].getValue2();
			rtrn[i]=num/denom;
		}
		
		return rtrn;
	}
	
	private static double[] getYScaledFractions(Map<String, Pair<Integer>>[] scores, Collection<String> junctions, Map<String, Double> junctionCounts) {
		double sumCount=0;
		for(String j: junctions) {
			sumCount+=junctionCounts.get(j);
		}
		
		double[] rtrn=new double[scores.length];
		for(int i=0; i<scores.length; i++) {
			for(String junc: junctions) {
				double weight=junctionCounts.get(junc)/sumCount;
				Pair<Integer> pair=get(scores[i],junc);
				double fraction=(double)pair.getValue1()/(double)(pair.getValue1()+pair.getValue2());
				rtrn[i]+=(weight*fraction);
			}
		}
		
		
		return rtrn;
	}
	
	
	private static double[] getYScaledFractions(Map<String, Pair<Integer>>[] scores, Collection<String> junctions, Map<String, Double> junctionCounts, Map<String, int[]>[] perms, int index) {
		double sumCount=0;
		for(String j: junctions) {
			sumCount+=junctionCounts.get(j);
		}
		
		double[] rtrn=new double[scores.length];
		for(int i=0; i<scores.length; i++) {
			for(String junc: junctions) {
				double weight=junctionCounts.get(junc)/sumCount;
				Pair<Integer> pair=get(scores[i],junc);
				double denom=(double)(pair.getValue1()+pair.getValue2());
				double num=perms[i].get(junc)[index];
				double fraction=num/denom;
				rtrn[i]+=(weight*fraction);
			}
		}
		
		
		return rtrn;
	}

	
	private static double[] getY(Map<String, Pair<Integer>>[] scores, String junction, Map<String, int[]>[] perms, int index) {
		double sumCount=0;
		
		
		double[] rtrn=new double[scores.length];
		for(int i=0; i<scores.length; i++) {
			Pair<Integer> pair=get(scores[i],junction);
			double denom=(double)(pair.getValue1()+pair.getValue2());
			double num=perms[i].get(junction)[index];
			double fraction=num/denom;
			rtrn[i]=(fraction);
			
		}
		
		
		return rtrn;
	}
	

	private static double[] getY(Map<String, Pair<Integer>>[] scores, Collection<String> junctions) {
		Pair<Integer>[] vals=new Pair[scores.length];
		
		for(int i=0; i<scores.length; i++) {
			Pair<Integer> sum=new Pair<Integer>(0,0);
			for(String junc: junctions) {
				Pair<Integer> pair=get(scores[i],junc);
				sum.setValue1(sum.getValue1()+pair.getValue1());
				sum.setValue2(sum.getValue2()+pair.getValue2());
			}
			vals[i]=sum;
		}
		
		double[] rtrn=new double[scores.length];
		for(int i=0; i<vals.length; i++) {
			double num=vals[i].getValue1();
			double denom=vals[i].getValue1()+vals[i].getValue2();
			rtrn[i]=num/denom;
		}
		
		return rtrn;
	}


	private static Pair<Integer> get(Map<String, Pair<Integer>> map, String junc) {
		if(map.containsKey(junc)) {return map.get(junc);}
		return new Pair<Integer>(0,0);
	}


	private static void add(Map<String, Collection<String>> junctionsByGene, String name, String gene) {
		if(!junctionsByGene.containsKey(name)) {junctionsByGene.put(name, new TreeSet<String>());}
		Collection<String> set=junctionsByGene.get(name);
		set.add(gene);
	}


	private static Map<String, Double> getMax(Map<String, Double> junctionCounts) {
		Map<String, List<Double>> list=new TreeMap<String, List<Double>>();
		for(String gene: junctionCounts.keySet()) {
			double count=junctionCounts.get(gene);
			String name=gene.split("_")[0];
			if(!list.containsKey(name)) {list.put(name, new ArrayList<Double>());}
			List<Double> set=list.get(name);
			set.add(count);
		}
		
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(String name: list.keySet()) {
			List<Double> vals=list.get(name);
			rtrn.put(name, Statistics.max(vals));
		}
		
		return rtrn;
	}


	private static Map<String, Double> getJunctionCount(Map<String, Pair<Integer>>[] scores) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(int i=0; i<scores.length; i++) {
			for(String gene: scores[i].keySet()) {
				int count=scores[i].get(gene).getValue1();
				double sum=0;
				if(rtrn.containsKey(gene)) {
					sum=rtrn.get(gene);
				}
				sum+=count;
				rtrn.put(gene, sum);
			}
		}
		
		return rtrn;
	}


	private static void writeTable(File[] files, Map<String, SingleInterval> coordinates,
			Map<String, Pair<Integer>>[] scores, String save, Map<String, Integer> sum, Map<String, Double> ratios,
			Map<String, Boolean> unique, Map<String, Double> junctionCounts, Map<String, String> linesByGene, Map<String, double[]>[] permScores) throws IOException {
		
		
		int numPerm=100;
		FileWriter writer=new FileWriter(save);
		
		FileWriter writerBG=new FileWriter(save+".bedgraph");
		FileWriter writerBGLow=new FileWriter(save+".low.bedgraph");
		FileWriter writerBGHigh=new FileWriter(save+".high.bedgraph");
		
		double r2cutoff=0.9;
		
		System.err.println(coordinates.size());
		//Collection<String> allGenes=new TreeSet<String>();
		
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		
		//double[] x= {0, 11.75, 16.92, 21.50, 26.92, 31.92, 47.63, 62.57, 77.50, 92.50, 122.50, 250.92};
		
		writer.write("Gene Name\t coordinates \t unique mapped junction \t junction counts \t sum of read counts \t exon sum \t intron sum\t min read count \t length \texpected ratio\ttop\tbottom\tC50\thill\tr2\tslope [0-1]\tN\taverage C50\tstdev C50\tmedian\t5%\t95%");
		
		for(int i=0; i<scores.length; i++) {
			writer.write("\t"+x[i]);
		}
		/*for(int i=0; i<scores.length; i++) {
			writer.write("\t"+x[i]);
		}
		for(int i=0; i<scores.length; i++) {
			writer.write("\t"+x[i]);
		}*/
		writer.write("\n");
		
		int counter=0;
		for(String gene: coordinates.keySet()) {
			SingleInterval coordinate=coordinates.get(gene);
		
			String line="";
			if(linesByGene.containsKey(gene)) {line=linesByGene.get(gene);}
			
			double[] y=getY(scores, gene);
			double minCount=minCount(scores, gene);
			if(minCount>10) {
			Pair<Integer> sumEI=sumExonIntron(scores, gene);
			
			/*double[] ymean=get(ratioEstimates, gene, 0);
			double[] yBottom=get(ratioEstimates, gene, 1);
			double[] yTop=get(ratioEstimates, gene, 2);*/
			
			List<Double> c50Perm=new ArrayList<Double>();
			for(int i=0; i<numPerm; i++) {
				double[] yperm=getYPerm(permScores, gene, i);
				double[] fitParms=fit4PL(x,yperm);
				if(fitParms!=null && fitParms[4]>r2cutoff && fitParms[2]<200) {c50Perm.add(fitParms[2]);}
			}
			
			//System.err.println(x.length+" "+y.length);
			
			double[] fitParms=fit4PL(x,y);
			if(fitParms!=null && fitParms[4]>r2cutoff) {
				double[] ynorm=norm(y, fitParms);
				double slope2=slope(fitParms, 0, 1);
				String name=gene.split("_")[0];
				Boolean isUnique=unique.get(gene);
				double counts=junctionCounts.get(gene);
				double ratio=ratios.get(gene);
				
				if(!line.isEmpty()) {
					writer.write(line+"\t");
				}
				
				writer.write(name+"\t"+coordinate.toUCSC()+"\t"+isUnique+"\t"+counts+"\t"+sum.get(gene)+"\t"+sumEI.getValue1()+"\t"+sumEI.getValue2()+"\t"+minCount+"\t"+coordinates.get(gene).getGenomicLength()+"\t"+ratio+"\t"+fitParms[0]+"\t"+fitParms[1]+"\t"+ fitParms[2]+"\t"+fitParms[3]+"\t"+fitParms[4]+"\t"+slope2);
				
				double mean=Statistics.mean(c50Perm);
				writer.write("\t"+c50Perm.size()+"\t"+mean);
				writer.write("\t"+Statistics.sem(c50Perm, mean));
				writer.write("\t"+Statistics.quantile(c50Perm, 0.5));
				writer.write("\t"+Statistics.quantile(c50Perm, 0.05));
				writer.write("\t"+Statistics.quantile(c50Perm, 0.95));
				
				for(int i=0; i<y.length; i++) {writer.write("\t"+y[i]);}
				//for(int i=0; i<ynorm.length; i++) {writer.write("\t"+ynorm[i]);}
				
				/*for(int i=0; i<ymean.length; i++) {writer.write("\t"+ymean[i]);}
				for(int i=0; i<yBottom.length; i++) {writer.write("\t"+yBottom[i]);}
				for(int i=0; i<yTop.length; i++) {writer.write("\t"+yTop[i]);}*/
				
				//for(Double c50: c50Perm) {writer.write("\t"+c50);}
				
				//System.out.println(gene+"\t"+derivative);
				writer.write("\n");
				
				writerBG.write(coordinate.toBedgraph(50-Statistics.quantile(c50Perm, 0.5))+"\n");
				writerBGHigh.write(coordinate.toBedgraph(50-Statistics.quantile(c50Perm, 0.95))+"\n");
				writerBGLow.write(coordinate.toBedgraph(50-Statistics.quantile(c50Perm, 0.05))+"\n");
			}
			}
			
			counter++;
			if(counter%1000 ==0) {System.err.println(counter);}
			//if(counter>100) {System.err.println("done"); break;}
		}
		
		
		writer.close();
		writerBG.close();
		writerBGLow.close();
		writerBGHigh.close();
		
	}


	private static double[] getYPerm(Map<String, double[]>[] scores, String gene, int index) {
		double[] yPerm=new double[scores.length];
		for(int i=0; i<scores.length; i++) {
			double[] vals=getVals(scores[i], gene);
			yPerm[i]=0;
			if(vals!=null) {
				yPerm[i]=vals[index];
			}
		}
		return yPerm;
	}

	private static double[] getVals(Map<String, double[]> map, String gene) {
		if(map.containsKey(gene)) {
			double[] vals=map.get(gene);
			return vals;
		}
		return null;
	}


	

	
	

	private static double[] get(Map<String, double[]>[] ratioEstimates, String gene, int index) {
		double[] rtrn=new double[ratioEstimates.length];
		
		for(int i=0; i<ratioEstimates.length; i++) {
			if(ratioEstimates[i].containsKey(gene)) {
				double[] vals=ratioEstimates[i].get(gene);
				rtrn[i]=vals[index];
			}
			else {rtrn[i]=-1;}
		}
		
		return rtrn;
	}


	private static Pair<Integer> sumExonIntron(Map<String, Pair<Integer>>[] scores, String gene) {
		Pair<Integer> rtrn=new Pair<Integer>(0,0);
		
		for(int i=0; i<scores.length; i++) {
			if(scores[i].containsKey(gene)) {
				Pair<Integer> pair=scores[i].get(gene);
				rtrn.setValue1(rtrn.getValue1()+pair.getValue1());
				rtrn.setValue2(rtrn.getValue2()+pair.getValue2());
			}
		}
		
		return rtrn;
	}


	private static Map<String, Integer> sum(Map<String, Pair<Integer>>[] scores) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(int i=0; i<scores.length; i++) {
			for(String key: scores[i].keySet()) {
				int sum=sum(scores[i].get(key));
				int count=0;
				if(rtrn.containsKey(key)) {count=rtrn.get(key);}
				count+=sum;
				rtrn.put(key, count);
			}
		}
		
		return rtrn;
	}


	private static int sum(Pair<Integer> pair) {
		return pair.getValue1()+pair.getValue2();
	}


	private static Map<String, Pair<Integer>>[] parseJunctionScores(File[] files) throws NumberFormatException, IOException {
		Map<String, Pair<Integer>>[] rtrn=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			
			rtrn[i]=parseJunctionScores(files[i]);
			
		}
		
		return rtrn;
	}


	private static Map<String, Pair<Integer>> parseJunctionScores(File file) throws NumberFormatException, IOException {
		Map<String, Pair<Integer>> rtrn=new TreeMap<String, Pair<Integer>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			
			for(int index=5; index<tokens.length; index+=3) {
				String region=tokens[index];
				int val1=Integer.parseInt(tokens[index+1]);
				int val2=Integer.parseInt(tokens[index+2]);
				Pair<Integer> pair=new Pair<Integer>(val1, val2);
				rtrn.put(region, pair);
			}
			
			
		}
		reader.close();
		
		
		return rtrn;
	}


	private static Map<String, Collection<String>> parseJunctions(File[] files) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getName());
			
			
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(files[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				String[] tokens=nextLine.split("\t");
				String name=tokens[0];
				Collection<String> regions=new TreeSet<String>();
				
				for(int index=5; index<tokens.length; index+=3) {
					regions.add(tokens[index]);
				}
				
				rtrn.put(name,  regions);
			}
			reader.close();
			
		}
		
		return rtrn;
	}
	
	
	private static Map<String, String> parse(String file) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			rtrn.put(tokens[0], nextLine);
		}
		reader.close();
		return rtrn;
	}

	
	public static void main(String[] args) throws NumberFormatException, IOException {
		if(args.length>1) {
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			Map<String, String> lineByGene=new TreeMap<String, String>();
			
			if(args.length>2) {
				lineByGene=parse(args[2]);	
			}
			run(files, save, 10, lineByGene);
		}
		else {System.err.println(usage);}
		
		
		
		
	}
	
	
	


	//TODO write bedgraph
	//TODO pull specific genes with all replicates

	static String usage="-Djava.awt.headless=true\n args[0]=files \n args[1]=save";






	
	
}

package guttmanlab.core.spidr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.distribution.HypergeometricDistributionImpl;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class PeakOverlapToEncode {

	public static void comparePeaks(File peakFile1, File peakFile2) {
		//compute hypergeometric p-value and enrichment
	}
	
	public static void computePeakEnrichment(File peakFile, File bamFile) {
		//for each peak, compute enrichment from BAM
	}
	
	
	public static void makePeakMatrix(List<File> bedgraphFiles, String save) throws IOException {
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		List<String> columns=new ArrayList<String>();
		
		for(File file: bedgraphFiles) {
			System.err.println(file.getName());
			Map<SingleInterval, Double> scores=getScore(file);
			regions.addAll(scores.keySet());
			columns.add(file.getName());
		}
		
		List<String> rows=new ArrayList<String>();
		for(SingleInterval region: regions) {rows.add(region.toUCSC(region.getOrientation()));}
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		for(File file: bedgraphFiles) {
			Map<SingleInterval, Double> scores=getScore(file);
			add(mwh, scores, file.getName());
		}
		
		mwh=filter(mwh);
		
		mwh.write(save);
	}

	private static MatrixWithHeaders filter(MatrixWithHeaders mwh) {
		Collection<String> list=new TreeSet<String>();
		for(String row: mwh.getRowNames()) {
			double[] vals=mwh.getRow(row);
			int count=count(vals);
			if(count>=2) {list.add(row);}
		}
		return mwh.submatrixByRowNames(list);
	}

	private static int count(double[] vals) {
		int counter=0;
		
		for(int i=0; i<vals.length; i++) {
			if(vals[i]>0) {counter++;}
		}
		
		return counter;
	}

	private static Map<SingleInterval, Double> getScore(File file) throws IOException {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		List<String> lines=BEDFileIO.loadLines(file.getAbsolutePath());
		for(String line: lines) {
			String[] tokens=line.split("\t");
			SingleInterval region=new SingleInterval(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
			region.setOrientation(Strand.fromString(tokens[5]));
			
			
			double score=Double.parseDouble(tokens[4]);
			
			rtrn.put(region, score);
		}
		
		return rtrn;
	}

	private static void add(MatrixWithHeaders mwh, Map<SingleInterval, Double> scores, String name) {
		for(SingleInterval region: scores.keySet()) {
			String row=region.toUCSC(region.getOrientation());
			double val=scores.get(region);
			mwh.set(row, name, val);
		}
	}
	
	private static List<File> parseFiles(String string) {
		File[] files=new File(string).listFiles();
		List<File> rtrn=new ArrayList<File>();
		
		for(int i=0; i<files.length; i++) {
			if(files[i].getName().endsWith(".bed")) {
				rtrn.add(files[i]);
			}
		}
		
		
		return rtrn;
	}
	
	public static void filterMatrix(File input, Collection<SingleInterval> regions) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			//rtrn.add(nextLine);
			if(counter==0) {System.out.println(nextLine);}
			else {
				SingleInterval region=new SingleInterval(nextLine.split("\t")[0]);
				if(regions.contains(region)) {
					double[] vals=parseDouble(nextLine);
					if(moreThan(vals,2)) {System.out.println(nextLine);}
				}
			}
			counter++;
		}
		reader.close();
	}
	
	
	private static double[] parseDouble(String nextLine) {
		String[] tokens=nextLine.split("\t");
		double[] rtrn=new double[tokens.length-1];
		
		for(int i=1; i<tokens.length; i++) {
			rtrn[i-1]=Double.parseDouble(tokens[i]);
		}
		
		return rtrn;
	}

	private static boolean moreThan(double[] vals, int min) {
		int counter=0;
		
		for(int i=0; i<vals.length; i++) {
			if(vals[i]>1) {counter++;}
		}
		
		return counter>=min;
	}
	
	private static void mergeAll(File[] files1, File[] files2) throws NumberFormatException, IOException {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		for(int i=0; i<files1.length; i++) {
			rtrn.addAll(parse(files1[i].getAbsolutePath()));
		}
		for(int i=0; i<files2.length; i++) {
			rtrn.addAll(parse(files2[i].getAbsolutePath()));
		}
		
		for(SingleInterval r: rtrn) {
			System.out.println(r.toBED());
		}
		
	}
	
	private static void hypergeometric(File[] files1, File[] files2, String universeFile, String save, Map<String, String> names) throws NumberFormatException, IOException {
		MatrixWithHeaders mwh=new MatrixWithHeaders(getNames(files1, names), getNames(files2, names));
		Collection<SingleInterval> universeSet=parse(universeFile);
		double alpha=0.05;
		
		for(int i=0; i<files1.length; i++) {
			String row=getName(files1[i].getName(), names);
			File file2=getFile(files2, row, names);
			//Collection<File> nonSelf=getNonself(files2, row, names);
			if(file2!=null) {
				double[] vals=hypergeometric(files1[i], file2, universeSet);
				System.out.println(row+"\t"+files1[i].getName()+"\t"+file2.getName()+"\t"+vals[0]+"\t"+vals[1]+"\t"+vals[2]+"\t"+vals[3]+"\t"+vals[4]);
			}
		}
		
		mwh.write(save);
		
	}
	
	private static void hypergeometricAverage(File[] files1, File[] files2, String universeFile, String save, Map<String, String> names) throws NumberFormatException, IOException {
		MatrixWithHeaders mwh=new MatrixWithHeaders(getNames(files1, names), getNames(files2, names));
		Collection<SingleInterval> universeSet=parse(universeFile);
		//double alpha=0.05;
		
		for(int i=0; i<files1.length; i++) {
			System.err.println(files1[i].getName());
			String row=getName(files1[i].getName(), names);
			File file2=getFile(files2, row, names);
			if(file2!=null) {
				double avg=average(files1[i], files2, row, names, universeSet);
				System.out.println(row+"\t"+files1[i].getName()+"\t"+avg);
			}
		}
		
		mwh.write(save);
		
	}
	
	private static double average(File file, File[] files2, String row, Map<String, String> names, Collection<SingleInterval> universeSet) throws NumberFormatException, IOException {
		List<Double> scores=new ArrayList<Double>();
		double alpha=0.01;
		for(int i=0; i<files2.length; i++) {
			File file2=files2[i];
			if(names.containsKey(files2[i].getName())) {
				String name=names.get(files2[i].getName());
				if(!name.equalsIgnoreCase(row)) {
					double[] vals=hypergeometric(file, file2, universeSet);
					double score=0;
					if(vals[3]<alpha) {score=vals[4];}
					scores.add(score);
				}
			}
		}
		return Statistics.quantile(scores, 0.5);
	}

	private static File getFile(File[] files2, String row, Map<String, String> names) {
		for(int i=0; i<files2.length; i++) {
			String name=names.get(files2[i].getName());
			if(name.equalsIgnoreCase(row)) {return files2[i];}
		}
		return null;
	}

	private static String getName(String name, Map<String, String> names) {
		if(names.containsKey(name)) {return names.get(name);}
		return name;
	}

	private static List<String> getNames(File[] files, Map<String, String> names) {
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<files.length; i++) {
			String name=names.get(files[i].getName());
			rtrn.add(name);
		}
		
		return rtrn;
	}

	private static double[] hypergeometric(File file1, File file2) throws NumberFormatException, IOException {
		Collection<SingleInterval> set1=parse(file1.getAbsolutePath());
		Collection<SingleInterval> set2=parse(file2.getAbsolutePath());
		//Collection<SingleInterval> universeSet=parse(universeFile);
		
		Collection<SingleInterval> universeSet=new TreeSet<SingleInterval>();
		universeSet.addAll(set1);
		universeSet.addAll(set2);
		
		int overlap=overlap(set1, set2);
		
		//int total=set1.size()+set2.size()-overlap;
		
		double prob1=(double)set1.size()/(double)universeSet.size();
		double prob2=(double)set2.size()/(double)universeSet.size();
		
		int numPerm=100;
		int[] randomOverlap=sample(universeSet.size(), prob1, prob2, numPerm, overlap);
		
		double percentile=Statistics.percentGreaterThan(overlap, randomOverlap);
		double expected=Math.max(1, Statistics.mean(randomOverlap));
		
		double enrich=(double)overlap/expected;
		System.out.println(file1.getName()+"\t"+file2.getName()+"\t"+universeSet.size()+"\t"+set1.size()+"\t"+set2.size()+"\t"+overlap+"\t"+percentile+"\t"+enrich);
		
		double[] rtrn= {set1.size(), set2.size(), overlap, percentile, enrich};
		/*for(int i=overlap; i<Math.min(set1.size(), set2.size()); i+=1000) {
			System.err.println(i+"\t"+(1-dist.cumulativeProbability(i)));
		}*/
		
		return rtrn;
	}
	
	private static double[] hypergeometric(File file1, File file2, Collection<SingleInterval> universeSet) throws NumberFormatException, IOException {
		Collection<SingleInterval> set1=parse(file1.getAbsolutePath());
		Collection<SingleInterval> set2=parse(file2.getAbsolutePath());
		//Collection<SingleInterval> universeSet=parse(universeFile);
		
		int overlap=overlap(set1, set2);
		
		//int total=set1.size()+set2.size()-overlap;
		
		double prob1=(double)set1.size()/(double)universeSet.size();
		double prob2=(double)set2.size()/(double)universeSet.size();
		
		int numPerm=100;
		int[] randomOverlap=sample(universeSet.size(), prob1, prob2, numPerm, overlap);
		
		double percentile=Statistics.percentGreaterThan(overlap, randomOverlap);
		double expected=Math.max(1, Statistics.mean(randomOverlap));
		
		double enrich=(double)overlap/expected;
		System.out.println(file1.getName()+"\t"+file2.getName()+"\t"+universeSet.size()+"\t"+set1.size()+"\t"+set2.size()+"\t"+overlap+"\t"+percentile+"\t"+enrich);
		
		double[] rtrn= {set1.size(), set2.size(), overlap, percentile, enrich};
		/*for(int i=overlap; i<Math.min(set1.size(), set2.size()); i+=1000) {
			System.err.println(i+"\t"+(1-dist.cumulativeProbability(i)));
		}*/
		
		return rtrn;
	}

	private static int[] sample(int size, double prob1, double prob2, int numPerm, int overlap) {
		int[] rtrn=new int[numPerm+1];
		rtrn[0]=overlap;
		for(int i=0; i< numPerm; i++) {
			rtrn[i+1]=sample(size, prob1, prob2);
			//System.err.println(i+" "+rtrn[i]);
		}
		return rtrn;
	}

	private static int sample(int size, double prob1, double prob2) {
		int counter=0;
		for(int i=0; i<size; i++) {
			double rand1=Math.random();
			double rand2=Math.random();
			if(rand1<prob1 && rand2<prob2) {counter++;}
		}
		return counter;
	}

	private static int overlap(Collection<SingleInterval> set1, Collection<SingleInterval> set2) {
		int counter=0;
		
		for(SingleInterval r: set1) {
			if(set2.contains(r)) {
				//System.out.println(r.toBED());
				counter++;
			}
		}
		
		return counter;
	}

	private static Collection<SingleInterval> parse(String input) throws NumberFormatException, IOException {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			SingleInterval r=new SingleInterval(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]), Strand.fromString(tokens[5]));
			rtrn.add(r);
		}
		reader.close();
		
		return rtrn;
	}
	
	private static Map<String, String> parseNames(String string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			rtrn.put(tokens[0], tokens[1]);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		//filterMatrix(new File(args[0]), parse(args[1]));
		
		//mergeAll(new File(args[0]).listFiles(), new File(args[1]).listFiles());
		
		//Map<String, String> names=parseNames(args[4]);
		
		//hypergeometricAverage(new File(args[0]).listFiles(), new File(args[1]).listFiles(), args[2], args[3], names);
		
		
		
		/*List<File> files=parseFiles(args[0]);
		String save=args[1];
		makePeakMatrix(files, save);*/
		
		//filterMatrix(new File(args[0]), parse(args[1]));
		
		List<String> rows=new ArrayList<String>();
		List<String> columns=new ArrayList<String>();
		
		List<String> lines=BEDFileIO.loadLines(args[0],1);
		Collection<SingleInterval> universe=BEDFileIO.loadSingleIntervalFromFile(args[1]);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			//System.err.println(line);
			rows.add(tokens[0]);
			columns.add(tokens[1]);
			double[] vals=hypergeometric(new File(tokens[0]), new File(tokens[1]), universe);
		}
		
		
		MatrixWithHeaders matrix=new MatrixWithHeaders(rows,columns);
		for(String row: rows) {
			for(String column: columns) {
				double[] vals=hypergeometric(new File(row), new File(column), universe);
				//System.err.println(row+" "+column+" "+vals[3]);
				matrix.set(row, column, vals[3]);
			}
		}
		
		matrix.write(args[2]);
		
		//makeUniverse(new File(args[0]).listFiles(), new File(args[1]).listFiles());
		
	}

	private static void makeUniverse(File[] listFiles, File[] listFiles2) throws IOException {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		for(int i=0; i<listFiles.length; i++) {
			rtrn.addAll(BEDFileIO.loadSingleIntervalFromFile(listFiles[i].getAbsolutePath()));
		}
		
		for(int i=0; i<listFiles2.length; i++) {
			rtrn.addAll(BEDFileIO.loadSingleIntervalFromFile(listFiles2[i].getAbsolutePath()));
		}
		
		for(SingleInterval region: rtrn) {System.out.println(region.toBED());}
		
	}

	

	

	
	
}

package guttmanlab.core.splicing.speckle;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.SimpleCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;

import flanagan.analysis.Regression;
import flanagan.physchem.ImmunoAssay;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class Utils {

	public static void filterJunctions(String input, String output, double fraction, int minRPKM) throws IOException {
		FileWriter writer=new FileWriter(output);
		
		List<String> lines=BEDFileIO.loadLines(input, 1);
		for(String line: lines) {
			String[] tokens=line.split("\t");
			SingleInterval region=new SingleInterval(tokens[2]);
			double p=Double.parseDouble(tokens[3]);
			double rpkm=Double.parseDouble(tokens[4]);
			double observed=Double.parseDouble(tokens[7]);
			int junctionCounts=Integer.parseInt(tokens[5]);
			if(p>=fraction && rpkm>=minRPKM && observed>=0 && junctionCounts>10) {
				List<Double> vals=getVals(tokens);
				double min=Statistics.quantile(vals,0.05);
				double max=Statistics.quantile(vals,0.95);
				
				double range=Math.abs(max-min);
				writer.write(tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+tokens[4]+"\t"+tokens[5]+"\t"+observed+"\t"+range+"\n");
				
				/*if(Math.abs(max-min)<delta) {
				
					writer.write(region.toBedgraph(observed)+"\n");
					writer.write(region.toBedgraph(Statistics.mean(vals))+"\n");
					writer.write(region.toBedgraph(Statistics.quantile(vals,0.5))+"\n");
					writer.write(region.toBedgraph(Statistics.quantile(vals,0.95))+"\n");
					writer.write(region.toBedgraph(Statistics.quantile(vals,0.05))+"\n");
				}*/
				
				
				
				
			}
			
		}
		
		writer.close();
	}
	
	
	private static Map<SingleInterval, Double> computeUniqueJunctions(String observed, String save) throws IOException {
		Map<SingleInterval, Integer> junctionCounts=parseJunctionCounts(observed);
		Map<String, Collection<SingleInterval>> junctionsByGene=parseJunctions(observed);
		
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		FileWriter writer=new FileWriter(save);
		
		for(String gene: junctionsByGene.keySet()) {
			
			Collection<SingleInterval> junctions=junctionsByGene.get(gene);
			for(SingleInterval junction: junctions) {
				double overlaps=fractionOfOverlappers(junction, junctions, junctionCounts);
				rtrn.put(junction, overlaps);
				writer.write(junction+"\n");
			}
		}
		writer.close();
		return rtrn;
	}
	
	private static double fractionOfOverlappers(SingleInterval junction1, Collection<SingleInterval> junctions, Map<SingleInterval, Integer> junctionCounts) {
		double o=junctionCounts.get(junction1);
		
		Collection<SingleInterval> overlappers=new TreeSet<SingleInterval>();
		
		for(SingleInterval junction2: junctions) {
			if(!junction1.equals(junction2)) {
				if(junction1.overlaps(junction2) || junction2.overlaps(junction1)) {
					overlappers.add(junction2);
				}
			}
		}
		
		System.out.println(junction1.toUCSC()+" "+overlappers.size());
		
		double sum=o;
		for(SingleInterval val: overlappers) {
			sum+=junctionCounts.get(val);
		}
		
		
		return o/sum;
	}



	private static Map<String, Collection<SingleInterval>> parseJunctions(String file) throws IOException {
		Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0];
			SingleInterval junction=new SingleInterval(tokens[1]);
			if(!rtrn.containsKey(name)) {rtrn.put(name, new TreeSet<SingleInterval>());}
			Collection<SingleInterval> list=rtrn.get(name);
			list.add(junction);
		}
		
		reader.close();
		return rtrn;
	}
	
	private static Map<SingleInterval, Integer> parseJunctionCounts(String file) throws IOException {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			//String name=tokens[0];
			SingleInterval junction=new SingleInterval(tokens[1]);
			Integer count=Integer.parseInt(tokens[4]);
			rtrn.put(junction, count);
		}
		
		reader.close();
		return rtrn;
	}
	
	
	public static void writeJunctions(String input, String output, double fraction) throws IOException {
		FileWriter writer=new FileWriter(output);
		
		List<String> lines=BEDFileIO.loadLines(input, 1);
		for(String line: lines) {
			String[] tokens=line.split("\t");
			SingleInterval region=new SingleInterval(tokens[1]);
			int count=Integer.parseInt(tokens[4]);
			double p=Double.parseDouble(tokens[2]);
			String name=tokens[0];
			if(p>=fraction) {
				//writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\tO="+count+"\n");
				writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+name+"\n");
			}
		}
		
		writer.close();
	}
	
	private static List<Double> getVals(String[] tokens) {
		List<Double> rtrn=new ArrayList<Double>();
		
		for(int i=8; i<tokens.length; i++) {
			double val=Double.parseDouble(tokens[i]);
			if(val>=0) {
				rtrn.add(val);
			}
		}
		
		return rtrn;
	}
	
	private static void filterRPKM(String string, double filter) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			double rpkm=Double.parseDouble(tokens[5]);
			if(rpkm>filter) {System.out.println(line);}
		}
		
	}
	
	private static void writeBEDgraph(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			double rpkm=Double.parseDouble(tokens[5]);
			System.out.println(tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+tokens[5]);
		}
		
	}
	
	private static void getJunctions(String gtfFile, String save) throws IOException {
		GTFToJunctions gtf=new GTFToJunctions(new File(gtfFile));
		
		Map<String, Collection<Gene>> geneToJunctions=gtf.getGeneToJunctions();
		
		FileWriter writer=new FileWriter(save);
		
		for(String gene: geneToJunctions.keySet()) {
			Collection<Gene> junctions=geneToJunctions.get(gene);
			for(Gene junction: junctions) {
				junction.setName(gene);
				writer.write(junction+"\n");
			}
		}
		
		writer.close();
	}
	
	private static void slidingWindow(String file, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Map<String, IntervalTree<Double>> trees=new TreeMap<String, IntervalTree<Double>>();
		
		Collection<SingleInterval> junctions=new TreeSet<SingleInterval>();
		
		List<String> lines=BEDFileIO.loadLines(file,1);
		for(String line: lines){
			String[] tokens=line.split("\t");
			//String name=tokens[0];
			SingleInterval junction=getJunction(tokens);
			double c50=getC50(tokens);
			double rpkm=getRPKM(tokens);
			double fraction=getFractionOfSpliced(tokens);
			if(fraction>0.75 && rpkm>5) {
				junctions.add(junction);
				add(junction, c50, trees);
			}
			
		}
		
		
		for(SingleInterval junction: junctions) {
			List<Double> vals=getWindows(junction, trees);
			writer.write(junction.toBedgraph(32-Statistics.quantile(vals, 0.5))+"\n");
		}
		
		
		writer.close();
	}
	
	
	private static void averagePerGene(String file, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		//Map<String, IntervalTree<Double>> trees=new TreeMap<String, IntervalTree<Double>>();
		
		Collection<SingleInterval> junctions=new TreeSet<SingleInterval>();
		
		Map<String, List<Double>> valsByGene=new TreeMap<String, List<Double>>();
		
		List<String> lines=BEDFileIO.loadLines(file,1);
		for(String line: lines){
			String[] tokens=line.split("\t");
			//String name=tokens[0];
			String gene=tokens[0];
			SingleInterval junction=getJunction(tokens);
			junction.setName(gene);
				double c50=getC50(tokens);
				double rpkm=getRPKM(tokens);
				double fraction=getFractionOfSpliced(tokens);
				List<Double> perms=getVals(tokens);
				double range=Math.abs(Statistics.quantile(perms, 0.95)-Statistics.quantile(perms, 0.05));
				if(fraction>0.75 && rpkm>1 && c50>=0 && range<10) {
					if(!valsByGene.containsKey(gene)) {valsByGene.put(gene, new ArrayList<Double>());}
					List<Double> list=valsByGene.get(gene);
					list.add(c50);
					System.out.println(gene+"\t"+rpkm+"\t"+junction.toBedgraph(c50));
					junctions.add(junction);
				}
			
		}
		
		
		for(String gene: valsByGene.keySet()) {
			List<Double> vals=valsByGene.get(gene);
			writer.write(gene+"\t"+Statistics.mean(vals)+"\t"+Statistics.quantile(vals, 0.5)+"\t"+Statistics.max(vals)+"\t"+Statistics.min(vals)+"\n");
		}
		
		for(SingleInterval junction: junctions) {
			List<Double> vals=valsByGene.get(junction.getName());
			double max=Statistics.min(vals);
			//System.out.println(junction.toBedgraph(1.0/max));
		}
		
		
		writer.close();
	}

	private static boolean has(Map<String, IntervalTree<SingleInterval>> trees, SingleInterval junction) {
		if(trees.containsKey(junction.getReferenceName())) {
			return trees.get(junction.getReferenceName()).hasOverlappers(junction.getReferenceStartPosition(), junction.getReferenceEndPosition());
		}
		return false;
	}


	private static double getFractionOfSpliced(String[] tokens) {
		return Double.parseDouble(tokens[3]);
	}


	private static double getRPKM(String[] tokens) {
		return Double.parseDouble(tokens[4]);
	}


	private static List<Double> getWindows(SingleInterval junction, Map<String, IntervalTree<Double>> trees) {
		List<Double> rtrn=new ArrayList<Double>();
		
		int start=junction.getReferenceStartPosition();
		int end=junction.getReferenceEndPosition();
		String chr=junction.getReferenceName();
		if(trees.containsKey(chr)) {
			IntervalTree<Double> tree=trees.get(chr);
			double v1=tree.findValue(start, end);
			double v2=-1;
			double v3=-1;
			if(tree.hasValueBeforeInterval(start, end)) {
				v2=tree.getValueBeforeInterval(start, end);
			}
			if(tree.hasValueAfterInterval(start, end)) {
				v3=tree.getValueAfterInterval(start, end);
			}
			if(v1>-1) {rtrn.add(v1);}
			if(v2>-1) {rtrn.add(v2);}
			if(v3>-1) {rtrn.add(v3);}
			//System.err.println(junction.toUCSC()+"\t"+v1+"\t"+v2+"\t"+v3);
		}
		
		return rtrn;
	}


	private static SingleInterval getJunction(String[] tokens) {
		return new SingleInterval(tokens[2]);
	}


	private static double getC50(String[] tokens) {
		return Double.parseDouble(tokens[7]);
	}


	private static void add(SingleInterval junction, double c50, Map<String, IntervalTree<Double>> trees) {
		String chr=junction.getReferenceName();
		if(!trees.containsKey(chr)) {
			trees.put(chr, new IntervalTree<Double>());
		}
		
		IntervalTree<Double> tree=trees.get(chr);
		tree.put(junction.getReferenceStartPosition(), junction.getReferenceEndPosition(), c50);
		
	}
	
	private static void countBams(File[] listFiles) {
		for(int i=0; i<listFiles.length; i++) {
			File file=listFiles[i];
			if(file.getName().endsWith(".bam")) {
				int count=count(file);
				System.out.println(file.getName()+"\t"+count);
			}
		}
		
	}
	
	
	private static void countFastqs(File[] listFiles) {
		for(int i=0; i<listFiles.length; i++) {
			File files=listFiles[i];
			if(files.isDirectory()) {
				File file=getFastq(files.listFiles());
				if(file!=null) {
					System.err.println(file.getAbsolutePath());
					int count=countFastq(file);
					System.out.println(file.getAbsolutePath()+"\t"+count);
				}
			}
		}
		
	}


	private static File getFastq(File[] listFiles) {
		File file=null;
		
		for(int i=0; i< listFiles.length; i++) {
			if(listFiles[i].getName().endsWith("fastq.gz")) {return listFiles[i];}
		}
		
		return file;
	}


	private static int countFastq(File file) {
		FastqReader parser=new FastqReader(file);
		Iterator<FastqRecord> iter=parser.iterator();
		int counter=0;
		while(iter.hasNext()) {
			iter.next();
			counter++;
		}
		parser.close();
		return counter;
	}


	private static int count(File bam) {
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(!record.getDuplicateReadFlag()) {
				counter++;
			}
		}
	
		reader.close();
		reads.close();
		return counter;
	}

	
	private static void downsample(File[] listFiles, String string, String saveDir) throws IOException {
		Map<String, Double> proportions=parseProportion(string);
		
		File[] files=order(listFiles);
		Map<SingleInterval, Integer>[] maps=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			File file=files[i];
			String name=file.getName();
			double p=proportions.get(name);
			System.err.println(file.getName()+"\t"+p);
			String save=saveDir+"/"+file.getName()+".downsample.bam";
			maps[i]=downsample(file, p, save);
		}
		
		String save=saveDir+"/intronCounts.matrix";
		write(maps, save);
	}
	
	
	private static void downsampleFastq(File fastq1, File fastq2, double proportion, String saveDir) throws IOException {
		String name=fastq1.getParentFile().getName();
		System.err.println(name);
		
		int numSamples=10;
		
		sample(fastq1, fastq2, proportion, saveDir, name, numSamples);
		
		/*for(int i=0; i<numSamples; i++) {
			String save1=saveDir+"/"+name+".s"+i+".R1.fq";
			String save2=saveDir+"/"+name+".s"+i+".R2.fq";
			
			sample(fastq1, fastq2, proportion, save1, save2);
			
			System.err.println(save1+" "+save2);
		}*/
		
		
		
	}
	
	
	private static void sample(File fastq1, File fastq2, double proportion, String saveDir, String name, int numSamples) throws IOException {
		FastqReader reader1=new FastqReader(fastq1);
		FastqReader reader2=new FastqReader(fastq2);
		
		Pair<FileWriter>[] writers=makeWriters(numSamples, saveDir, name);
		
		Iterator<FastqRecord> iter1=reader1.iterator();
		Iterator<FastqRecord> iter2=reader2.iterator();
		
		int counter=0;
		while(iter1.hasNext()) {
			FastqRecord record1=iter1.next();
			FastqRecord record2=iter2.next();
			
			double[] rand=getRandom(numSamples);
			for(int i=0; i<rand.length; i++) {
				if(rand[i]<proportion) {
					write(writers[i], record1, record2);
				}
			}
			
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
			
		}
		
		reader1.close();
		reader2.close();
		close(writers);
	}
	
	
	private static void close(Pair<FileWriter>[] writers) throws IOException {
		for(int i=0; i<writers.length; i++) {
			writers[i].getValue1().close();
			writers[i].getValue2().close();
		}
		
	}


	private static Pair<FileWriter>[] makeWriters(int numSamples, String saveDir, String name) throws IOException {
		Pair<FileWriter>[] rtrn=new Pair[numSamples];
		
		for(int i=0; i<numSamples; i++) {
			FileWriter writer1=new FileWriter(saveDir+"/"+name+".s"+i+".R1.fq");
			FileWriter writer2=new FileWriter(saveDir+"/"+name+".s"+i+".R2.fq");
			rtrn[i]=new Pair<FileWriter>(writer1, writer2);
		}
		
		return rtrn;
	}


	private static void write(Pair<FileWriter> pair, FastqRecord record1, FastqRecord record2) throws IOException {
		write(record1, record2, pair.getValue1(), pair.getValue2());
	}


	private static double[] getRandom(int numSamples) {
		double[] rtrn=new double[numSamples];
		
		for(int i=0; i<numSamples; i++) {
			rtrn[i]=Math.random();
		}
		
		return rtrn;
	}


	private static void sample(File fastq1, File fastq2, double proportion, String save1, String save2) throws IOException {
		FastqReader reader1=new FastqReader(fastq1);
		FastqReader reader2=new FastqReader(fastq2);
		
		FileWriter writer1=new FileWriter(save1);
		FileWriter writer2=new FileWriter(save2);
		
		Iterator<FastqRecord> iter1=reader1.iterator();
		Iterator<FastqRecord> iter2=reader2.iterator();
		
		while(iter1.hasNext()) {
			FastqRecord record1=iter1.next();
			FastqRecord record2=iter2.next();
			
			double rand=Math.random();
			if(rand<proportion) {
				write(record1, record2, writer1, writer2);
			}
			
		}
		
		reader1.close();
		reader2.close();
		writer1.close();
		writer2.close();
	}


	private static void write(FastqRecord record1, FastqRecord record2, FileWriter writer1, FileWriter writer2) throws IOException {
		write(record1, writer1);
		write(record2, writer2);
	}


	private static void write(FastqRecord record1, FileWriter writer1) throws IOException {
		writer1.write("@"+record1.getReadHeader()+"\n");
		writer1.write(record1.getReadString()+"\n");
		writer1.write("+\n");
		writer1.write(record1.getBaseQualityString()+"\n");
		
	}


	private static void spliceCounts(File[] listFiles, String save) throws IOException {
		//Map<String, Double> proportions=parseProportion(string);
		
		File[] files=order(listFiles);
		Map<SingleInterval, Integer>[] maps=new Map[files.length];
		
		for(int i=0; i<files.length; i++) {
			File file=files[i];
			String name=file.getName();
			//double p=proportions.get(name);
			System.err.println(file.getName());
			//String save=saveDir+"/"+file.getName()+".downsample.bam";
			maps[i]=spliceCounts(file);
		}
		
		//String save=saveDir+"/intronCounts.matrix";
		write(maps, save);
	}

	private static void write(Map<SingleInterval, Integer>[] maps, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		for(int i=0; i<maps.length; i++) {regions.addAll(maps[i].keySet());}
		
		for(SingleInterval region: regions) {
			System.out.println(region.toShortBED());
			writer.write(region.toUCSC());
			for(int i=0; i<maps.length; i++) {
				int score=get(maps[i], region);
				writer.write("\t"+score);
			}
			writer.write("\n");
		}
		
		writer.close();
	}


	private static int get(Map<SingleInterval, Integer> map, SingleInterval region) {
		int val=0;
		if(map.containsKey(region)) {val=map.get(region);}
		return val;
	}


	private static File[] order(File[] allFiles) {
		List<File> files=new ArrayList<File>();
		for(int i=0; i<allFiles.length; i++) {
			File file=allFiles[i];
			if(file.getName().endsWith(".bam")) {files.add(file);}
		}
		
		
		File[] rtrn=new File[files.size()];
		
		int counter=0;
		for(File f: files) {
			rtrn[counter]=f;
			counter++;
		}
		
		return order2(rtrn);
	}
	
	private static File[] order2(File[] files) {
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
	

	private static Map<String, Double> parseProportion(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			rtrn.put(tokens[0], Double.parseDouble(tokens[2]));
		}
		
		return rtrn;
	}

	
	private static Map<SingleInterval, Integer> spliceCounts(File bam) {
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			//if(!record.getDuplicateReadFlag()) {
				//double rand=Math.random();
				if(SAMFragment.isSpliced(record)) {
					add(record, scores);
				}
				counter++;
			//}
		}
		
		//System.err.println(bam.getName()+"\t"+counter);	
		//print(bam.getName(), scores);
		reader.close();
		reads.close();
		//writer.close();
		return scores;
	}

	private static Map<SingleInterval, Integer> downsample(File bam, double p, String save) {
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save));
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(!record.getDuplicateReadFlag()) {
				double rand=Math.random();
				if(rand<p) {
					if(SAMFragment.isSpliced(record)) {
						writer.addAlignment(record);
						add(record, scores);
					}
					counter++;
				}
			}
		}
		
		System.err.println(bam.getName()+"\t"+counter);	
		//print(bam.getName(), scores);
		reader.close();
		reads.close();
		writer.close();
		return scores;
	}


	private static void print(String name, Map<SingleInterval, Integer> scores) {
		for(SingleInterval region: scores.keySet()) {
			int score=scores.get(region);
			System.out.println(name+"\t"+region.toBedgraph(score));
		}
		
	}


	private static void add(SAMRecord record, Map<SingleInterval, Integer> scores) {
		
		SAMFragment frag=new SAMFragment(SAMFragment.compressCIGAR(record));
		//SAMFragment frag=new SAMFragment(record);
		
		Collection<Annotation> introns=frag.getIntrons();
		for(Annotation intron: introns) {
			int count=0;
			SingleInterval intronSI=intron.getSingleInterval();
			if(scores.containsKey(intronSI)) {
				count=scores.get(intronSI);
			}
			count++;
			scores.put(intronSI, count);
		}	
	}
	
	
	private static void computeSlopes(String input, String save) throws IOException {
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String geneName=tokens[0];
			String intron=tokens[1];
			String rpkm=tokens[2];
			String junctionCount=tokens[3];
			//double[] vals=vals(tokens);
			if(Integer.parseInt(junctionCount)>10 && Double.parseDouble(tokens[5])>-1.0 && Double.parseDouble(rpkm)>=1.0) {
				double[] scale=parseVals(x, tokens);
				if(Statistics.min(scale)>=0) {
					//double[] scale=scale(vals, weightFactors);
					double[] fit=fit4PL(x,scale);
					if(fit[4]>0.9) {
						write(writer, intron, scale, fit, geneName, x, rpkm, junctionCount);
					}
				}
			}
			
			counter++;
			if(counter%1000==0) {System.err.println(counter);}
		}
		reader.close();
		writer.close();
		
	}
	
	
	private static void compute4PL(String input) throws IOException {
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			if(counter>0) {
				String[] tokens=nextLine.split("\t");
				String geneName=tokens[0];
				String name=tokens[1];
				
				double[] y=vals(tokens, 2);
				
				double[] fit=fit4PL(x,y);
				System.out.println(geneName+"\t"+name+"\t"+fit[0]+"\t"+fit[1]+"\t"+fit[2]+"\t"+fit[3]+"\t"+fit[4]);
				
				
			}
			
			counter++;
			if(counter%1000==0) {System.err.println(counter);}
		}
		reader.close();
		
	}
	
	
	private static void localMaximum(String input) throws IOException {
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			if(counter>0) {
				String[] tokens=nextLine.split("\t");
				String geneName=tokens[0];
				String name=tokens[1];
				
				double[] y=vals(tokens, 2);
				
				double max=Statistics.max(y);
				
				System.out.print(geneName+"\t"+name);
				
				boolean started=false;
				for(int i=0; i<y.length; i++) {
					if(y[i]>=max) {started=true;}
					if(started) {System.out.print("\t"+y[i]);}
					else {System.out.print("\t");}
				}
				
				System.out.println();
				
			}
			
			counter++;
			if(counter%1000==0) {System.err.println(counter);}
		}
		reader.close();
		
	}
	
	
	private static void ratioAtMax(String inputSR, String inputIntron) throws IOException {
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		
		Map<String, Integer> genePositionMax=getPositionMax(inputIntron);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(inputSR)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			if(counter>0) {
				String[] tokens=nextLine.split("\t");
				String geneName=tokens[0];
				String name=tokens[1];
				
				double[] y=vals(tokens, 2);
				int maxPosition=genePositionMax.get(geneName);
				
				double ratio=y[maxPosition];
				
				System.out.println(geneName+"\t"+name+"\t"+maxPosition+"\t"+ratio);
				//double max=Statistics.max(y);
				
				
				
			}
			
			counter++;
			if(counter%1000==0) {System.err.println(counter);}
		}
		reader.close();
		
	}
	
	private static Map<String, Integer> getPositionMax(String input) throws IOException {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			if(counter>0) {
				String[] tokens=nextLine.split("\t");
				String geneName=tokens[0];
				//String name=tokens[1];
				
				double[] y=vals(tokens, 2);
				int maxPosition=maxPosition(y);
				
				rtrn.put(geneName, maxPosition);
				
				
			}
			
			counter++;
			if(counter%1000==0) {System.err.println(counter);}
		}
		reader.close();
		
		return rtrn;
	}


	private static int maxPosition(double[] y) {
		int pos=0;
		double max=Statistics.max(y);
		
		for(int i=0; i<y.length; i++) {
			if(y[i]==max) {return i;}
		}
		
		return pos;
	}


	private static double[] vals(String[] tokens, int index) {
		double[] rtrn=new double[tokens.length-index];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=Double.parseDouble(tokens[index+i]);
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

	
	private static double[] parseVals(double[] x, String[] tokens) {
		double[] rtrn=new double[x.length];
		
		for(int i=0; i<x.length; i++) {
			rtrn[i]=Double.parseDouble(tokens[tokens.length-rtrn.length+i]);
		}
		
		return rtrn;
	}

	
	private static void sumByGene(String input, double[] weights, GTFToJunctions gtf) throws IOException {
		Map<String, double[]> scores=new TreeMap<String, double[]>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String intron=tokens[0];
			double[] vals=vals(tokens);
			Collection<String> genes=gtf.getGene(new SingleInterval(intron));
			if(genes.isEmpty()) {System.err.println("Skipped "+intron);}
			
			for(String gene: genes) {
				if(!scores.containsKey(gene)) {scores.put(gene, new double[vals.length]);}
				double[] oldVals=scores.get(gene);
				for(int i=0; i<oldVals.length; i++) {
					oldVals[i]=oldVals[i]+vals[i];
				}
				scores.put(gene, oldVals);
			}
			
			counter++;
			if(counter%1000==0) {System.err.println(counter);}
		}
		reader.close();
		
		write(scores, weights);
	}

	private static void write(Map<String, double[]> scores, double[] weights) {
		for(String gene: scores.keySet()) {
			double[] vals=scores.get(gene);
			double[] scaled=mean(vals, weights);
			double max=max(vals, weights);
			if(max>100) {
				System.out.print(gene);
				for(int i=0; i<scaled.length; i++) {
					System.out.print("\t"+scaled[i]);
				}
				System.out.println();
			}
		}
		
	}


	private static void downsampleCounts(String input, double[] weightFactors, String save, GTFToJunctions gtf) throws IOException {
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String intron=tokens[0];
			double[] vals=vals(tokens);
			double[] scale=mean(vals, weightFactors);
			//double[] scale=scale(vals, weightFactors);
			if(Statistics.sum(vals)>100) {
				Collection<String> genes=gtf.getGene(new SingleInterval(intron));
				if(genes.isEmpty()) {System.err.println("Skipped "+intron);}
				double[] fit=fit4PL(x,scale);
				//double derivative=derivative(x, scale);
				if(fit[2]>-1 && fit[4]>0.9) {
					write(writer, intron, scale, fit, genes, x);
				}
				//if(fit[3]!=-1) {System.out.println(new SingleInterval(intron).toBedgraph(fit[3]));}
			}
			counter++;
			if(counter%1000==0) {System.err.println(counter);}
		}
		reader.close();
		writer.close();
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
	
	private static double[] derivatives(double[] x, double[] y) {
		Collection<Double> vals=new ArrayList<Double>();
		for(int i=1; i<x.length; i++) {
			double num=y[i]-y[i-1];
			double denom=x[i]-x[i-1];
			double rate=num/denom;
			vals.add(rate);
		}
		
		double[] rtrn=new double[vals.size()];
		int i=0;
		for(Double val: vals) {
			rtrn[i]=val;
			i++;
		}
		return rtrn;
	}

	private static double[] mean(double[] vals, double[] weightFactors) {
		double[] rtrn=new double[vals.length];
		
		for(int i=0; i<vals.length; i++) {
			double observed=vals[i];
			double p=weightFactors[i];
			rtrn[i]=observed*p;
		}
		
		double max=Statistics.max(rtrn);
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=rtrn[i]/max;
		}
		
		return rtrn;
	}
	
	
	private static double max(double[] vals, double[] weightFactors) {
		double[] rtrn=new double[vals.length];
		
		for(int i=0; i<vals.length; i++) {
			double observed=vals[i];
			double p=weightFactors[i];
			rtrn[i]=observed*p;
		}
		
		double max=Statistics.max(rtrn);
		
		return max;
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
	
	private static double[] scale(double[] vals, double[] weightFactors) {
		double[] rtrn=new double[vals.length];
		
		for(int i=0; i<vals.length; i++) {
			double observed=vals[i];
			double p=weightFactors[i];
			rtrn[i]=sample(observed, p);
		}
		
		return rtrn;
	}


	private static double sample(double observed, double p) {
		double count=0;
		for(int i=0; i<observed; i++) {
			double rand=Math.random();
			if(rand<p) {count++;}
		}
		return count;
	}


	private static double[] vals(String[] tokens) {
		double[] rtrn=new double[tokens.length-1];
		for(int i=1; i<tokens.length; i++) {
			rtrn[i-1]=Double.parseDouble(tokens[i]);
		}
		return rtrn;
	}

	
	private static void write(FileWriter writer, String intron, double[] y, double[] fit, String gene, double[] x, String rpkm, String jc) throws IOException {
		double slope=slope(fit);
		double slope2=slope(fit, 0,1);
		double derivative=derivative(x,y);
		double[] slopes=derivatives(x,y);
		double[] scaled=norm(y, fit);
		double derivative2=derivative(x, scaled);
		Arrays.sort(slopes);
		writer.write(gene+"\t"+intron+"\t"+rpkm+"\t"+jc+"\t"+fit[2]+"\t"+fit[3]+"\t"+fit[4]+"\t"+derivative+"\t"+slopes[slopes.length-2]+"\t"+derivative2+"\t"+slope+"\t"+slope2);
		for(int i=0; i<y.length; i++) {writer.write("\t"+y[i]);}
		for(int i=0; i<scaled.length; i++) {writer.write("\t"+scaled[i]);}
		writer.write("\n");
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


	private static void write(FileWriter writer, String intron, double[] y, double[] fit, Collection<String> genes, double[] x) throws IOException {
		double slope=slope(fit);
		double slope2=slope(fit,0,1);
		double derivative=derivative(x,y);
		//double[] slopes=derivatives(x,y);
		//Arrays.sort(slopes);
		
		double[] scaled=norm(y, fit);
		double derivative2=derivative(x,scaled);
		
		for(String gene: genes) {
			writer.write(gene+"\t"+intron+"\t"+fit[2]+"\t"+fit[3]+"\t"+fit[4]+"\t"+derivative+"\t"+derivative2+"\t"+slope+"\t"+slope2);
			//for(int i=0; i<slopes.length; i++) {writer.write("\t"+slopes[i]);}
			for(int i=0; i<y.length; i++) {writer.write("\t"+y[i]);}
			for(int i=0; i<scaled.length; i++) {writer.write("\t"+scaled[i]);}
			writer.write("\n");
		}
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

	private static double[] weightFactors(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		double[] rtrn=new double[lines.size()];
		int i=0;
		for(String line: lines) {
			rtrn[i]=Double.parseDouble(line.split("\t")[2]);
			i++;
		}
		return rtrn;
	}
	
	
	private static void intronCounts(File[] bamFiles, GTFToJunctions gtf) throws IOException {
		File[] files=order(bamFiles);
		Map<Gene, JunctionScore>[] maps=new Map[files.length];
		Collection<Gene> junctions=new TreeSet<Gene>();
		Map<String, IntervalTree<Gene>> geneTree=gtf.getJunctionTree();
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i]);
			File bam=files[i];
			maps[i]=readsOverlappingSpliceJunction(bam, geneTree);
			junctions.addAll(maps[i].keySet());
		}
		
		for(Gene junction: junctions) {
			SingleInterval junctionSI=junction.getIntrons().iterator().next().getSingleInterval();
			System.out.print(junction.getName()+"\t"+junctionSI.toUCSC());
			for(int i=0; i<maps.length; i++) {
				int count=0;
				if(maps[i].containsKey(junction)) {
					JunctionScore score=maps[i].get(junction);
					count=score.getExonJunctionCounts();
					//count=score.getIntronCounts();
					//count=score.getAllExonCounts();
				}
				System.out.print("\t"+count);
			}
			System.out.println();
		}
		
	}
	
	
	private static Map<Gene, JunctionScore> readsOverlappingSpliceJunction(File bamFile, Map<String, IntervalTree<Gene>> geneTree) {
		Map<Gene, JunctionScore> rtrn=new TreeMap<Gene, JunctionScore>();
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(!record.getDuplicateReadFlag()) {
				SAMFragment read=new SAMFragment(record);
				Collection<Gene> genes=findGene(read, geneTree);
				Collection<Gene> overlappingGenes=overlapsGene(read, genes);
				for(Gene g: overlappingGenes) {
					if(!rtrn.containsKey(g)) {
						rtrn.put(g, new JunctionScore(g));
					}
					JunctionScore score=rtrn.get(g);
					int[] pair=assign(g, read); 
					score.increment(pair[0], pair[1], pair[2]);
				}
			}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		
		reader.close();
		reads.close();
		return rtrn;
	}
	
	
	private static String assignState(SAMFragment read, Gene gene) {
		//If read is spliced, check if splice junction matches gene junction, if so --> spliced
		if(read.isSpliced()) {
			if(junctionMatch(read, gene)) {return "splicedJunction";}
			else {return "amb";}
		}
		
		for(Annotation intron: gene.getIntrons()) {
			if(read.overlaps(intron)) {
				//TODO Trim read by 3 bases
				return "unspliced";
			}
		}
		//return "spliced";
		return "splicedNoJunction";
	}

	private static boolean junctionMatch(SAMFragment read, Gene gene) {
		if(read.getOrientation().equals(gene.getOrientation())) {
			Collection<SingleInterval> readJunctions=getJunctions(read);
			Collection<SingleInterval> geneJunctions=getJunctions(gene);
			for(SingleInterval readJunction: readJunctions) {
				if(geneJunctions.contains(readJunction)){return true;}
			}	
		}
		return false;
	}

	
	private static int[] assign(Gene g, SAMFragment read) {
		
		int[] rtrn=new int[3];
		
		String state=assignState(read, g);
		if(state.equals("splicedJunction")) {
			rtrn[0]++;
		}
		if(state.equals("splicedNoJunction")) {rtrn[1]++;}
		if(state.equals("unspliced")) {
			rtrn[2]++;
		}
		//Pair<Integer> rtrn=new Pair<Integer>(exonCount, intronCount);
		
		return rtrn;
		
		
	}

	
	private static Collection<Gene> overlapsGene(SAMFragment read, Collection<Gene> genes) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		for(Gene gene: genes){
			if(overlaps(gene, read)) {rtrn.add(gene);}
		}
		return rtrn;
		
	}
	
	private static Collection<SingleInterval> getJunctions(Annotation gene) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		Collection<Annotation> blocks=gene.getIntrons();
		
		for(Annotation intron: blocks) {
			SingleInterval intronSI=new SingleInterval(intron.getReferenceName(), intron.getReferenceStartPosition(), intron.getReferenceEndPosition());
			intronSI.setOrientation(gene.getOrientation());
			rtrn.add(intronSI);
		}
		
		return rtrn;
	}
	
	private static boolean overlaps(Gene gene, SAMFragment read) {
		if(read.getOrientation().equals(gene.getOrientation())) {
			
			Collection<SingleInterval> readJunctions=getJunctions(read);
			Collection<SingleInterval> geneJunctions=getJunctions(gene);
			for(SingleInterval readJunction: readJunctions) {
				if(geneJunctions.contains(readJunction)){return true;}
			}
			
			if(read.getReferenceName().equals(gene.getReferenceName()) && read.getReferenceStartPosition()>=gene.getReferenceStartPosition() && read.getReferenceEndPosition()<=gene.getReferenceEndPosition()) {
				return true;
			}
		}
		return false;
	}

	
	private static Collection<Gene> findGene(Annotation fragment, Map<String, IntervalTree<Gene>> geneTree) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		IntervalTree<Gene> tree=geneTree.get(fragment.getReferenceName());
		if(tree!=null){
			Iterator<Gene> genes=tree.overlappingValueIterator(fragment.getReferenceStartPosition(), fragment.getReferenceEndPosition());
			while(genes.hasNext()){
				Gene gene=genes.next();
				if(fragment.getOrientation().equals(gene.getOrientation())){
					if(gene.getNumberOfBlocks()>1){
						rtrn.add(gene);
					}
				}
			}
		}
		return rtrn;
		
	}


	private static Collection<SingleInterval> parseIntrons(String introns) throws IOException {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(introns)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			rtrn.add(new SingleInterval(nextLine.split("\t")[0]));
		}
		reader.close();
		return rtrn;
	}

	private static void scaleNorm(String input, double minVal) throws IOException {
		System.out.println("gene\tjunction\t0\t10\t15\t20\t25\t30\t45\t60\t75\t90\t120\t240");
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String gene=tokens[0];
			String junction=tokens[1];
			double[] vals=getVals2(tokens);
			double max=Statistics.max(vals);
			double sum=Statistics.sum(vals);
			if(sum>minVal) {
				System.out.print(gene+"\t"+junction);
				for(int i=0; i<vals.length; i++) {
					double scaled=vals[i]/max;
					System.out.print("\t"+scaled);
				}
				System.out.println();
			}
		}
		reader.close();
		
	}
	
	
	private static double[] getVals2(String[] tokens) {
		double[] rtrn=new double[tokens.length-2];
		
		for(int i=2; i<tokens.length; i++) {
			rtrn[i-2]=Double.parseDouble(tokens[i]);
		}
		
		return rtrn;
	}

	
	private static void splitFileIntoChunks(String input, int chunkSize, String saveDir, String saveResults) throws IOException, InterruptedException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		
		int counter=0;
		int currentChunk=0;
		String nextLine;
		FileWriter writer=new FileWriter(saveDir+"/chunk"+currentChunk+".csv");
		while ((nextLine = reader.readLine()) != null) {
			if(counter>0) {
				int chunk=(counter-1)/chunkSize;
				if(chunk>currentChunk) {
					writer.close();
					writer=new FileWriter(saveDir+"/chunk"+chunk+".csv");
					currentChunk=chunk;
				}
				writer.write(nextLine+"\n");
			}
			counter++;
		}
		writer.close();
		reader.close();
		
		
		File[] files=new File(saveDir).listFiles();
		for(int i=0; i<files.length; i++) {
			String cmd="java -Xmx16000m -jar /central/groups/guttman/mguttman/Splicing/FitODE.jar "+files[i].getAbsolutePath()+" "+saveResults+"/"+files[i].getName();
			String name=files[i].getName();
			cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
			Process p=Runtime.getRuntime().exec(cmd);
			System.err.println(cmd);
		}
		
	}
	
	private static void merge(File[] listFiles) throws IOException {
		for(int i=0; i<listFiles.length; i++) {
			
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(listFiles[i])));
			String nextLine;
			while ((nextLine = reader.readLine()) != null) {
				System.out.println(nextLine);
			}
			reader.close();
		}
		
	}
	
	
	private static void fit4PL(String input) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		//FileWriter writer=new FileWriter(save);
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			if(counter>0) {
				
				double[] y=parseY(nextLine, x);
				double[] fit=fit4PL(x,y);
				System.out.println(nextLine+"\t"+fit[0]+"\t"+fit[1]+"\t"+fit[2]+"\t"+fit[3]+"\t"+fit[4]);
				
			}
			counter++;
		}
		reader.close();
		//writer.close();
	}
	
	private static double[] parseY(String nextLine, double[] x) {
		String[] tokens=nextLine.split("\t");
		double[] rtrn=new double[x.length];
		
		for(int i=9; i<tokens.length; i++) {
			rtrn[i-9]=Double.parseDouble(tokens[i]);
		}
		
		return rtrn;
	}


	private static void addDistanceToTSS(String input, GTFToJunctions gtfToJunctions) throws IOException {
		//Collection<Gene>junctions=gtfToJunctions.getAllJunctions();
		Map<SingleInterval, Collection<Integer>> junctionDistance=gtfToJunctions.getJunctionDistance();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		//int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			//if(counter>0) {
				String[] tokens=nextLine.split("\t");
				SingleInterval junction=new SingleInterval(tokens[1]);
				junction.setOrientation(Strand.POSITIVE);
				Collection<Integer> distances=junctionDistance.get(junction);
				
				if(distances==null) {
					junction.setOrientation(Strand.NEGATIVE);
					distances=junctionDistance.get(junction);
				}
				double min=Statistics.min(distances);
				double max=Statistics.max(distances);
				
				System.out.println(min+"\t"+max+"\t"+nextLine);
			//}
			//counter++;
		}
		reader.close();
		
		
		
		
	}

	
	

	public static void main(String[] args) throws IOException, NumberFormatException, InterruptedException {
		//match(args[0], args[1]);
		//filterJunctions(args[0], args[1], 0.75, 1);
		//writeJunctions(args[0], args[1], 0.75);
		//computeUniqueJunctions(args[0], args[1]);
		//getJunctions(args[0], args[1]);
		//filterRPKM(args[0], Double.parseDouble(args[1]));
		//writeBEDgraph(args[0]);
		//slidingWindow(args[0], args[1]);
		//averagePerGene(args[0], args[1]);
		//countBams(new File(args[0]).listFiles());
		
		//countFastqs(new File(args[0]).listFiles());
		
		//downsampleFastq(new File(args[0]), new File(args[1]), Double.parseDouble(args[2]), args[3]);
		
		
		//downsample(new File(args[0]).listFiles(), args[1], args[2]);
		//spliceCounts(new File(args[0]).listFiles(), args[1]);
		//downsampleCounts(args[0], weightFactors(args[1]), args[2], new GTFToJunctions(new File(args[3])));
		
		//computeSlopes(args[0], args[1]);
		
		//sumByGene(args[0], weightFactors(args[1]), new GTFToJunctions(new File(args[2])));
		
		//intronCounts(new File(args[0]).listFiles(), new GTFToJunctions(new File(args[1])));
		
		//scaleNorm(args[0], Double.parseDouble(args[1]));
		
		//TODO Local per gene downsample
		//TODO compute slope
		//compute4PL(args[0]);
		
		//localMaximum(args[0]);
	
		//ratioAtMax(args[0], args[1]);
		
		
		//splitFileIntoChunks(args[0], Integer.parseInt(args[1]), args[2], args[3]);
		
		//merge(new File(args[0]).listFiles());
		
		//addDistanceToTSS(args[0], new GTFToJunctions(new File(args[1])));
		
		
		fit4PL(args[0]);
	}



	

	


	



	
	


	

	

	


	


	

	


	


	
	
}

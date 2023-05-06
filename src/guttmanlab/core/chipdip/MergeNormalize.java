package guttmanlab.core.chipdip;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.datastructures.Pair;



public class MergeNormalize {

	static String sizes="/groups/guttman/mguttman/scripts/mm10.chrom.sizes";
	
	
	private static Map<String, Collection<File>> getFilesByName(File[] files) {
		Map<String, Collection<File>> rtrn=new TreeMap<String, Collection<File>>();
		
		for(File file: files) {
			String name=getName(file.getName()); //strip chr
			System.err.println(file.getName()+" "+name);
			if(!rtrn.containsKey(name)) {rtrn.put(name, new ArrayList<File>());}
			Collection<File> list=rtrn.get(name);
			list.add(file);
			rtrn.put(name, list);
		}
		
		return rtrn;
	}
	
	
	private static String getName(String name) {
		String rtrn="";
		String[] tokens=name.split("\\.");
		for(int i=1; i<tokens.length; i++) {
			rtrn+=tokens[i]+".";
		}
		rtrn=rtrn.substring(0, rtrn.length()-1);
		return rtrn;
	}
	
	
	private static void writeMerge(Map<String, Collection<File>> filesByName, String saveDir) throws IOException, InterruptedException {
		for(String name: filesByName.keySet()) {
			System.err.println("writing "+name);
			String fileName=saveDir+"/"+name;
			FileWriter writer=new FileWriter(fileName);
			
			//double normFactor=Double.MAX_VALUE;
			
			double sum=0;
			double count=0;
			
			for(File file: filesByName.get(name)) {
				Pair<Double> f=write(file, writer);
				sum+=f.getValue1();
				count+=f.getValue2();
			}
			
			double normFactor=sum/count;
			
			writer.close();
			
			System.err.println("normalizing "+name+" "+normFactor);
			fileName=normalize(fileName, normFactor);
			
			System.err.println("converting "+fileName+" to bigwig");
			convertToBigWig(fileName);
		}
		
	}


	private static String normalize(String fileName, double normFactor) throws IOException {
		String newFileName=fileName.replace(".bedgraph", ".norm.bedgraph");
		
		FileWriter writer=new FileWriter(newFileName);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			double val=Double.parseDouble(tokens[3]);
			double norm=(val/normFactor);
			writer.write(tokens[0]+"\t"+tokens[1]+"\t"+tokens[2]+"\t"+norm+"\n");
		}
		reader.close();
		writer.close();
		
		return newFileName;
	}


	private static void convertToBigWig(String input) throws IOException, InterruptedException {
		String output=input.replace(".bedgraph", ".bw");
		String cmd="/groups/guttman/software/kentUtils/bin/bedGraphToBigWig "+input+" "+sizes+" "+output;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
	}


	


	private static Pair<Double> write(File file, FileWriter writer) throws IOException {
		//double max=0;
		double sum=0;
		double count=0;
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			writer.write(nextLine+"\n");
			sum+=Double.parseDouble(nextLine.split("\t")[3]);
			count++;
		}
		reader.close();
		
		Pair<Double> rtrn=new Pair<Double>();
		rtrn.setValue1(sum);
		rtrn.setValue2(count);
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException, InterruptedException {
		File[] files=new File(args[0]).listFiles();
		Map<String, Collection<File>> filesByName=getFilesByName(files);
		String saveDir=args[1];
		
		writeMerge(filesByName, saveDir);
	}
	
}

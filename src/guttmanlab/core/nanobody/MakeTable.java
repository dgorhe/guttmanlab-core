package guttmanlab.core.nanobody;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.inference.ChiSquareTest;

import guttmanlab.core.annotation.io.BEDFileIO;

public class MakeTable {

	public MakeTable(String[] files, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		List<String> fullLines=BEDFileIO.loadLines(files[0]);
		List<Integer>[] scores=parse(files);
		
		long[] totals=total(scores);
		
		for(int position=0; position<fullLines.size(); position++) {
			String line=fullLines.get(position);
			double[] score=getScores(scores, position);
			writer.write(line.split("\t")[0]);
			for(int i=0; i<score.length; i++) {writer.write("\t"+score[i]);}
			writer.write("\n");
		}
		
		writer.close();
	}

	private long[] total(List<Integer>[] scores) {
		long[] totals=new long[scores.length];
		
		for(int i=0; i<scores.length; i++) {
			totals[i]=sum(scores[i]);
		}
		return totals;
	}

	private long sum(List<Integer> list) {
		long sum=0;
		
		for(Integer val: list) {sum+=val;}
		
		return sum;
	}

	private List<Integer>[] parse(String[] files) throws NumberFormatException, IOException {
		List<Integer>[] rtrn=new List[files.length];
		
		for(int i=0; i<files.length; i++) {
			rtrn[i]=parse(files[i]);
		}
		
		return rtrn;
	}

	private List<Integer> parse(String file) throws NumberFormatException, IOException {
		List<Integer> rtrn=new ArrayList<Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String line;
		int counter=0;
		while ((line = reader.readLine()) != null) {
			String[] tokens=line.split("\t");
			rtrn.add(Integer.parseInt(tokens[tokens.length-1]));
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		reader.close();
		return rtrn;
	}

	private double[] getScores(List<Integer>[] scores, int position) {
		double[] rtrn=new double[scores.length];
		
		for(int i=0; i<scores.length; i++) {
			rtrn[i]=scores[i].get(position);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		String[] files=new String[args.length-1];
		for(int i=0; i<files.length; i++) {
			files[i]=args[i];
		}
		String save=args[args.length-1];
		new MakeTable(files, save);
	}
	
}

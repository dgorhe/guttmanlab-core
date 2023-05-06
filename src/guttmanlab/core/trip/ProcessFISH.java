package guttmanlab.core.trip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.stat.StatUtils;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class ProcessFISH {
	
	double quantile=0.9;
	int numberOfPerms=100;

	public ProcessFISH(File input, String saveDir) throws IOException {
		FileWriter writer=new FileWriter(saveDir);
		Map<String, Collection<String>> linesByCell=parseByCell(input);
		
		Map<String, Integer> cellToNumber=cellToNumbers(linesByCell);
		
		for(String cell: linesByCell.keySet()) {
			Collection<String> lines=linesByCell.get(cell);
			scores(lines, writer, cell, cellToNumber);
			
		}
		
		writer.close();
	}
	
	
	private Map<String, Integer> cellToNumbers(Map<String, Collection<String>> linesByCell) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		Collection<String> nameSet=new TreeSet<String>();
		for(String cell: linesByCell.keySet()) {
			String condition=cell.split("\\.")[0];
			nameSet.add(condition);
		}
		
		int index=0;
		for(String cell: nameSet) {
			rtrn.put(cell, index);
			index++;
		}
		
		return rtrn;
	}


	private void scores(Collection<String> lines, FileWriter writer, String cell, Map<String, Integer> cellToNumber) throws IOException {
		
		List<Double> spot=new ArrayList<Double>();
		List<Double> notSpot=new ArrayList<Double>();
		
		int cellNumber=cellToNumber.get(cell.split("\\.")[0]);
		
		for(String line: lines) {
			String[] tokens=line.split(",");
			int isSpot=Integer.parseInt(tokens[8]);
			//String x=tokens[7];
			//String y=tokens[6];
			//double FISH=Double.parseDouble(tokens[12]);
			double mcherry=Double.parseDouble(tokens[11]);
			//writer.write(FISH+"\t"+mcherry+"\t"+isSpot+"\n");
			//rtrn.set(x, y, FISH);
			if(isSpot==1) {spot.add(mcherry);}
			if(isSpot==0) {notSpot.add(mcherry);}
		}
		
		for(Double score: spot) {
			double percentile2=Statistics.percentLessThan(score, notSpot);
			writer.write(cell+"\t"+cellNumber+"\t"+score+"\t"+percentile2+"\n");
		}
		
		/*double[] randomVals=sample(notSpot, spot.size());
		double score=score(spot);
		double percentile=Statistics.percentLessThan(score, randomVals);
		double percentile2=Statistics.percentLessThan(score, notSpot);
		double zscore=Statistics.zScore(score, randomVals);		
		writer.write(cell+"\t"+score+"\t"+percentile+"\t"+percentile2+"\t"+zscore);
		for(int i=0; i<randomVals.length; i++) {
			writer.write("\t"+randomVals[i]);
		}
		writer.write("\n");*/
	}
	
	


	private Map<String, Collection<String>> parseByCell(File input) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		List<String> lines=BEDFileIO.loadLines(input.getAbsolutePath(), 1);
		
		for(String line: lines) {
			line=line.replaceAll("\"", "");
			String[] tokens=line.split(",");
			String cell=tokens[0]+"."+tokens[1]+"."+tokens[2]+"."+tokens[3];
			if(!rtrn.containsKey(cell)) {rtrn.put(cell, new ArrayList<String>());}
			
			Collection<String> list=rtrn.get(cell);
			list.add(line);
		}
		
		return rtrn;
	}
	
	
	private static Map<String, MatrixWithHeaders> matrixByCell(File input, String name, int type) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		List<String> lines=BEDFileIO.loadLines(input.getAbsolutePath(), 1);
		
		for(String line: lines) {
			line=line.replaceAll("\"", "");
			String[] tokens=line.split(",");
			String cell=getCell(tokens, type, name);
			if(!rtrn.containsKey(cell)) {rtrn.put(cell, new ArrayList<String>());}
			
			Collection<String> list=rtrn.get(cell);
			list.add(line);
		}
		
		
		Map<String, MatrixWithHeaders> map=new TreeMap<String, MatrixWithHeaders>();
		for(String cell: rtrn.keySet()) {
			System.err.println(cell);
			Collection<String> list=rtrn.get(cell);
			MatrixWithHeaders mwh=makeMatrix(list, cell, type);
			map.put(cell, mwh);
		}
		
		
		return map;
	}

	
	private static String getCell(String[] tokens, int type, String name) {
		if(type==0) {return tokens[0]+"."+tokens[1]+"."+tokens[2]+"."+tokens[3];}
		return name+"."+tokens[0]+"."+tokens[1]+"."+tokens[2];
	}


	private static void writeDispersion(String name, List<Double> list) {
		double[] vals=makeArray(list);
		double[] mode=StatUtils.mode(vals);
		double percentile99=Statistics.quantile(vals, 0.99);
		double percentile95=Statistics.quantile(vals, 0.95);
		double percentile90=Statistics.quantile(vals, 0.90);
		double total=Statistics.sum(vals);
		double mean=Statistics.mean(vals);
		double median=Statistics.quantile(vals, 0.5);
		double max=Statistics.max(vals);
		double ratio=percentile99/mode[0];
		System.out.println(name+"\t"+total+"\t"+mean+"\t"+median+"\t"+max+"\t"+mode[0]+"\t"+percentile99+"\t"+percentile90+"\t"+percentile95+"\t"+ratio+"\t"+mode.length);
	}
	
	private static double[] makeArray(List<Double> list) {
		double[] rtrn=new double[list.size()];
		for(int i=0; i<list.size(); i++) {rtrn[i]=list.get(i);}
		return rtrn;
	}


	private static MatrixWithHeaders makeMatrix(Collection<String> list, String name, int parseType) {
		Collection<String> rows=new TreeSet<String>();
		Collection<String> columns=new TreeSet<String>();
		
		
		for(String line: list) {
			String[] tokens=line.split(",");
			String x=getX(tokens, parseType);
			String y=getY(tokens, parseType);
			rows.add(x);
			columns.add(y);
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(new ArrayList<String>(rows), new ArrayList<String>(columns));
		
		List<Double> vals=new ArrayList<Double>();
		
		for(String line: list) {
			String[] tokens=line.split(",");
			String x=getX(tokens, parseType);
			String y=getY(tokens, parseType);
			double val=getmCherry(tokens, parseType);
			rtrn.set(x, y, val);
			vals.add(val);
		}
		
		writeDispersion(name, vals);
		//scoreSquares(name, rtrn);
		
		return rtrn;
	}


	private static void scoreSquares(String name, MatrixWithHeaders mwh) {
		int size=3;
		List<Double> vals=new ArrayList<Double>();
		for(int i=0; i<mwh.getNumberRows(); i++) {
			for(int j=0; j<mwh.getNumberColumns(); j++) {
				double sum=getValuesInSquare(mwh, i,j, size);
				if(sum>0) {
					vals.add(sum);
				}
			}
		}
		
		writeDispersion(name, vals);
		
		
		
		
	}


	private static double getValuesInSquare(MatrixWithHeaders mwh, int i, int j, int size) {
		double sum=0;
		for(int row=i-size; row<=i+size; row++) {
			for(int column=j-size; column<=j+size; column++) {
				if(row>0 && column>0 && row<mwh.getNumberRows() && column<mwh.getNumberColumns()) {
					double score=mwh.get(row, column);
					sum+=score;
					if(score==0) {return 0;}
				}
				else {return 0;}
			}
		}
		return sum;
	}


	private static String getX(String[] tokens, int parseType) {
		if(parseType==0) {return tokens[7];}
		return tokens[5];
	}
	
	private static String getY(String[] tokens, int parseType) {
		if(parseType==0) {return tokens[6];}
		return tokens[4];
	}
	
	private static double getmCherry(String[] tokens, int parseType) {
		if(parseType==0) {return Double.parseDouble(tokens[11]);}
		return Double.parseDouble(tokens[8]);
	}


	public static void main(String[] args) throws IOException {
		File input=new File(args[0]);
		String saveDir=args[1];
		int parseType=0;
		String name="";
		if(args.length>2) {
			name=args[2];
			parseType=1;
		}
		
		//double quantile=Double.parseDouble(args[2]);
		Map<String, MatrixWithHeaders> map=matrixByCell(input, name, parseType);
		
		for(String cell: map.keySet()) {
			String save=saveDir+"/"+cell+".matrix";
			MatrixWithHeaders mwh=map.get(cell);
			mwh.write(save);
		}
		
	}
	
}

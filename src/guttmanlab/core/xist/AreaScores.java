package guttmanlab.core.xist;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.stat.StatUtils;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class AreaScores {

	public static double getRelativeIntensityScores(File file) throws IOException {
		Map<String, List<Double>> valsByName=parse(file);
		
		System.out.println("Cell\tTotal\tAvg\tMedian\tMax\tMode\t99%\t90%\t95%\tRatio\t# of modes");
		
		for(String name: valsByName.keySet()) {
			System.err.println(name+" "+valsByName.get(name).size());
			double[] vals=makeArray(valsByName.get(name));
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
		
		return 0.0;
	}

	
	private static double[] makeArray(List<Double> list) {
		double[] rtrn=new double[list.size()];
		
		for(int i=0; i<list.size(); i++) {rtrn[i]=list.get(i);}
		
		return rtrn;
	}


	/*private static double[] parse(File file) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file.getAbsolutePath(), 1);
		double[] rtrn=new double[lines.size()];
		
		for(int i=0; i<lines.size(); i++) {
			rtrn[i]=Double.parseDouble(lines.get(i));
		}
		
		return rtrn;
	}*/
	
	private static Map<String, List<Double>> parse(File file) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file.getAbsolutePath());
		
		Map<String, List<Double>> rtrn=new TreeMap<String, List<Double>>();
		
		String[] names=lines.get(0).split(",");
		
		for(int i=0; i<names.length; i++) {
			rtrn.put(names[i], new ArrayList<Double>());
		}
		
		for(int i=1; i<lines.size(); i++) {
			String[] vals=lines.get(i).split(",");
			for(int j=0; j<vals.length; j++) {
				String name=names[j];
				if(!vals[j].isEmpty()) {
					double val=Double.parseDouble(vals[j]);
					rtrn.get(name).add(val);
				}
			}
		}
		
		return rtrn;
	}


	private static double[] getVals(MatrixWithHeaders mwh) {
		double[] rtrn=new double[mwh.getNumberRows()*mwh.getNumberColumns()];
		
		int counter=0;
		for(int i=0; i<mwh.getNumberRows(); i++) {
			for(int j=0; j<mwh.getNumberColumns(); j++) {
				rtrn[counter]=mwh.get(i, j);
				counter++;
			}
		}
		return rtrn;
	}


	public static double getAreaScores(File file, int size) throws IOException {
		MatrixWithHeaders mwh=XYToMatrix(file);
		
		double total=0;
		List<Double> vals=new ArrayList<Double>();
		for(int i=0; i<mwh.getNumberRows(); i++) {
			for(int j=0; j<mwh.getNumberColumns(); j++) {
				total+=mwh.get(i, j);
				double sum=getValuesInSquare(mwh, i,j, size);
				if(sum>0) {
					vals.add(sum);
					//writer.write(sum+"\n");
				}
				//square.set(i, j, sum);
			}
		}
		//writer.close();
		double ratio01=Statistics.quantile(vals, 0.999)/Statistics.quantile(vals, 0.001);
		double ratio1=Statistics.quantile(vals, 0.99)/Statistics.quantile(vals, 0.01);
		double ratio5=Statistics.quantile(vals, 0.95)/Statistics.quantile(vals, 0.05);
		double ratio10=Statistics.quantile(vals, 0.90)/Statistics.quantile(vals, 0.1);
		System.out.println(total+"\t"+ratio01+"\t"+ratio1+"\t"+ratio5+"\t"+ratio10);
		return ratio10;
		//square.write(save);
	}
	
	
	public static MatrixWithHeaders XYToMatrix(File file) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file.getAbsolutePath(),1);
		
		List<String> rows=getRows(lines);
		List<String> columns=getColumns(lines);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			if(tokens.length>5) {
				mwh.set(tokens[0], tokens[1], Double.parseDouble(tokens[5]));
			}
		}
		
		return mwh;
	}
	
	
	private static List<String> getColumns(List<String> lines) {
		Set<String> rtrn=new TreeSet<String>();
			
		for(String line: lines) {
			if(line.split("\t").length>2) {
				rtrn.add(line.split("\t")[1]);
			}
		}
			
		ArrayList<String> list=new ArrayList<String>();
		list.addAll(rtrn);
			
		return list;
	}

	private static List<String> getRows(List<String> lines) {
		Set<String> rtrn=new TreeSet<String>();
		
		for(String line: lines) {
			if(line.split("\t").length>2) {
				rtrn.add(line.split("\t")[0]);
			}
		}
		
		ArrayList<String> list=new ArrayList<String>();
		list.addAll(rtrn);
		
		return list;
	}

	private static double getValuesInSquare(MatrixWithHeaders mwh, int i, int j, int size) {
		double sum=0;
		for(int row=i-size; row<=i+size; row++) {
			for(int column=j-size; column<=j+size; column++) {
				//System.err.println(i+" "+j+" "+row+" "+column);
				if(row>0 && column>0 && row<mwh.getNumberRows() && column<mwh.getNumberColumns()) {
					sum+=mwh.get(row, column);
				}
				else {return 0;}
			}
		}
		
		return sum;
	}
	
	
	public static void main(String[] args) throws IOException {
		/*File[] files=new File(args[0]).listFiles();
		int size=Integer.parseInt(args[1]);
		
		for(int i=0; i<files.length; i++) {
			MatrixWithHeaders mwh=new MatrixWithHeaders(files[i]);
			double ratio=getAreaScores(mwh, size);
			//System.out.println(files[i].getName()+"\t"+ratio);
		}*/
		
		File file=new File(args[0]);
		getRelativeIntensityScores(file);
		
	}
	
}

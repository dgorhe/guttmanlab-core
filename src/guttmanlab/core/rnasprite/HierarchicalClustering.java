package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.ml.distance.EuclideanDistance;

import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;

public class HierarchicalClustering {

	public HierarchicalClustering(MatrixWithHeaders mwh, String save) throws IOException{
		List<String> order=cluster(mwh);
		write(save, order);
	}
	
	
	public static List<String> clusterColumns(MatrixWithHeaders mwh) {
		List<String> rtrn=new ArrayList<String>();
		
		MatrixWithHeaders distanceMatrix=correlationMatrix(mwh, false);
		
		while(distanceMatrix.getRowNames().size()>1){
			int size=distanceMatrix.getRowNames().size();
			System.err.println(size);
			//Get minimum and merge
			Pair<String> pair=getMinimum(distanceMatrix);
			//System.err.println(pair.getValue1()+" "+pair.getValue2());
			add(pair, rtrn);
			//connect nodes
			distanceMatrix=updateMatrix(distanceMatrix, pair);
		}
		
		
		
		return rtrn;
	}
	
	
	public static List<String> clusterRows(MatrixWithHeaders mwh) {
		List<String> rtrn=new ArrayList<String>();
		
		MatrixWithHeaders distanceMatrix=correlationMatrix(mwh, true);
		
		while(distanceMatrix.getRowNames().size()>1){
			int size=distanceMatrix.getRowNames().size();
			System.err.println(size);
			//Get minimum and merge
			Pair<String> pair=getMinimum(distanceMatrix);
			//System.err.println(pair.getValue1()+" "+pair.getValue2());
			add(pair, rtrn);
			//connect nodes
			distanceMatrix=updateMatrix(distanceMatrix, pair);
		}
		
		
		
		return rtrn;
	}
	
	public static List<String> cluster(MatrixWithHeaders mwh) {
		List<String> rtrn=new ArrayList<String>();
		
		MatrixWithHeaders correlationMatrix=correlationMatrix(mwh, true);
		
		while(correlationMatrix.getRowNames().size()>1){
			System.err.println(correlationMatrix.getRowNames().size());
			//Get minimum and merge
			Pair<String> pair=getMaximum(correlationMatrix);
			System.err.println(pair.getValue1()+" "+pair.getValue2());
			add(pair, rtrn);
			//connect nodes
			correlationMatrix=updateMatrix(correlationMatrix, pair);
		}
		
		
		
		return rtrn;
	}
	
	private static Pair<String> getMaximum(MatrixWithHeaders correlationMatrix) {
		double maxVal=-Double.MAX_VALUE;
		Pair<String> maxPair=new Pair<String>();
		for(String gene1: correlationMatrix.getRowNames()){
			for(String gene2: correlationMatrix.getColumnNames()){
				if(!gene1.equals(gene2)){
					double val=correlationMatrix.get(gene1, gene2);
					if(val>maxVal){
						maxVal=val;
						maxPair.setValue1(gene1);
						maxPair.setValue2(gene2);
					}
				}
			}
		}
		
		return maxPair;
	}
	
	private static Pair<String> getMinimum(MatrixWithHeaders correlationMatrix) {
		double minVal=Double.MAX_VALUE;
		Pair<String> minPair=new Pair<String>();
		for(String gene1: correlationMatrix.getRowNames()){
			for(String gene2: correlationMatrix.getColumnNames()){
				if(!gene1.equals(gene2)){
					double val=correlationMatrix.get(gene1, gene2);
					if(val<minVal){
						minVal=val;
						minPair.setValue1(gene1);
						minPair.setValue2(gene2);
					}
				}
			}
		}
		
		return minPair;
	}

	private static MatrixWithHeaders updateMatrix(MatrixWithHeaders correlationMatrix, Pair<String> pair) {
		correlationMatrix=updateRows(correlationMatrix, pair);
		correlationMatrix=updateColumns(correlationMatrix, pair);
		return correlationMatrix;	
	}


	private static MatrixWithHeaders updateColumns(MatrixWithHeaders correlationMatrix, Pair<String> pair) {
		List<String> updatedColumns=correlationMatrix.getColumnNames();
		
		String newColumn=pair.getValue1()+"-"+pair.getValue2();
		updatedColumns.add(newColumn);
		updatedColumns.remove(pair.getValue1());
		updatedColumns.remove(pair.getValue2());
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(correlationMatrix.getRowNames(), updatedColumns);
		
		//merge values
		double[] vals1=correlationMatrix.getColumn(pair.getValue1());
		double[] vals2=correlationMatrix.getColumn(pair.getValue2());
		double[] mergedVals=merge(vals1, vals2);
		
		for(String column: updatedColumns){
			if(column.equals(newColumn)){
				rtrn.setColumn(mergedVals, newColumn);
			}
			
			else{
				rtrn.setColumn(correlationMatrix.getColumn(column), column);
			}
		}
		
		return rtrn;
	}


	private static MatrixWithHeaders updateRows(MatrixWithHeaders correlationMatrix, Pair<String> pair) {
		List<String> updatedRows=correlationMatrix.getRowNames();
		
		String newRow=pair.getValue1()+"-"+pair.getValue2();
		updatedRows.add(newRow);
		updatedRows.remove(pair.getValue1());
		updatedRows.remove(pair.getValue2());
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(updatedRows, correlationMatrix.getColumnNames());
		
		//merge values
		double[] vals1=correlationMatrix.getRow(pair.getValue1());
		double[] vals2=correlationMatrix.getRow(pair.getValue2());
		double[] mergedVals=merge(vals1, vals2);
		
		for(String row: updatedRows){
			if(row.equals(newRow)){
				rtrn.setRow(newRow, mergedVals);
			}
			
			else{
				rtrn.setRow(row, correlationMatrix.getRow(row));
			}
		}
		
		return rtrn;
	}


	private static double[] merge(double[] vals1, double[] vals2) {
		double[] rtrn=new double[vals1.length];
		
		for(int i=0; i<vals1.length; i++){
			//rtrn[i]=(vals1[i]+vals2[i]/2.0);
			rtrn[i]=Math.max(vals1[i], vals2[i]);
		}
		
		return rtrn;
	}


	private static void add(Pair<String> pair, List<String> rtrn) {
		String string1=pair.getValue1();
		String string2=pair.getValue2();
		if(!string1.contains("-")){rtrn.add(string1);}
		if(!string2.contains("-")){rtrn.add(string2);}
	}


	private static MatrixWithHeaders correlationMatrix(MatrixWithHeaders mwh, boolean rows) {
		//PearsonsCorrelation p=new PearsonsCorrelation();
		EuclideanDistance d=new EuclideanDistance();
		
		List<String> list=mwh.getRowNames();
		if(!rows){list=mwh.getColumnNames();}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(list, list);
		int counter=0;
		for(String gene1: list){
			//System.err.println(gene1);
			double[] vals1;
			if(rows){vals1=mwh.getRow(gene1);}
			else{vals1=mwh.getColumn(gene1);}
			for(String gene2: list){
				double[] vals2;
				if(rows){vals2=mwh.getRow(gene2);}
				else{vals2=mwh.getColumn(gene2);}
				//double correlation=p.correlation(vals1, vals2);
				double distance=d.compute(vals1, vals2);
				//if(correlation>0.7){System.out.println(gene1+"\t"+gene2+"\t"+correlation);}
				rtrn.set(gene1, gene2, distance);
			}
			counter++;
			if(counter%100 ==0){System.err.println("Building correlation matrix ... "+counter+" "+mwh.getRowNames().size());}
			//System.err.println(gene1+" "+counter+" "+mwh.getRowNames().size());
		}
		return rtrn;
	}


	private void write(String save, List<String> order) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("ID\tOrder\n");
		
		int index=0;
		for(String name: order){
			writer.write(name+"\t"+index+"\n");
			index++;
		}
		
		writer.close();
	}


	public static void main(String[] args) throws IOException{
		MatrixWithHeaders mwh=new MatrixWithHeaders(new File(args[0]));
		String save=args[1];
		new HierarchicalClustering(mwh, save);
	}
	
}

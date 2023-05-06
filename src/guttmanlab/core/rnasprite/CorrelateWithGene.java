package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.stat.correlation.PearsonsCorrelation;

import guttmanlab.core.datastructures.MatrixWithHeaders;

public class CorrelateWithGene {
	
	int minClusterCount=10;

	public CorrelateWithGene(MatrixWithHeaders mwh, String save) throws IOException{
		/*mwh=filterByInputCount(mwh);
		mwh=removeInput(mwh);*/
		
		MatrixWithHeaders correlationMatrix=correlationMatrix(mwh);
		
		correlationMatrix.write(save);
		
		//MatrixWithHeaders
		
		/*FileWriter writer=new FileWriter(save);
		PearsonsCorrelation p=new PearsonsCorrelation();
		SpearmansCorrelation s=new SpearmansCorrelation();
		double[] referencePattern=mwh.getRow(gene);
		for(String gene1: mwh.getRowNames()){
			double[] comparePattern=mwh.getRow(gene1);
			double correlation=p.correlation(referencePattern, comparePattern);
			double spearman=s.correlation(referencePattern, comparePattern);
			System.err.println(gene1+" "+correlation+" "+spearman);
			writer.write(gene1+"\t"+correlation+"\t"+spearman+"\n");
		}
		writer.close();*/
	}
	
	private MatrixWithHeaders correlationMatrix(MatrixWithHeaders mwh) {
		PearsonsCorrelation p=new PearsonsCorrelation();
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getRowNames());
		int counter=0;
		for(String gene1: mwh.getRowNames()){
			//System.err.println(gene1);
			double[] vals1=mwh.getRow(gene1);
			for(String gene2: mwh.getRowNames()){
				double[] vals2=mwh.getRow(gene2);
				double correlation=p.correlation(vals1, vals2);
				if(correlation>0.7){System.out.println(gene1+"\t"+gene2+"\t"+correlation);}
				rtrn.set(gene1, gene2, correlation);
			}
			counter++;
			System.err.println(gene1+" "+counter+" "+mwh.getRowNames().size());
		}
		return rtrn;
	}

	private MatrixWithHeaders filterByInputCount(MatrixWithHeaders mwh) {
		List<String> geneNames=new ArrayList<String>();
		for(String gene: mwh.getRowNames()){
			double inputCount=mwh.get(gene, "Input");
			if(inputCount>this.minClusterCount){
				geneNames.add(gene);
			}
		}
		
		return mwh.submatrixByRowNames(geneNames);
		
	}

	private MatrixWithHeaders removeInput(MatrixWithHeaders mwh) {
		List<String> rows=new ArrayList<String>();
		rows.addAll(mwh.getRowNames());
		for(String row: mwh.getRowNames()){
			if(row.contains("Input")){rows.remove(row);}
		}
		
		
		List<String> columns=new ArrayList<String>();
		columns.addAll(mwh.getColumnNames());
		for(String column: mwh.getColumnNames()){
			if(column.contains("Input")){columns.remove(column);}
			if(column.contains("chr17:39000000")){	
				System.err.println(column);
				columns.remove(column);
			}
		}
		
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(String row: rows){
			for(String column: columns){
				rtrn.set(row, column, mwh.get(row, column));
			}
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		MatrixWithHeaders mwh=new MatrixWithHeaders(new File(args[0]));
		//MatrixWithHeaders transpose=mwh.transpose();
		String save=args[1];
		new CorrelateWithGene(mwh, save);
	}
	
}

package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class XChromosomeSingleCell {

	double minMedian=2.0;
	
	public XChromosomeSingleCell(File dataFile, File featuresToCoordinates, String save, String geneName) throws IOException{
		Map<String, String> geneToChr=parse(featuresToCoordinates);
		
		MatrixWithHeaders data=new MatrixWithHeaders(dataFile, geneToChr);
		
		data=normalize(data);
		
		
		//data.writeGCT(save+".norm.gct");
		
		//System.err.println(data.rowDimension()+" "+filtered.rowDimension());
		
		//filtered.writeGCT(save+".filtered.gct");
		
		MatrixWithHeaders XGenes=pullXGenes(data, geneToChr);
		
		//XGenes.writeGCT(save+".x.gct");
		
		
		MatrixWithHeaders XistPositive=getXist(XGenes);
		MatrixWithHeaders XistNegative=getXistNeg(XGenes);
		MatrixWithHeaders XistZero=getXistZero(XGenes);
		
		XistZero.appendColumns(XistPositive);
		XistZero=filterGenes(XistZero);
		MatrixWithHeaders transpose=XistZero.transpose();
		
		transpose.writeGCT(save+".merged.gct");
		
		
		/*double[] pos=XistPositive.getRow(geneName);
		double[] neg=XistZero.getRow(geneName);
		
		FileWriter writer=new FileWriter(save);
		writer.write("Pos");
		for(int i=0; i<pos.length; i++){
			writer.write("\t"+pos[i]);
		}
		writer.write("\n");
		
		writer.write("Neg");
		for(int i=0; i<neg.length; i++){
			writer.write("\t"+neg[i]);
		}
		writer.write("\n");
		writer.close();*/
		
		/*XistNegative.writeGCT(save+".XistNeg.gct");
		XistPositive.writeGCT(save+".XistPos.gct");
		XistZero.writeGCT(save+".XistZero.gct");*/
		
	}
	
	
	private MatrixWithHeaders normalize(MatrixWithHeaders data) {
		Map<String, Double> normFactor=new TreeMap<String, Double>();
		for(String cell: data.getColumnNames()){
			double[] vals=data.getColumn(cell);
			double sum=1/Statistics.sum(vals)*1000000;
			normFactor.put(cell, sum);
			
		}
		
		return data.multiplyColumnsWithConstants(normFactor);
	}

	private MatrixWithHeaders getXistZero(MatrixWithHeaders xGenes) {
		List<String> list=new ArrayList<String>();
		for(String sample: xGenes.getColumnNames()){
			double val=xGenes.get("Xist", sample);
			if(val==0){list.add(sample);}
		}
		return xGenes.submatrixByColumnNames(list);
	}

	private MatrixWithHeaders getXistNeg(MatrixWithHeaders xGenes) {
		List<String> list=new ArrayList<String>();
		for(String sample: xGenes.getColumnNames()){
			double val=xGenes.get("Xist", sample);
			if(val>0 && val<5){list.add(sample);}
		}
		return xGenes.submatrixByColumnNames(list);
	}
	

	private MatrixWithHeaders getXist(MatrixWithHeaders xGenes) {
		List<String> list=new ArrayList<String>();
		for(String sample: xGenes.getColumnNames()){
			double val=xGenes.get("Xist", sample);
			if(val>1000){list.add(sample);}
		}
		return xGenes.submatrixByColumnNames(list);
	}



	private MatrixWithHeaders pullXGenes(MatrixWithHeaders filtered, Map<String, String> geneToChr) {
		List<String> list=new ArrayList<String>();
		for(String gene: filtered.getRowNames()){
			String chr=geneToChr.get(gene);
			System.err.println(gene+" "+chr);
			if(chr!=null && chr.equalsIgnoreCase("chrX")){list.add(gene);}
		}
		return filtered.submatrixByRowNames(list);
	}

	private Map<String, String> parse(File featuresToCoordinates) throws IOException {
		Map<String, String> geneToChr=new TreeMap<String, String>();
		List<String> lines=BEDFileIO.loadLines(featuresToCoordinates.getAbsolutePath());
		
		for(String line: lines){
			if(!line.startsWith("#")){
				String[] tokens=line.split("\t");
				geneToChr.put(tokens[3], tokens[0]);
			}
			
		}
		return geneToChr;
	}

	private MatrixWithHeaders filterGenes(MatrixWithHeaders data) {
		List<String> list=new ArrayList<String>();
		for(String gene: data.getRowNames()){
			double[] vals=data.getRow(gene);
			double percentZero=fractionZero(vals);
			if(percentZero<0.9){
				list.add(gene);
			}
		}
		return data.submatrixByRowNames(list);
	}

	private double fractionZero(double[] vals) {
		int countZero=0;
		for(int i=0; i<vals.length; i++){
			if(vals[i]==0){countZero++;}
		}
		return (double)countZero/(double)vals.length;
	}


	public static void main(String[] args) throws IOException{
		File dataFile=new File(args[0]);
		File features=new File(args[1]);
		String save=args[2];
		String geneName=args[3];
		new XChromosomeSingleCell(dataFile, features, save, geneName);
	}
	
}

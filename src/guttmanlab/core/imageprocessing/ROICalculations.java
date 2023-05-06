package guttmanlab.core.imageprocessing;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.math.Statistics;

public class ROICalculations {

	public ROICalculations(File[] files, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		writer.write("Name\tMean\tMedian\t95 Percentile\tMax\n");
		
		for(int i=0; i<files.length; i++){
			List<String> lines=BEDFileIO.loadLines(files[i].getAbsolutePath());
			double[] vals=getVals(lines);
			writer.write(files[i].getName()+"\t"+Statistics.mean(vals)+"\t"+Statistics.quantile(vals, 0.5)+"\t"+Statistics.quantile(vals, 0.95)+"\t"+Statistics.max(vals)+"\n");
		}
		
		writer.close();
	}

	private double[] getVals(List<String> lines) {
		double[] rtrn=new double[lines.size()-1];
		int counter=0;
		for(String line: lines){
			if(!line.split("\t")[1].equals("Value")){
				double val=new Double(line.split("\t")[1]).doubleValue();
				rtrn[counter]=val;
				counter++;
			}
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		File[] files=new File("/Users/mguttman/Downloads/Individual_ROI_Histograms_880_11-8-18/1K").listFiles();
		String save="/Users/mguttman/Downloads/Individual_ROI_Histograms_880_11-8-18/1K.txt";
		new ROICalculations(files, save);
	}
	
}

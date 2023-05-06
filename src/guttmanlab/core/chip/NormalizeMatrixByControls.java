package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class NormalizeMatrixByControls {
	double minVal=1.0;

	public NormalizeMatrixByControls(MatrixWithHeaders mwh, List<String> controls, String save) throws IOException {
		List<String> samples=pullSamples(mwh.getColumnNames(), controls);
		/*FileWriter writerMedian=new FileWriter(save+".median.bedgraph");
		FileWriter writerMean=new FileWriter(save+".mean.bedgraph");
		FileWriter writer95=new FileWriter(save+".95.bedgraph");
		FileWriter writer75=new FileWriter(save+".75.bedgraph");
		FileWriter writerMax=new FileWriter(save+".max.bedgraph");*/
		
		MatrixWithHeaders norm=new MatrixWithHeaders(mwh.getRowNames(), samples);
		
		
		for(String row: norm.getRowNames()) {
			double[] controlVals=getVals(mwh, row, controls);
			/*SingleInterval region=new SingleInterval(row);
			writerMedian.write(region.toBedgraph(Statistics.quantile(controlVals, 0.5))+"\n");
			writerMean.write(region.toBedgraph(Statistics.mean(controlVals))+"\n");
			writer95.write(region.toBedgraph(Statistics.quantile(controlVals, 0.95))+"\n");
			writer75.write(region.toBedgraph(Statistics.quantile(controlVals, 0.75))+"\n");
			writerMax.write(region.toBedgraph(Statistics.max(controlVals))+"\n");*/
			
			for(String col: norm.getColumnNames()) {
				double val=mwh.get(row, col);
				double normVal=val-Statistics.max(controlVals);
				norm.set(row, col, normVal);
			}
		}
		
		/*writerMedian.close();
		writerMean.close();
		writer95.close();
		writer75.close();
		writerMax.close();*/
		//writeBedgraphs(norm, save);
		//norm.write(save);
	}
	
	
	private void writeBedgraphs(MatrixWithHeaders norm, String saveDir) throws IOException {
		for(String col: norm.getColumnNames()) {
			System.err.println(col);
			write(saveDir+"/"+col+".bedgraph", norm, col);
		}
		
	}
	
	private void write(String string, MatrixWithHeaders data, String col) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(String row: data.getRowNames()) {
			SingleInterval region=new SingleInterval(row);
			double val=data.get(row, col);
			writer.write(region.toBedgraph(val)+"\n");
		}
		
		writer.close();
		
	}


	private List<String> pullSamples(List<String> columnNames, List<String> controls) {
		List<String> rtrn=new ArrayList<String>();
		
		for (String c: columnNames) {
			if(!controls.contains(c)) {rtrn.add(c);}
		}
		
		return rtrn;
	}


	private double[] getVals(MatrixWithHeaders mwh, String row, List<String> controls) {
		double[] rtrn=new double[controls.size()];
		
		int i=0;
		for(String c: controls) {
			rtrn[i]=mwh.get(row, c);
			i++;
		}
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException {
		MatrixWithHeaders data=new MatrixWithHeaders(new File(args[0]));
		List<String> controlList=BEDFileIO.loadLines(args[1]);
		String save=args[2];
		new NormalizeMatrixByControls(data, controlList, save);
	}
	
}

package guttmanlab.core.splicing.speckle.kineticmodel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

/**
 * 
 * @author mguttman
 * Take observed values and try and fit best params that fit the data using simulation model written by Ben Yeh 
 */

public class FitBestParams {

	static int numIter=2;
	static double elongationRate=50;
	static int timeShift=600;
	
	/**
	 * Given a set of params, generate data and compute distance to observed
	 * @param initiationRate rate of initiation
	 * @param splicingRate rate of splicing
	 * @param decayRate rate of decay (turnover)
	 * @param elongationRate rate of elongation
	 * @param annotation a file specifying the annotation structure
	 * @return MatrixWithHeaders expected count values for exons, introns, and junctions
	 * @throws IOException 
	 */
	private static MatrixWithHeaders simulate(double initiationRate, double splicingRate, double decayRate, double elongationRate, File annotation) throws IOException {
		//initiation, decay, splicing, elongation
		String cmd="python3 scripts/simulate_main.py --file_annot "+annotation+" --params "+initiationRate+" "+decayRate +" "+splicingRate +" "+elongationRate +" --n_time_steps 15000 -n 10 --t_wash "+timeShift;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		
		InputStream stream=p.getInputStream();
		BufferedReader reader=new BufferedReader(new InputStreamReader(stream));
		
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(reader, "\t", 0);
		rtrn=removeLastIntron(rtrn);
		
		addFractionSpliced(rtrn);
		
		//rtrn=columnNormalize(rtrn);
		
		return rtrn;
	}
	
	
	private static void addFractionSpliced(MatrixWithHeaders rtrn) {
		List<String> newColumns=new ArrayList<String>();
		for(String col: rtrn.getColumnNames()) {
			if(col.startsWith("spliced")) {
				String num=col.split("_")[1];
				rtrn.addColumn("ratio_"+num);
				newColumns.add("ratio_"+num);
			}
		}
		
		for(String row: rtrn.getRowNames()) {
			for(String column: newColumns) {
				String num=column.split("_")[1];
				double val1=rtrn.get(row, "spliced_"+num);
				double val2=rtrn.get(row, "exon_"+num);
				double ratio=val1/val2;
				System.err.println(row+" "+column+" "+val1 +" "+val2);
				rtrn.set(row, column, ratio);
			}
		}
		
	}


	private static MatrixWithHeaders columnNormalize(MatrixWithHeaders rtrn) {
		MatrixWithHeaders norm=new MatrixWithHeaders(rtrn.getRowNames(), rtrn.getColumnNames());
		for(String col: rtrn.getColumnNames()) {
			double max=Statistics.max(rtrn.getColumn(col));
			System.err.println(col+" "+max);
			for(String row: rtrn.getRowNames()) {
				double val=rtrn.get(row, col)/max;
				norm.set(row, col, val);
			}
		}
		return norm;
	}


	private static MatrixWithHeaders simulate(double[] paramFit, File annotation) throws IOException {
		return simulate(paramFit[0], paramFit[1], paramFit[2], elongationRate, annotation);
	}
	
	private static MatrixWithHeaders removeLastIntron(MatrixWithHeaders rtrn) {
		List<String> columns=new ArrayList<String>();
		columns.addAll(rtrn.getColumnNames());
		List<String> introns=new ArrayList<String>();
		
		for(String c: rtrn.getColumnNames()) {
			if(c.startsWith("intron")) {introns.add(c);}
		}
		
		columns.remove(introns.get(introns.size()-1));
		
		rtrn=rtrn.submatrixByColumnNames(columns);
		return rtrn;
	}

	/**
	 * 
	 * @param observed actual observed values for exons, introns, junctions per unit time
	 * @param expected simulated values for defined params
	 * @return double distance between these two matrices
	 */
	private static double distance(MatrixWithHeaders observed, MatrixWithHeaders expected) {
		//Start with simple least squares distance
		//take times from observed and pull from expected
		
		List<String> rows=observed.getRowNames();
		expected=expected.submatrixByRowNames(rows);
		
		
		Matrix data1=observed.getData();
		Matrix data2=expected.getData();
		
		double dist=0;
		double count=0;
		for(int i=0; i<data1.getRowDimension(); i++) {
			for(int j=0; j<data1.getColumnDimension(); j++) {
				double val1=data1.get(i, j);
				double val2=data2.get(i, j);
				double diff=val1-val2;
				dist+=Math.sqrt((Math.pow(diff, 2)));
				count++;
			}
		}
		return dist/count;
	}
	
	
	/**
	 * Align the matrices by the maximum within each and then compare point-by-point distances
	 * @param observed
	 * @param expected
	 * @return
	 */
	private static double distance2(MatrixWithHeaders observed, MatrixWithHeaders expected) {
		for(String col: observed.getColumnNames()) {
			int timeAtMax_O=getTimeMax(observed, col);
			int timeAtMax_E=getTimeMax(expected, col);
			System.err.println(col+" "+timeAtMax_O+" "+timeAtMax_E);
		}
		return 0;
	}
	
	
	private static int getTimeMax(MatrixWithHeaders expected, String col) {
		String name=null;
		double max=-Double.MAX_VALUE;
		for(String row: expected.getRowNames()) {
			double val=expected.get(row, col);
			if(val>max) {
				max=val;
				name=row;
			}
		}
		return Integer.parseInt(name);
	}


	private static void BEDToAnnotationFormat(Gene transcript) {
		System.out.println("# "+transcript.getName());
		//3640
		System.out.println(transcript.getGenomicLength());
		for(Annotation junction: transcript.getIntrons()) {
			int relativeStart=junction.getReferenceStartPosition()-transcript.getReferenceStartPosition();
			int relativeEnd=junction.getReferenceEndPosition()-transcript.getReferenceStartPosition();
			System.out.println(relativeStart+" "+relativeEnd+" "+junction.getReferenceName()+":"+junction.getReferenceStartPosition()+"-"+junction.getReferenceEndPosition());
		}
	}
	
	
	private static double[] fit(MatrixWithHeaders observed, double[] initiationRange, double[] splicingRange, double[] decayRange, File annotation) throws IOException {
		double absoluteMin=Double.MAX_VALUE;
		double[] minParams= {random(initiationRange), random(splicingRange), random(decayRange), Double.MAX_VALUE};
		
		for(int i=0; i<numIter; i++) {
			double[] minFit=optimize(observed, initiationRange, splicingRange, decayRange, minParams, annotation);
			//MatrixWithHeaders expected=simulate(minFit, annotation);
			double distance=minFit[3];
			System.err.println("iter "+i+" "+minFit[0]+" "+minFit[1]+" "+minFit[2]+" "+distance);
			if(distance<absoluteMin) {
				minParams=minFit;
				absoluteMin=distance;
			}
		}
		
		return minParams;
	}
	
	
	
	
	private static double[] scaleAndRange(MatrixWithHeaders observed, File annotation) throws IOException {
		double absoluteMin=Double.MAX_VALUE;
		
		
		//start range in integer space
		//find mins then subdivide
		
		double scale=1;
		double[] initiationRange= {1,2,3,4,5,6,7,8,9,10};
		double[] splicingRange= {1,2,3,4,5,6,7,8,9,10};
		double[] decayRange= {1,2,3,4,5,6,7,8,9,10};
		double[] minParams= {random(initiationRange), random(splicingRange), random(decayRange), Double.MAX_VALUE};
		
		for(int i=0; i<numIter; i++) {
			scale=scale/10.0;
			double[] minFit=optimize(observed, initiationRange, splicingRange, decayRange, minParams, annotation);
			initiationRange=updateRange(minFit[0], scale);
			splicingRange=updateRange(minFit[1], scale);
			decayRange=updateRange(minFit[2], scale);
			double distance=minFit[3];
			System.err.println("iter "+i+" "+minFit[0]+" "+minFit[1]+" "+minFit[2]+" "+distance);
			if(distance<absoluteMin) {
				minParams=minFit;
				absoluteMin=distance;
			}
		}
		
		return minParams;
	}
	
	
	private static double[] updateRange(double min, double scale) {
		double[] vals=new double[20];
		
		int counter=0;
		for(int i=-10; i<10; i++) {
			double val=min+(i*scale);
			System.err.println(min+" "+val);
			vals[counter]=val;
			counter++;
		}
		
		return vals;
	}


	private static double random(double[] initiationRange) {
		int index=new Double(Math.random()*initiationRange.length).intValue();
		return initiationRange[index];
	}


	private static double[] optimize(MatrixWithHeaders observed, double[] initiationRange, double[] splicingRange, double[] decayRange, double[] initialConditions, File annotation) throws IOException {
		
		double min=Double.MAX_VALUE;
		for(int i=0; i<initiationRange.length; i++) {
			MatrixWithHeaders expected=simulate(initiationRange[i], initialConditions[1], initialConditions[2], elongationRate, annotation);
			double distance=Double.MAX_VALUE;
			if(expected!=null) {distance=distance(observed, expected);}
			System.out.println("initiation "+initiationRange[i]+" "+initialConditions[1]+ " "+initialConditions[2]+" "+distance);
			if(distance<min) {
				min=distance;
				initialConditions[0]=initiationRange[i];
			}
		}
		initialConditions[3]=min;
		
		min=Double.MAX_VALUE;
		for(int i=0; i<decayRange.length; i++) {
			MatrixWithHeaders expected=simulate(initialConditions[0], decayRange[i], initialConditions[2], elongationRate, annotation);
			double distance=Double.MAX_VALUE;
			if(expected!=null) {
				distance=distance(observed, expected);
			}
			System.out.println("decay "+initialConditions[0]+" "+ decayRange[i]+ " "+initialConditions[2]+" "+distance);
			if(distance<min) {
				min=distance;
				initialConditions[1]=decayRange[i];
			}
		}
		initialConditions[3]=min;
		
		min=Double.MAX_VALUE;
		for(int i=0; i<splicingRange.length; i++) {
			MatrixWithHeaders expected=simulate(initialConditions[0], initialConditions[1], splicingRange[i], elongationRate, annotation);
			double distance=Double.MAX_VALUE;
			if(expected!=null) {
				distance=distance(observed, expected);
			}
			System.out.println("splicing "+initialConditions[0]+" "+initialConditions[1]+" "+ splicingRange[i]+" "+distance);
			if(distance<min) {
				min=distance;
				initialConditions[2]=splicingRange[i];
			}
		}
		initialConditions[3]=min;
		
		if(initialConditions[3]==Double.MAX_VALUE) {
			initialConditions[0]= initiationRange[0];
			initialConditions[1]=splicingRange[0];
			initialConditions[2]=decayRange[0];
		}
		
		return initialConditions;
		
	}


	public static void main(String[] args) throws IOException {
		/**
		 * Simulate data matrix
		 */
		if(args.length>1) {
		double initiationRate=Double.parseDouble(args[0]);
		double splicingRate=Double.parseDouble(args[1]);
		double decayRate=Double.parseDouble(args[2]);
		double elongationRate=Double.parseDouble(args[3]);
		File annot=new File(args[4]);
		
		FitBestParams.timeShift=Integer.parseInt(args[5]);
		
		MatrixWithHeaders simulated=simulate(initiationRate, splicingRate, decayRate, elongationRate, annot);
		
		simulated.write(args[6]);
		
		MatrixWithHeaders observed=new MatrixWithHeaders(new File(args[7]));
		distance2(observed, simulated);
		
		}else{System.err.println(" args[0]=1.0 0.001 0.005 50.0 \n args[1]=annotation \n args[2]=label time \n args[3]=save");}
		
		
		/**
		 * Generate annotation files
		 */
		
		
		/*if(args.length>1) {
		MatrixWithHeaders observed=new MatrixWithHeaders(new File(args[0])); //1.0, 0.001, 0.005, 5.0
		
		File annot=new File(args[1]);
		double[] initiationRate= {0.001, 0.01, 0.1};
		double[] decay= {0.000001, 0.00001, 0.0001, 0.001,0.01, 0.1};
		double[] splicing= {0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1};
		
		double[] fit=fit(observed, initiationRate, decay, splicing, annot);
		
		
		
		System.err.println(fit[0]+" "+fit[1]+" "+fit[2]+" "+fit[3]);
		
		
		MatrixWithHeaders simulated=simulate(fit[0], fit[1], fit[2], elongationRate, annot);
	
		simulated.write("sim.matrix");
		}else {System.err.println(" args[0]=observed matrix \n args[1]=annotation");}*/
		
	}
	
	
}

package guttmanlab.core.test;

import flanagan.physchem.ImmunoAssay;

public class FitTest {

	public static void main(String[] args) {
		double[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		double[] y= {0.148,0,0.054,0.02,0.286,0.188,0.36,0.632,1,0.88,0.973,1};
		
		
		
		double[] fit=fit4PL(x,y);
		System.err.println(fit[0]+"\t"+fit[1]+"\t"+fit[2]+"\t"+fit[3]+"\t"+fit[4]);
	}
	
	private static double[] fit4PL(double[] x, double[] y) {
		try {
		ImmunoAssay assay = new ImmunoAssay(x, y);
		
		assay.suppressPrint();
		assay.suppressYYplot();
		
		try {
			assay.fourParameterLogisticFit();
		}catch(java.awt.HeadlessException ex) {}
		double[] vals=assay.getModelParameterValues(); ///top, bottom, C50, HillSlope
		
		//double ss=assay.getTotalSumOfWeightedSquares();
		double r2=assay.getCoefficientOfDetermination();
		
		double[] rtrn=new double[vals.length+1];
		for(int i=0; i<vals.length; i++) {
			rtrn[i]=vals[i];
		}
		rtrn[vals.length]=r2;
		return rtrn;
		}catch(java.lang.IllegalArgumentException ex) {return null;}
	}
	
}

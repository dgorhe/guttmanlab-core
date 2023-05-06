package guttmanlab.core.splicing.speckle;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.stat.inference.ChiSquareTest;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.math.Statistics;

public class FitODE {
	
	private static final FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-14, 100.0, 1.0e-12, 1.0e-12);
	//private static final FirstOrderIntegrator dp853 = new GillIntegrator(1.0e-12);

	
	//TODO Make more robust to individual outlier timepoint
	public FitODE(String input, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		List<String> lines=BEDFileIO.loadLines(input);
		int[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		for(String line: lines) {
			String gene=line.split(",")[0];
			double[] S=parseSpliced(line, x);
			double[] U=parseUnspliced(line, x);
			
			
			double sum=Statistics.sum(S)+Statistics.sum(U);
			if(Statistics.min(S)>0 && Statistics.min(U)>0) {
				double[] params=optimizeParams(x, S, U);
				System.err.println(gene+"\t"+sum+"\t"+params[5]+"\t"+params[6]);
				write(writer, gene, params, sum, line.split(",")[2]);
			}
			else {System.err.println("skipped "+gene+" "+sum);}
		}
		writer.close();
	}
	
	
	private double[] parseUnspliced(String line, int[] x) {
		String[] tokens=line.split(",");
		double[] rtrn=new double[x.length];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=Double.parseDouble(tokens[i+53])*1000000;
		}
		
		return rtrn;
	}


	private double[] parseSpliced(String line, int[] x) {
		String[] tokens=line.split(",");
		double[] rtrn=new double[x.length];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=(Double.parseDouble(tokens[i+42])*1000000);
		}
		
		return rtrn;
	}


	private void write(FileWriter writer, String gene, double[] params, double sum, String geneLength) throws IOException {
		writer.write(gene+"\t"+geneLength+"\t"+sum);
		for(int i=0; i<params.length; i++) {
			writer.write("\t"+params[i]);
		}
		writer.write("\n");
	}


	private static double[] optimizeParams(int[] x, double[] S, double[] U) {
		return optimizeParams(x, x, S, U, S[0], U[0]);
	}
	
	
	//TODO Do we need to update S0 and U0?
	private static double[] optimizeParams(int[] termTime, int[] x, double[] S, double[] U, double S0, double U0) {
		double[] initialParams= {0.1,0.1,0.02}; //TODO how to estimate
		
		int numIter=10;
		double prevSS=Double.MAX_VALUE;
		boolean started=false;
			
		int timeAlpha=15;
		int tBeta=75;
		
		for(int i=0; i<numIter; i++) {
			double alpha=optimizeAlpha(initialParams, timeAlpha, tBeta, x, S, U, S0, U0);
			initialParams[0]=alpha;
				
			double beta=optimizeBeta(initialParams, timeAlpha, tBeta, x, S, U, S0, U0);
			initialParams[1]=beta;
				
			double gamma=optimizeGamma(initialParams, timeAlpha, tBeta, x, S, U, S0, U0);
			initialParams[2]=gamma;
			
			//TODO optimize tAlpha, tBeta
			timeAlpha=optimizeTAlpha(initialParams, termTime, tBeta, x,S, U, S0, U0);
			tBeta=optimizeTBeta(initialParams, timeAlpha, termTime, x,S, U, S0, U0);
			
			
			double ss=getSumOfSquares(initialParams, termTime, tBeta, x, S, U, S0, U0);
			
			double diff=prevSS-ss;
			double percentError=(diff/prevSS)*100;
			//System.err.println(i+" "+ss +" "+alpha+" "+beta+" "+gamma+" "+diff+" "+percentError+" "+optimalSS);
			
			if(started && percentError<0.01) {break;}
			prevSS=ss;
			started=true;
		}
			
		//int time=getOptimalTime(initialParams, termTime, tBeta, x, S, U);
		
		double[] chiSquare=goodnessOfFit(initialParams, timeAlpha, tBeta, x, U, S, U0, S0);
		
			
		//getModelVals(initialParams, time, tBeta, U0, S0);
		
		
		double[] rtrn= {initialParams[0], initialParams[1], initialParams[2], timeAlpha, tBeta, prevSS, chiSquare[1], chiSquare[0]};
		return rtrn;
		
		
	}
	
	
	private static int optimizeTAlpha(double[] initialParams, int[] termTime, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		double minSS=Double.MAX_VALUE;
		int rtrn=Integer.MAX_VALUE;
		for(int k=0; k<termTime.length; k++) {
			int tAlpha=termTime[k];
			double ss=0;
			FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], tAlpha, tBeta);
			for(int i=0; i<x.length; i++) {
				double[] y = new double[] {S0, U0}; // initial state
				dp853.integrate(ode, 0, y, x[i], y);
				double diff1=Math.abs(y[0]-S[i]);
				double diff2=Math.abs(y[1]-U[i]);
				ss+=diff1;
				ss+=diff2;
			}
			if(ss<minSS) {
				minSS=ss;
				rtrn=tAlpha;
			}
			
		}
		return rtrn;
	}
	
	private static int optimizeTAlpha(double[] initialParams, List<Double> options, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		double minSS=Double.MAX_VALUE;
		int rtrn=Integer.MAX_VALUE;
		for(Double time: options) {
			int tAlpha=new Double(time).intValue();
			double ss=0;
			FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], tAlpha, tBeta);
			for(int i=0; i<x.length; i++) {
				double[] y = new double[] {S0, U0}; // initial state
				dp853.integrate(ode, 0, y, x[i], y);
				double diff1=Math.abs(y[0]-S[i]);
				double diff2=Math.abs(y[1]-U[i]);
				ss+=diff1;
				ss+=diff2;
			}
			if(ss<minSS) {
				minSS=ss;
				rtrn=tAlpha;
			}
			
		}
		return rtrn;
	}
	
	private static int optimizeTBeta(double[] initialParams, int tAlpha, List<Double> options, int[] x, double[] S, double[] U, double S0, double U0) {
		double minSS=Double.MAX_VALUE;
		int rtrn=Integer.MAX_VALUE;
		for(Double time: options) {
			int tBeta=new Double(time).intValue();
			double ss=0;
			FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], tAlpha, tBeta);
			for(int i=0; i<x.length; i++) {
				double[] y = new double[] {S0, U0}; // initial state
				dp853.integrate(ode, 0, y, x[i], y);
				double diff1=Math.abs(y[0]-S[i]);
				double diff2=Math.abs(y[1]-U[i]);
				ss+=diff1;
				ss+=diff2;
			}
			if(ss<minSS) {
				minSS=ss;
				rtrn=tBeta;
			}
			
		}
		return rtrn;
	}
	
	private static int optimizeTBeta(double[] initialParams, int tAlpha, int[] termTime, int[] x, double[] S, double[] U, double S0, double U0) {
		double minSS=Double.MAX_VALUE;
		int rtrn=Integer.MAX_VALUE;
		for(int k=0; k<termTime.length; k++) {
			int tBeta=termTime[k];
			double ss=0;
			FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], tAlpha, tBeta);
			for(int i=0; i<x.length; i++) {
				double[] y = new double[] {S0, U0}; // initial state
				dp853.integrate(ode, 0, y, x[i], y);
				double diff1=Math.abs(y[0]-S[i]);
				double diff2=Math.abs(y[1]-U[i]);
				ss+=diff1;
				ss+=diff2;
			}
			if(ss<minSS) {
				minSS=ss;
				rtrn=tBeta;
			}
			
		}
		return rtrn;
	}


	private static double getSumOfSquares(int[] x, double[] s, double[] u) {
		Map<Integer,Double> SMean=getReplicates(x, s);
		Map<Integer,Double> UMean=getReplicates(x, u);
		
		
		double ss=0;
		for(int i=0; i<x.length; i++) {
			int xVal=x[i];
			double diff1=Math.abs(s[i]-SMean.get(xVal));
			double diff2=Math.abs(u[i]-UMean.get(xVal));
			ss+=diff1;
			ss+=diff2;
		}
		
		return ss;
	}


	private static Map<Integer, Double> getReplicates(int[] x, double[] s) {
		Map<Integer, List<Double>> replicates=new TreeMap<Integer, List<Double>>();
		for(int i=0; i<x.length; i++) {
			int val=x[i];
			if(!replicates.containsKey(val)) {replicates.put(val, new ArrayList<Double>());}
			List<Double> vals=replicates.get(val);
			vals.add(s[i]);
			replicates.put(val, vals);
		}
		
		Map<Integer, Double> avg=new TreeMap<Integer, Double>();
		
		for(int val: replicates.keySet()) {
			List<Double> vals=replicates.get(val);
			double mean=Statistics.mean(vals);
			avg.put(val, mean);
		}
		
		return avg;
	}


	//TODO Do we need to update S0 and U0?
	private static double[] optimizeParams(int time, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		int[] termTime= {time};
		return optimizeParams(termTime, x, S, U, S0, U0);
			
	}
	
	private static double[] goodnessOfFit(double[] initialParams, int time, int tBeta, int[] x, double[] U, double[] S, double U0, double S0) {
		double ss=0;
		FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], time, tBeta);
		
		double[] observedS=new double[x.length-1];
		double[] observedU=new double[x.length-1];
		double[] expectedS=new double[x.length-1];
		double[] expectedU=new double[x.length-1];
		
		for(int i=1; i<x.length; i++) {
			double[] y = new double[] {S0, U0}; // initial state
			dp853.integrate(ode, 0, y, x[i], y);
			
			observedU[i-1]=y[1];
			observedS[i-1]=y[0];
			expectedU[i-1]=U[i];
			expectedS[i-1]=S[i];
			
			double diff1=Math.pow(y[0]-S[i], 2)/S[i];
			double diff2=Math.pow(y[1]-U[i], 2)/U[i];
			ss+=diff1;
			ss+=diff2;
		}
		
		//double df=(x.length*2)-1;
		//ChiSquaredDistribution chi=new ChiSquaredDistribution(df);
		
		double[] expected=mergeD(expectedS, expectedU);
		long[] observed=merge(observedS, observedU);
		
		double[]rtrn= {1.0, 0.0};
		
		if((Statistics.sum(observedU)+Statistics.sum(observedS))>0 && Statistics.sum(expected)>0) {
			ChiSquareTest test=new ChiSquareTest();
			double p=test.chiSquareTest(expected, observed);
			//double p2=1-chi.cumulativeProbability(ss);
		
			//System.err.println(ss+" "+df+" "+p+" "+p2);
			rtrn[0]= p;
			rtrn[1]=ss;
		}
		
		return rtrn;
	}


	private static double[] mergeD(double[] expectedS, double[] expectedU) {
		double[] rtrn=new double[expectedU.length+expectedS.length];
		
		for(int i=0; i<expectedS.length; i++) {
			rtrn[i]=expectedS[i];
		}
		for(int i=expectedS.length; i<rtrn.length; i++) {
			rtrn[i]=expectedU[i-expectedS.length];
		}
		
		return rtrn;
	}


	private static long[] merge(double[] observedS, double[] observedU) {
		long[] rtrn=new long[observedU.length+observedS.length];
		
		for(int i=0; i<observedS.length; i++) {
			rtrn[i]=new Double(observedS[i]).longValue();
		}
		for(int i=observedS.length; i<rtrn.length; i++) {
			rtrn[i]=new Double(observedU[i-observedS.length]).longValue();
		}
		
		return rtrn;
	}


	private static void getModelVals(double[] initialParams, int termTime, int tBeta, double U0, double S0) {
		FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], termTime, tBeta);
		double[] t= getStep(0, 300);
		
		for(int i=1; i<t.length; i++) {
			double[] y = new double[] {S0, U0}; // initial state
			dp853.integrate(ode, 0, y, t[i], y);
			System.out.println(t[i]+"\t"+y[0]+"\t"+y[1]+"\t"+(initialParams[2]*(y[1]+y[0])));
		}
		
	}
	
	private static double[] getStep(int start, int end) {
		double[] rtrn=new double[end-start];
		
		for(int i=start; i<end; i++) {
			rtrn[i-start]=i;
		}
		
		return rtrn;
	}
	
	private static int[] getStepInt(int start, int end) {
		int[] rtrn=new int[end-start];
		
		for(int i=start; i<end; i++) {
			rtrn[i-start]=i;
		}
		
		return rtrn;
	}
	
	private static int getOptimalTime(double[] initialParams, int[] termTime, int tBeta, int[] x, double[] S, double[] U) {
		//FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		double minSS=Double.MAX_VALUE;
		int rtrn=Integer.MAX_VALUE;
		for(int time=0; time<termTime.length; time++) {
				double ss=0;
				FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], termTime[time], tBeta);
				for(int i=1; i<x.length; i++) {
					double[] y = new double[] {S[0], U[0]}; // initial state
					dp853.integrate(ode, 0, y, x[i], y);
					double diff1=Math.abs(y[0]-S[i]);
					double diff2=Math.abs(y[1]-U[i]);
					ss+=diff1;
					ss+=diff2;
				}
				if(ss<minSS) {
					minSS=ss;
					rtrn=termTime[time];
				}
			}
		
		return rtrn;
	}
	
	
	private static double getSumOfSquares(double[] initialParams, int[] termTime, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		//FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		double minSS=Double.MAX_VALUE;
		
			for(int time=0; time<termTime.length; time++) {
				double ss=0;
				FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], termTime[time], tBeta);
				for(int i=0; i<x.length; i++) {
					double[] y = new double[] {S0, U0}; // initial state
					dp853.integrate(ode, 0, y, x[i], y);
					double diff1=Math.abs(y[0]-S[i]);
					double diff2=Math.abs(y[1]-U[i]);
					ss+=diff1;
					ss+=diff2;
				}
				if(ss<minSS) {
					minSS=ss;
				}
			}
	
		return minSS;
	}
	
	private static double optimizeBeta(double[] initialParams, int tAlpha, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		List<Double> options=getInitialGuesses(0,100,10);
		double beta=optimizeBeta(options, initialParams, tAlpha, tBeta, x, S, U, S0, U0);
		
		double prevStep=10;
		for(int i=0; i<=4; i++) {
			double step=1.0/Math.pow(10, i);
			//System.err.println(i+" "+step+" "+beta);
			options=update(beta, prevStep, step);
			beta=optimizeBeta(options, initialParams, tAlpha, tBeta, x, S, U, S0, U0);
			prevStep=step;
		}
		
		return beta;
	}
	
	private static double optimizeBeta(List<Double> betas, double[] initialParams, int tAlpha, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		for(Double b: betas) {
			double ss=0;
				FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], b, initialParams[2], tAlpha, tBeta);
				for(int i=0; i<x.length; i++) {
					double[] y = new double[] {S0, U0}; // initial state
					dp853.integrate(ode, 0, y, x[i], y);
					double diff1=Math.abs(y[0]-S[i]);
					double diff2=Math.abs(y[1]-U[i]);
					ss+=diff1;
					ss+=diff2;
				}
				if(ss<minSS) {
					minSS=ss;
					rtrn=b;
				}
			}
		return rtrn;
	}
	
	
	private static double optimizeGamma(double[] initialParams, int tAlpha, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		List<Double> options=getInitialGuesses(0,100,10);
		double gamma=optimizeGamma(options, initialParams, tAlpha, tBeta, x, S, U, S0, U0);
		
		double prevStep=10;
		for(int i=0; i<=4; i++) {
			double step=1.0/Math.pow(10, i);
			//System.err.println(i+" "+step+" "+gamma);
			options=update(gamma, prevStep, step);
			gamma=optimizeGamma(options, initialParams, tAlpha, tBeta, x, S, U, S0, U0);
			prevStep=step;
		}
		
		return gamma;
	}
	
	private static double optimizeTAlpha(double[] initialParams, int tAlpha, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		List<Double> options=getInitialGuesses(0,300,10);
		tAlpha=optimizeTAlpha(initialParams, options, tBeta, x, S, U, S0, U0);
		
		double prevStep=10;
		for(int i=0; i<=0; i++) {
			double step=1.0/Math.pow(10, i);
			options=update(tAlpha, prevStep, step);
			tAlpha=optimizeTAlpha(initialParams, options, tBeta, x, S, U, S0, U0);
			prevStep=step;
		}
		
		return tAlpha;
	}
	
	private static double optimizeTBeta(double[] initialParams, int tAlpha, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		List<Double> options=getInitialGuesses(0,300,10);
		tBeta=optimizeTBeta(initialParams, tAlpha, options, x, S, U, S0, U0);
		
		double prevStep=10;
		for(int i=0; i<=0; i++) {
			double step=1.0/Math.pow(10, i);
			options=update(tAlpha, prevStep, step);
			tBeta=optimizeTBeta(initialParams, tAlpha, options, x, S, U, S0, U0);
			prevStep=step;
		}
		
		return tBeta;
	}
	
	
	private static double optimizeGamma(List<Double> gammas, double[] initialParams, int tAlpha, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		for(Double g: gammas) {
			double ss=0;
				FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], g, tAlpha, tBeta);
				for(int i=0; i<x.length; i++) {
					double[] y = new double[] {S0, U0}; // initial state
					dp853.integrate(ode, 0, y, x[i], y);
					double diff1=Math.abs(y[0]-S[i]);
					double diff2=Math.abs(y[1]-U[i]);
					ss+=diff1;
					ss+=diff2;
				}
				if(ss<minSS) {
					minSS=ss;
					rtrn=g;
				}
			}
		return rtrn;
	}

	private static double optimizeAlpha(double[] initialParams, int tAlpha, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		List<Double> options=getInitialGuesses(0,100,10);
		double alpha=optimizeAlpha(options, initialParams, tAlpha, tBeta, x, S, U, S0, U0);
		
		double prevStep=10;
		for(int i=0; i<=4; i++) {
			double step=1.0/Math.pow(10, i);
			//System.err.println(i+" "+step+" "+alpha);
			options=update(alpha, prevStep, step);
			alpha=optimizeAlpha(options, initialParams, tAlpha, tBeta, x, S, U, S0, U0);
			prevStep=step;
		}
		
		return alpha;
	}
	
	private static double optimizeAlpha(List<Double> alphas, double[] initialParams, int tAlpha, int tBeta, int[] x, double[] S, double[] U, double S0, double U0) {
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		for(Double a: alphas) {
			double ss=0;
			FirstOrderDifferentialEquations ode = new DecayODE(a, initialParams[1], initialParams[2], tAlpha, tBeta);
			for(int i=0; i<x.length; i++) {
				double[] y = new double[] {S0, U0}; // initial state
				dp853.integrate(ode, 0, y, x[i], y);
				double diff1=Math.abs(y[0]-S[i]);
				double diff2=Math.abs(y[1]-U[i]);
				ss+=diff1;
				ss+=diff2;
			}
			if(ss<minSS) {
				minSS=ss;
				rtrn=a;
			}
			
		}
		return rtrn;
	}
	
	
	private static List<Double> getInitialGuesses(double start, double end, double step) {
		List<Double> rtrn=new ArrayList<Double>();
		for(double val=start; val<=end; val+=step) {
			rtrn.add(val);
		}
		return rtrn;
	}

	private static List<Double> update(double alpha, double prevStep, double newStep) {
		List<Double> rtrn=new ArrayList<Double>();
		for(double val=Math.max(0,alpha-prevStep); val<=alpha+prevStep; val+=newStep) {
			rtrn.add(val);
		}
		return rtrn;
	}

	
	
	
	private static class DecayODE implements FirstOrderDifferentialEquations {

		double alpha1;
		int transcriptionTermTime;
		double beta1;
		double gamma;
		int tBeta;
		
		public DecayODE(double alpha1, double beta, double gamma, int termTime, int tBeta) {
			this.alpha1=alpha1;
			this.transcriptionTermTime=termTime;
			this.beta1=beta;
			this.gamma=gamma;
			this.tBeta=tBeta;
		}
		
		@Override
		public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
			double alpha=0;
			double beta=0;
			if(t<transcriptionTermTime) {alpha=alpha1;}
			if(t>=tBeta) {beta=beta1;}
			yDot[0]=gamma*y[1]-(beta*y[0]); //dS/dt
			yDot[1]=alpha-(gamma*y[1]); //dU/dt
		}

		@Override
		public int getDimension() {
			return 2;
		}	
	}
	
	private static int[] repeat(int[] x, int n) {
		int[] rtrn=new int[x.length*n];
		
		int counter=0;
		for(int i=0; i<x.length; i++) {
			for(int j=0; j<n; j++) {
				rtrn[counter]=x[i];
				counter++;
			}
		}
		
		return rtrn;
	}
	
	private static double[] parse(String line) throws IOException {
		String[] tokens=line.split("\t");
		double[] rtrn=new double[tokens.length-11];
		
		for(int i=11; i<tokens.length; i++) {
			rtrn[i-11]=Double.parseDouble(tokens[i]);
		}
		
		return rtrn;
	}
	
	private static double parseInitial(String line) throws IOException {
		String[] tokens=line.split("\t");
		double[] rtrn=new double[10];
		
		for(int i=1; i<=10; i++) {
			rtrn[i-1]=Double.parseDouble(tokens[i]);
		}
		
		return Statistics.mean(rtrn);
	}
	
	public static void main(String[] args) throws IOException {
		/*int[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		double[] S= {19.4,20.1,22.7,18,18.3,23.7,46.1,54.7,68.4,70.6,77.3,112.181};
		double[] U= {15.7,19.7,29.1,30.5,24.4,28,14.7,6.92,4.39,2.62,1.15,2.59};
		double[] params=optimizeParams(x, S, U);
		System.err.println(params[0]+"\t"+params[1]+"\t"+params[2]+"\t"+params[3]+"\t"+params[4]+"\t"+params[5]+"\t"+params[6]);*/
		
		
		List<String> splicedLines=BEDFileIO.loadLines(args[0],1);
		List<String> unsplicedLines=BEDFileIO.loadLines(args[1],1);
		
		int[] truncatedTermTime= {10,15,20,25,30,45,60,75,90,120,240};
		//int[] termTime= {0, 10,15,20,25,30,45,60,75,90,120,240};
		int[] termTime= getStepInt(0, 300);
		int[] x=repeat(truncatedTermTime, 10);
		
		for(int i=0; i<splicedLines.size(); i++) {
			String splicedLine=splicedLines.get(i);
			String unsplicedLine=unsplicedLines.get(i);
			String name=splicedLine.split("\t")[0];
			double S0=parseInitial(splicedLine);
			double U0=parseInitial(unsplicedLine);
			double[] S=parse(splicedLine);
			double[] U=parse(unsplicedLine);
			if(Statistics.min(S)>0 && Statistics.min(U)>0) {
				double[] params=optimizeParams(termTime, x, S, U, S0, U0);
				double optimalSS=getSumOfSquares(x,S, U);
				System.out.println(name+"\t"+params[0]+"\t"+params[1]+"\t"+params[2]+"\t"+params[3]+"\t"+params[4]+"\t"+params[5]+"\t"+params[6]+"\t"+optimalSS);
			}
			else {System.err.println("skipped "+name);}
		}
		
		//getModelVals(params, new Double(params[3]).intValue(), new Double(params[4]).intValue(), U0, S0);
	}


	


	
	
}

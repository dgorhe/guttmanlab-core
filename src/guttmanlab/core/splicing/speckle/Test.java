package guttmanlab.core.splicing.speckle;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

public class Test {

	
	private static class DecayODE implements FirstOrderDifferentialEquations {

		double alpha1;
		int transcriptionTermTime;
		double beta;
		double gamma;
		
		public DecayODE(double alpha1, double beta, double gamma, int termTime) {
			this.alpha1=alpha1;
			this.transcriptionTermTime=termTime;
			this.beta=beta;
			this.gamma=gamma;
		}
		
		@Override
		public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
			double alpha=0;
			if(t<transcriptionTermTime) {alpha=alpha1;}
			yDot[0]=gamma*y[1]-(beta*y[0]); //dS/dt
			yDot[1]=alpha-(gamma*y[1]); //dU/dt
		}

		@Override
		public int getDimension() {
			return 2;
		}
		
	}
	
	
	private static class MatureODE implements FirstOrderDifferentialEquations {

		double beta;
		double gamma;
		
		public MatureODE(double beta, double gamma) {
			this.beta=beta;
			this.gamma=gamma;
		}
		
		@Override
		public void computeDerivatives(double t, double[] y, double[] yDot) throws MaxCountExceededException, DimensionMismatchException {
			yDot[0]=gamma*y[1]-(beta*y[0]); //dS/dt
			//yDot[1]=alpha-(gamma*y[1]); //dU/dt
		}

		@Override
		public int getDimension() {
			return 1;
		}
		
	}
	
	private static double[] getStep(int start, int end) {
		double[] rtrn=new double[end-start];
		
		for(int i=start; i<end; i++) {
			rtrn[i-start]=i;
		}
		
		return rtrn;
	}
	
	private static List<Double> getStep(double start, double end, double step) {
		List<Double> rtrn=new ArrayList<Double>();
		
		for(double i=start; i<end; i+=step) {
			rtrn.add(i);
		}
		
		return rtrn;
	}

	private static double sumOfSquares(int[] x, double[] s, double[] u, double[] y) {
		double rtrn=0;
		for(int i=1; i<x.length; i++) {
			double diff1=Math.abs(y[0]-s[i]);
			double diff2=Math.abs(y[1]-u[i]);
			rtrn+=diff1;
			rtrn+=diff2;
		}
		return rtrn;
	}
	
	

	private static double optimizeAlpha(double[] initialParams, int[] termTime, int[] x, double[] S, double[] U) {
		List<Double> options=getInitialGuesses(0,100,10);
		double alpha=optimizeAlpha(options, initialParams, termTime, x, S, U);
		
		double prevStep=10;
		for(int i=0; i<=4; i++) {
			double step=1.0/Math.pow(10, i);
			options=update(alpha, prevStep, step);
			alpha=optimizeAlpha(options, initialParams, termTime, x, S, U);
			prevStep=step;
		}
		
		return alpha;
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

	private static double optimizeAlpha(List<Double> alphas, double[] initialParams, int[] termTime, int[] x, double[] S, double[] U) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		for(Double a: alphas) {
			for(int time=1; time<termTime.length; time++) {
				double ss=0;
				FirstOrderDifferentialEquations ode = new DecayODE(a, initialParams[1], initialParams[2], termTime[time]);
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
					rtrn=a;
				}
			}
		}
		return rtrn;
	}
	
	private static int getOptimalTime(double[] initialParams, int[] termTime, int[] x, double[] S, double[] U) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		double minSS=Double.MAX_VALUE;
		int rtrn=Integer.MAX_VALUE;
		for(int time=1; time<termTime.length; time++) {
				double ss=0;
				FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], termTime[time]);
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
	
	
	private static double optimizeBeta(List<Double> betas, double[] initialParams, int[] termTime, int[] x, double[] S, double[] U) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		//List<Double> beta=getStep(0.0001,1,.001); //TODO How to estimate
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		for(Double b: betas) {
			for(int time=1; time<termTime.length; time++) {
				double ss=0;
				FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], b, initialParams[2], termTime[time]);
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
					rtrn=b;
				}
			}
		}
		return rtrn;
	}
	
	
	private static double optimizeBeta(double[] initialParams, int[] termTime, int[] x, double[] S, double[] U) {
		List<Double> options=getInitialGuesses(0,100,10);
		double beta=optimizeBeta(options, initialParams, termTime, x, S, U);
		
		double prevStep=10;
		for(int i=0; i<=4; i++) {
			double step=1.0/Math.pow(10, i);
			options=update(beta, prevStep, step);
			beta=optimizeBeta(options, initialParams, termTime, x, S, U);
			prevStep=step;
		}
		
		return beta;
	}
	
	private static double optimizeGamma(double[] initialParams, int[] termTime, int[] x, double[] S, double[] U) {
		List<Double> options=getInitialGuesses(0,100,10);
		double gamma=optimizeGamma(options, initialParams, termTime, x, S, U);
		
		double prevStep=10;
		for(int i=0; i<=4; i++) {
			double step=1.0/Math.pow(10, i);
			options=update(gamma, prevStep, step);
			gamma=optimizeGamma(options, initialParams, termTime, x, S, U);
			prevStep=step;
		}
		
		return gamma;
	}
	
	
	private static double optimizeBeta(double[] initialParams, int[] x, double[] S) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		List<Double> beta=getStep(0.0001,1,.001); //TODO How to estimate
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		for(Double b: beta) {
			double ss=0;
			FirstOrderDifferentialEquations ode = new MatureODE(b, initialParams[1]);
				for(int i=1; i<x.length; i++) {
					double[] y = new double[] {S[0]}; // initial state
					dp853.integrate(ode, 0, y, x[i], y);
					double diff1=Math.abs(y[0]-S[i]);
					//double diff2=Math.abs(y[1]-U[i]);
					ss+=diff1;
					//ss+=diff2;
				}
				if(ss<minSS) {
					minSS=ss;
					rtrn=b;
				}
			
		}
		return rtrn;
	}
	
	private static double optimizeGamma(List<Double> gammas, double[] initialParams, int[] termTime, int[] x, double[] S, double[] U) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		//List<Double> gamma=getStep(0.001,1,0.001); //TODO How to estimate
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		for(Double g: gammas) {
			for(int time=1; time<termTime.length; time++) {
				double ss=0;
				FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], g, termTime[time]);
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
					rtrn=g;
				}
			}
		}
		//System.err.println(minSS);
		return rtrn;
	}
	
	private static double optimizeGamma(double[] initialParams, int[] x, double[] S) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		List<Double> gamma=getStep(0.01,5,0.01); //TODO How to estimate
		double minSS=Double.MAX_VALUE;
		double rtrn=Double.MAX_VALUE;
		for(Double g: gamma) {
			//for(int time=1; time<termTime.length; time++) {
				double ss=0;
				FirstOrderDifferentialEquations ode = new MatureODE(initialParams[0], g);
				for(int i=1; i<x.length; i++) {
					double[] y = new double[] {S[0]}; // initial state
					dp853.integrate(ode, 0, y, x[i], y);
					double diff1=Math.abs(y[0]-S[i]);
					//double diff2=Math.abs(y[1]-U[i]);
					ss+=diff1;
					//ss+=diff2;
				}
				if(ss<minSS) {
					minSS=ss;
					rtrn=g;
				}
			//}
		}
		//System.err.println(minSS);
		return rtrn;
	}
	
	
	private static double getSumOfSquares(double[] initialParams, int[] termTime, int[] x, double[] S, double[] U) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		double minSS=Double.MAX_VALUE;
		
			for(int time=1; time<termTime.length; time++) {
				double ss=0;
				FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], termTime[time]);
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
				}
			}
	
		return minSS;
	}
	
	private static double getSumOfSquares(double[] initialParams, int[] x, double[] S) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		
		double ss=0;
		FirstOrderDifferentialEquations ode = new MatureODE(initialParams[0], initialParams[1]);
			for(int i=1; i<x.length; i++) {
				double[] y = new double[] {S[0]}; // initial state
				dp853.integrate(ode, 0, y, x[i], y);
				double diff1=Math.abs(y[0]-S[i]);
				//double diff2=Math.abs(y[1]-U[i]);
				ss+=diff1;
				//ss+=diff2;
			}
				
			
	
		return ss;
	}
	
	private static double[] optimizeParams(int[] x, double[] S, double[] U) {
		int[] termTime= x;
		double[] initialParams= {0.1,0.1,0.02}; //TODO how to estimate
			
		int numIter=10;
		double prevSS=Double.MAX_VALUE;
			
		for(int i=0; i<numIter; i++) {
			
			double alpha=optimizeAlpha(initialParams, termTime, x, S, U);
			initialParams[0]=alpha;
				
			double beta=optimizeBeta(initialParams, termTime, x, S, U);
			initialParams[1]=beta;
				
			double gamma=optimizeGamma(initialParams, termTime, x, S, U);
			initialParams[2]=gamma;
				
			double ss=getSumOfSquares(initialParams, termTime, x, S, U);
			System.err.println(i+" "+initialParams[0]+" "+initialParams[1]+" "+initialParams[2]+" "+ss);
				
			if(Math.abs(prevSS-ss)<1) {break;}
			prevSS=ss;
		}
			
		int time=getOptimalTime(initialParams, termTime, x, S, U);
			
		getModelVals(initialParams, time, x, U, S);
			
		return initialParams;
		
	}
	
	
	public static void main(String[] args) {
		int[] x= {0,10,15,20,25,30,45,60,75,90,120,240};
		int[] termTime= x;
		
		
		double[] S= {
				19.4,
				20.1,
				22.7,
				18,
				18.3,
				23.7,
				46.1,
				54.7,
				68.4,
				70.6,
				77.3,
				112.181};
		
		
		double[] U= {
				15.7,
				19.7,
				29.1,
				30.5,
				24.4,
				28,
				14.7,
				6.92,
				4.39,
				2.62,
				1.15,
				2.59};
		
		
		double[] initialParams= {0.1,0.1,0.02}; //TODO how to estimate
		
		int numIter=10;
		double prevSS=Double.MAX_VALUE;
		
		for(int i=0; i<numIter; i++) {
		
			double alpha=optimizeAlpha(initialParams, termTime, x, S, U);
			initialParams[0]=alpha;
			
			double beta=optimizeBeta(initialParams, termTime, x, S, U);
			initialParams[1]=beta;
			
			double gamma=optimizeGamma(initialParams, termTime, x, S, U);
			initialParams[2]=gamma;
			
			double ss=getSumOfSquares(initialParams, termTime, x, S, U); //TODO Leave one point out fit
			System.err.println(i+" "+initialParams[0]+" "+initialParams[1]+" "+initialParams[2]+" "+ss);
			
			if(Math.abs(prevSS-ss)<0.1) {break;}
			prevSS=ss;
		}
		
		int time=getOptimalTime(initialParams, termTime, x, S, U);
		
		getModelVals(initialParams, time, x, U, S);
		
	}
	

	
	private static void getModelVals(double[] initialParams, int[] x, double[] S) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		FirstOrderDifferentialEquations ode = new MatureODE(initialParams[0], initialParams[1]);
		double[] t= getStep(0, 300);
		
		for(int i=1; i<t.length; i++) {
			double[] y = new double[] {S[0]}; // initial state
			dp853.integrate(ode, 0, y, t[i], y);
			System.out.println(t[i]+"\t"+y[0]);
		}
		
	}
	
	private static void getModelVals(double[] initialParams, int termTime, int[] x, double[] U, double[] S) {
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-10, 100.0, 1.0e-10, 1.0e-10);
		FirstOrderDifferentialEquations ode = new DecayODE(initialParams[0], initialParams[1], initialParams[2], termTime);
		double[] t= getStep(0, 300);
		
		for(int i=1; i<t.length; i++) {
			double[] y = new double[] {S[0], U[0]}; // initial state
			dp853.integrate(ode, 0, y, t[i], y);
			System.out.println(t[i]+"\t"+y[0]+"\t"+y[1]);
		}
		
	}
	
	
	
	
	
	
}

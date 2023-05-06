package guttmanlab.core.clap.old;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;
import jsc.distributions.Binomial;

public class ComputeSignificanceWindow {

	public ComputeSignificanceWindow(int sampleTotal, int inputTotal, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		double p=(double)sampleTotal/((double)inputTotal+(double)sampleTotal);
			
		//double p=.5;
		
		System.err.println(p);
			
		Map<Double, Double> xy=new TreeMap<Double, Double>();
			
			for(int n=1; n<10000; n++){
				
				Binomial b=new Binomial(n, p);
				double criticalVal=b.inverseCdf(0.99);
				double pval=b.cdf(criticalVal);
				double x=n;
				double y=criticalVal;
				double inputCount=n-criticalVal;	
				double enrichment=y/inputCount;
				
				double iNorm=(double)inputCount/(double)inputTotal;
				double sNorm=(double)criticalVal/(double)sampleTotal;
				double eNorm=sNorm/iNorm;
				double sumNorm=iNorm+sNorm;
				
				writer.write(x+"\t"+sumNorm+"\t"+eNorm+"\n");
				xy.put(x, y);
				//System.out.println(sampleCount+" "+inputCount+" "+n+" "+x+" "+y+" "+sc);
			}
			
			writer.close();
			
			//iterate k=sample count
			
			//

	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
		new ComputeSignificanceWindow(new Integer(args[0]),new Integer(args[1]), args[2]);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=sample total \n args[1]=input total \n args[2]=save";
}

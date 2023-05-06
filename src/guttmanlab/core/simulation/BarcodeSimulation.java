package guttmanlab.core.simulation;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class BarcodeSimulation {

	double ligationFreq=0.75;
	int numberOfRounds=8;
	int numSamples=100000;
	int clusterSize=10;
	
	public BarcodeSimulation(String save) throws IOException{
		FileWriter writer=new FileWriter(save);

		double fractionGreater=0;
		double all=0;
		
		for(int i=0; i<numSamples; i++){
			writer.write("p"+i);
			int counter=0;
			for(int j=0; j<clusterSize; j++){
				ArrayList<Integer> list=ligation();
				writer.write("\t"+list+"\t"+list.size());
				if(list.size()>=numberOfRounds){counter++;}
			}
			double ratio=(double)counter/(double)clusterSize;
			//System.err.println(counter+" "+clusterSize);
			if(ratio>0.6){fractionGreater++;}
			all++;
		
			writer.write("\t"+counter);
			writer.write("\n");
		}
		
		
		
		writer.close();
		System.err.println(fractionGreater/all);
		
	}
	
	private ArrayList<Integer> ligation(){
		ArrayList<Integer> rtrn=new ArrayList<Integer>();
		for(int i=0; i<numberOfRounds; i++){
			//for each round choose from uniform
			double random=Math.random();
			if(random<ligationFreq){rtrn.add(i);}
			//else{i++;} //skip next one
		}
		return rtrn;
	}
	
	
	/*private ArrayList<Integer> ligation(){
		ArrayList<Integer> rtrn=new ArrayList<Integer>();
		int last=-1;
		for(int i=0; i<numberOfRounds; i++){
			//for each round choose from uniform
			double random=Math.random();
			if(random<ligationFreq){
				rtrn.add(i);
				last=i;
			}
			
			//else{i++;} //skip next one
		}
		return rtrn;
	}*/
	
	public static void main(String[] args) throws IOException{
		new BarcodeSimulation("/Users/mguttman/Desktop/freq95.8rounds.txt");
	}
	
}

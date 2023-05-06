package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import guttmanlab.core.annotation.SingleInterval;

public class PullRNARNAInteraction {

	public PullRNARNAInteraction(BarcodingDataStreaming data, String rna, String save, boolean weighted) throws IOException{
		Map<SingleInterval, Double> rnaVals=data.getRNACounts(rna, weighted);
		//Map<String, Double> inputVals=data.getRNACounts("Input", weighted);
		
		write(save, rnaVals);
	}

	private void write(String save, Map<SingleInterval, Double> rnaVals) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		double sum=0;
		for(SingleInterval rna: rnaVals.keySet()){
			String name=rna.getName();
			double val=rnaVals.get(rna);
			writer.write(name+"\t"+rna.toUCSC()+"\t"+val+"\n");
			sum+=val;
		}
		
		System.err.println(sum);
		writer.close();
	}
	
	
	public static void main (String[] args) throws IOException{
		if(args.length>3){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String rna=args[1];
			String save=args[2];
			boolean weighted=new Boolean(args[3]);
			new PullRNARNAInteraction(data, rna, save, weighted);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=rna \n args[2]=save \n args[3]=weight by cluster size";
	
}

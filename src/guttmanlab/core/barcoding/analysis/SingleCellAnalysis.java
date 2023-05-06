package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.SingleInterval;

public class SingleCellAnalysis {
	
	static int resolution=1000000;

	public static void main (String[] args) throws IOException{
		
		if(args.length>2){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String saveDir=args[1];
			int minClusterSize=new Integer(args[2]);
			new File(saveDir).mkdir();
			
			
			int i=0;
			while(data.hasNext()){
				Cluster c=data.next();
				c.bin(resolution);
				if(c.getClusterSize()>minClusterSize){
					write(saveDir, i, c.getAllIntervals());
					Collection<String> chromosomes=c.getAllChromosomes();
					System.err.println(c.getClusterSize()+" "+chromosomes);
				}
				i++;
			}
			data.close();
		}
		else{
			System.err.println(line);
		}
	}
	
	static String line=" args[0]=clusters \n args[1]=save directory \n args[2]=min cluster size to consider";

	private static void write(String saveDir, int i, Collection<SingleInterval> allIntervals) throws IOException {
		FileWriter writerHuman=new FileWriter(saveDir+"/"+"F"+i+".human.bed");
		FileWriter writerMouse=new FileWriter(saveDir+"/"+"F"+i+".mouse.bed");
		FileWriter writerAll=new FileWriter(saveDir+"/"+"F"+i+".all.bed");
		
		for(SingleInterval interval: allIntervals){
			String chr=interval.getReferenceName();
			FileWriter writer=null;
			if(chr.contains("human")){
				writer=writerHuman;
				chr=chr.replaceAll("_human", "");
			}
			else if(chr.contains("mouse")){
				writer=writerMouse;
				chr=chr.replaceAll("_mouse", "");
			}
			else{
				writer=writerAll;
			}
			
			writer.write(chr+"\t"+interval.getReferenceStartPosition()+"\t"+interval.getReferenceEndPosition()+"\n");
			
			
		}
		writerHuman.close();
		writerMouse.close();
		writerAll.close();
	}	
}
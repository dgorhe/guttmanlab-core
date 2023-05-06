package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.barcoding.analysis.BarcodingData;
import guttmanlab.core.barcoding.analysis.Cluster;

public class ClusterToBEDGraph {

	public ClusterToBEDGraph(File barcodeFile, String save) throws IOException{
		BarcodingData data=new BarcodingData(barcodeFile);
		data=data.bin(100);
		FileWriter writer=new FileWriter(save);
		
		for(Cluster c: data.getClusters()){
			if(!c.isInterchromosomal() && c.getClusterSize()>1){
				Annotation bed=c.toAnnotation("C");
				if(bed.getGenomicLength()<1000000){
					writer.write(bed.toString()+"\n");
				}
			}
			//else{System.err.println(c);}
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		File barcodeFile=new File(args[0]);
		String save=args[1];
		new ClusterToBEDGraph(barcodeFile, save);
	}
	
}

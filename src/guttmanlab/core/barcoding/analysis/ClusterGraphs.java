package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.SingleInterval;

public class ClusterGraphs {

	public ClusterGraphs(File barcodeFile, int resolution, String save, int minClusterSize) throws IOException{
		BarcodingData data=new BarcodingData(barcodeFile);
		data=data.bin(resolution);
		
		//for each barcode
		//write a cytoscape file with weighted edges
		for(String barcode:data.getBarcodes()){
			Collection<SingleInterval> positions=data.getPositionsWithBarcode(barcode);
			if(positions.size()>minClusterSize){
			writeSif(save+"."+barcode, positions, data);
			}
		}
	}

	private void writeSif(String save, Collection<SingleInterval> positions, BarcodingData data) throws IOException {
		FileWriter writerSif=new FileWriter(save+".sif");
		FileWriter writer=new FileWriter(save+".score");
		FileWriter writer2=new FileWriter(save+".chr");
		writer.write("name\tscore\n");
		writer2.write("Node\t=\tchromosome\n");
		for(SingleInterval pos1: positions){
			writer2.write(pos1.toUCSC()+"\t=\t"+pos1.getReferenceName()+"\n");
			for(SingleInterval pos2: positions){
				InteractionScore score=data.getInteractionScore(pos1, pos2);
				writerSif.write(pos1.toUCSC()+"\tpp\t"+pos2.toUCSC()+"\n");
				writer.write(pos1.toUCSC()+" (pp) "+pos2.toUCSC()+"\t"+score.getScore()+"\n");
				
			}
		}
		writerSif.close();
		writer.close();
		writer2.close();
	}
	
	public static void main(String[] args)throws IOException{
		File data=new File(args[0]);
		int resolution=new Integer(args[1]);
		String save=args[2];
		int minClusterSize=new Integer(args[3]);
		new ClusterGraphs(data, resolution, save, minClusterSize);
	}
	
	
}

package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.SingleInterval;

public class MakeBEDGraphUnnormalized {

	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		
		String in=args[0];
		String gene=args[1];
		String out=in+"."+gene+".clusters";
		String cmd="grep "+gene+" > "+in+"."+gene+".clusters";
		
		Process p=Runtime.getRuntime().exec(cmd);
		p.getInputStream();
		p.waitFor();
		
		
		
		
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(out));
		
		String save=args[2];
		FileWriter writer=new FileWriter(save);
		Collection<Cluster> clusters=data.getRNAClusters(gene);
		
		for(Cluster c: clusters){
			for(SingleInterval region: c.getAllDNAIntervals()){
				writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\n");
			}
		}
		
		writer.close();
		data.close();
		
		
		cmd="sort -k 1,1 "+save+"> "+save+".sorted.bed";
		p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		cmd=" bedtools genomecov -bg -i "+ save+".sorted.bed /groups/guttman/scripts/mm10.chrom.sizes > "+save+".bedgraph";
		p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
	}
	
}

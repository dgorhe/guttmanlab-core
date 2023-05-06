package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class MergeTables {

	public MergeTables(File bed1, File bed2, File bed3, String save) throws IOException{
		Map<SingleInterval, Double> map1=BEDFileIO.loadbedgraph(bed1);
		Map<SingleInterval, Double> map2=BEDFileIO.loadbedgraph(bed2);
		TreeMap<SingleInterval, Double> map3=BEDFileIO.loadbedgraph(bed3);
		
		write(map1, map2, map3, save);
		
	}
	
	private void write(Map<SingleInterval, Double> map1, Map<SingleInterval, Double> map2, TreeMap<SingleInterval, Double> map3, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map1.keySet()){
			if(map2.containsKey(region)){
				double score1=map1.get(region);
				double score2=map2.get(region);
				SingleInterval mbRegion=region.bin(1000000);
				if(map3.containsKey(mbRegion)){
					double score3=map3.get(mbRegion);
					writer.write(region.toUCSC()+"\t"+score1+"\t"+score2+"\t"+score3+"\n");
				}
				//else{System.err.println(mbRegion.toUCSC());}
			}
		}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		new MergeTables(new File(args[0]), new File(args[1]), new File(args[2]), args[3]);
	}
	
}

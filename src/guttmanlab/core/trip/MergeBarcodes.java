package guttmanlab.core.trip;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.sequence.Sequence;

public class MergeBarcodes {

	public static void main(String[] args) throws IOException {
		if(args.length>3) {
			/*List<String> cDNA=BEDFileIO.loadLines(args[0]);
			List<String> gDNA=BEDFileIO.loadLines(args[1]);
			
			Collection<String> allBarcodes=getAll(cDNA, gDNA);
			Collection<String> mergedBarcodes=CollapseCDNA.collapseBarcodes(allBarcodes, Integer.parseInt(args[3]));
			
			for(String mergedBarcode: mergedBarcodes) {
				Collection<SingleInterval> positions=new TreeSet<SingleInterval>();
				String[] barcodes=mergedBarcode.split(" ");
				for(int i=0; i<barcodes.length; i++) {
					String barcode
				}
				
			}*/
			
			
			
			List<String> cDNA=BEDFileIO.loadLines(args[0]);
			List<String> gDNA=BEDFileIO.loadLines(args[1]);
			
			Map<SingleInterval, Double> speckleScores=BEDFileIO.loadbedgraph(new File(args[2]));
			int resolution=1000;
			
			System.err.println("Resolution "+resolution);
			
			int counter=0;
			for(String genomic: gDNA) {
				String[] gDNABarcodes=genomic.split("\t")[0].split(" ");
				Map<SingleInterval, Double> regions=parse(genomic.split("\t")[2], resolution, speckleScores);
				
				for(String c: cDNA) {
					String[] cDNABarcodes=c.split("\t")[0].split(" ");
					if(overlaps(gDNABarcodes, cDNABarcodes)) {
						String line=genomic+"\t"+c+"\t"+regions.size();
						for(SingleInterval region: regions.keySet()) {
							line+="\t"+region.toUCSC()+" "+regions.get(region);
						}
						System.out.println(line);
					}
					
					
				}
				counter++;
				if(counter%100==0) {System.err.println(counter+" "+gDNA.size());}
			}
			
			
			
		}else {System.err.println(usage);}
	}

	static String usage=" args[0]=cDNA \n args[1]=gDNA \n args[2]=speckle distances (bedgraph, 100Kb) \n args[3]=hamming distance";
	
	
	
	
	
	private static Map<SingleInterval, Double> parse(String pos, int resolution, Map<SingleInterval, Double> speckleScores) {
		SingleInterval r=new SingleInterval(pos);
		
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		
		SingleInterval bin=r.bin(resolution);
		double val=-1;
		if(speckleScores.containsKey(bin)) {val=speckleScores.get(bin);}
		rtrn.put(bin, val);
		
		return rtrn;
		
	}

	private static boolean overlaps(String[] gDNABarcodes, String[] cDNABarcodes) {
		for(int i=0; i<gDNABarcodes.length; i++) {
			String barcode1=gDNABarcodes[i];
			for(int j=0; j<cDNABarcodes.length; j++) {
				String barcode2=cDNABarcodes[j];
				String antisense=Sequence.reverseComplement(barcode2);
				if(barcode1.equals(antisense)) {return true;}
			}
		}
		
		return false;
			
	}
	
}

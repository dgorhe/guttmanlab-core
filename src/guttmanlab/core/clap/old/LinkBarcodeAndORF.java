package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class LinkBarcodeAndORF {

	public LinkBarcodeAndORF(File bamFile, String save) throws IOException{
		
		SamReader lines=SamReaderFactory.makeDefault().open((bamFile));
		
		Map<String, Collection<SingleInterval>> map=new TreeMap<String, Collection<SingleInterval>>();
		
		for(SAMRecord line: lines){
			String barcode=line.getReadName().split("_")[1];
			SingleInterval region=new SingleInterval(line.getReferenceName(), line.getAlignmentStart(), line.getAlignmentEnd());
			Collection<SingleInterval> list=new TreeSet<SingleInterval>();
			if(map.containsKey(barcode)){
				list=map.get(barcode);
			}
			list.add(region);
			map.put(barcode, list);
		}
		
		write(map, save);
		
	}	
	
	private void write(Map<String, Collection<SingleInterval>> map, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: map.keySet()){
			writer.write(barcode);
			Collection<SingleInterval> list=map.get(barcode);
			for(SingleInterval region: list){
				writer.write("\t"+region.toUCSC());
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		File bam=new File(args[0]);
		String save=args[1];
		new LinkBarcodeAndORF(bam, save);
	}
	
}

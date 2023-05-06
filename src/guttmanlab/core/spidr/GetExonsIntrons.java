package guttmanlab.core.spidr;

import java.io.IOException;
import java.util.Collection;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class GetExonsIntrons {

	private static void getExons(String file) throws IOException {
		Collection<SingleInterval> exons=BEDFileIO.getExons(file);
		
		
		
		for(SingleInterval exon: exons) {
			System.out.println(exon.toShortBED());
		}
	}
	
	
	private static void getGenePosition(String file) throws IOException {
		Collection<SingleInterval> exons=BEDFileIO.getGeneInterval(file);
		exons=BEDFileIO.collapse(exons);
		
		for(SingleInterval exon: exons) {
			if(!exon.getReferenceName().contains("_")) {
				System.out.println(exon.toShortBED());
			}
		}
	}
	
	private static void getIntrons(String file) throws IOException {
		Collection<SingleInterval> introns=BEDFileIO.getIntrons(file);
		introns=BEDFileIO.collapse(introns);
		
		
		for(SingleInterval intron: introns) {
			if(!intron.getReferenceName().contains("_")) {
				System.out.println(intron.toShortBED());
			}
		}
	}
	
	public static void main(String[] args) throws IOException {
		getGenePosition(args[0]);
	}
	
}

package guttmanlab.core.spidr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class PullFeaturesFromMatrix {
	
	private static void pullExons(MatrixWithHeaders matrix, Collection<Gene> genes, int binSize, String save) throws IOException {
		Collection<SingleInterval> bins=new TreeSet<SingleInterval>();
		Collection<String> list=new TreeSet<String>();
		
		for(Gene gene: genes) {
			Collection<Annotation> introns= gene.getIntrons();
			
			
			//Iterator<SingleInterval> iter=gene.getBlocks();
			for(Annotation intron: introns) {
				SingleInterval exon=intron.getSingleInterval();
				Collection<SingleInterval> set=exon.allBins(binSize);
				
				
				for(SingleInterval bin: set) {
					bin.setName(gene.getName());
					String name=bin.toUCSC(Strand.antisense(bin.getOrientation()));
					//System.err.println(name);
					list.add(name);
					bins.add(bin);
				}
				
				
			}
		}
		
		//TODO write annotation file
		System.out.println("ID\tGeneName\tChr\tStart");
		for(SingleInterval bin: bins) {
			String id=bin.toUCSC(Strand.antisense(bin.getOrientation()));
			System.out.println(id+"\t"+bin.getName()+"\t"+bin.getReferenceName()+"\t"+bin.getReferenceStartPosition());
		}
		
		MatrixWithHeaders submatrix=matrix.submatrixByRowNames(list);
		submatrix.write(save);
	}

	

	private static void classify(List<String> rowNames, Collection<Gene> genes, int binSize, String save) throws IOException {
		Collection<String> exons=new TreeSet<String>();
		Collection<String> introns=new TreeSet<String>();
		Collection<String> utr5=new TreeSet<String>();
		Collection<String> utr3=new TreeSet<String>();
		Collection<String> cds=new TreeSet<String>();
		
		int count=0;
		for(Gene gene: genes) {
			utr3.addAll(bin(gene.get3UTR(), binSize));
			utr5.addAll(bin(gene.get5UTR(), binSize));
			cds.addAll(bin(gene.getCodingRegion(), binSize));
			
			exons.addAll(bin(gene, binSize));
			
			for(Annotation intron: gene.getIntrons()) {
				introns.addAll(bin(intron, binSize));
			}
			count++;
			if(count%100==0) {System.err.println(count);}
		}
		
		FileWriter writer=new FileWriter(save);
		writer.write("ID\tchromosome\tstart\tExon/Intron\tUTR/CDS\n");
		for(String row: rowNames) {
			String exonIntron=exonIntron(row, exons, introns);
			String utr=utr(row, utr5, utr3, cds);
			SingleInterval region=new SingleInterval(row);
			writer.write(row+"\t"+region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+exonIntron+"\t"+utr+"\n");
		}
		writer.close();
	}

	private static String utr(String row, Collection<String> utr5, Collection<String> utr3, Collection<String> cds) {
		if(utr5.contains(row)) {return "5UTR";}
		if(utr3.contains(row)) {return "3UTR";}
		if(cds.contains(row)) {return "CDS";}
		return "none";
	}

	private static String exonIntron(String row, Collection<String> exons, Collection<String> introns) {
		if(exons.contains(row) && introns.contains(row)) {return "both";}
		if(exons.contains(row)) {return "exon";}
		if(introns.contains(row)) {return "intron";}
		return "neither";
	}

	private static Collection<String> bin(Annotation region, int binSize) {
		Collection<String> rtrn=new TreeSet<String>();
		if(region!=null) {
			Iterator<SingleInterval> blocks=region.getBlocks();
			while(blocks.hasNext()) {
				SingleInterval exon=blocks.next();
				rtrn.addAll(binSI(exon, binSize));
			}
		}
		return rtrn;
	}

	private static Collection<String> binSI(SingleInterval exon, int binSize) {
		Collection<SingleInterval> set=exon.allBins(binSize);
		
		Collection<String> rtrn=new TreeSet<String>();
		for(SingleInterval bin: set) {
			String name=bin.toUCSC(Strand.antisense(bin.getOrientation()));
			rtrn.add(name);
		}
		return rtrn;
	}

	public static void main(String[] args) throws IOException {
		MatrixWithHeaders matrix=new MatrixWithHeaders(new File(args[0]));
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[1]);
		String save=args[2];
		
		int binSize=100;
		classify(matrix.getRowNames(), genes, binSize, save);
		
	}
	
	
}

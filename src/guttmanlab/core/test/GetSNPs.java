package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class GetSNPs {

	
	public GetSNPs(VCFFileReader vcf, Collection<Gene> genes, String sp1, String sp2, String save, int readLength) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Gene gene: genes){
			double fraction=getFractionSNPs(gene, vcf, sp1, sp2, readLength);
			System.err.println(gene.getName()+" "+gene.getGenomicLength()+" "+gene.toUCSC()+" "+fraction);
			writer.write(gene.getName()+"\t"+fraction+"\n");
		}
		
		/*CloseableIterator<VariantContext> iter=vcf.iterator();
		while(iter.hasNext()){
			VariantContext variant=iter.next();
			boolean isHet=isHetero(variant, sp1, sp2);
			if(isHet){
				writer.write(variant.getChr()+"\t"+variant.getStart()+"\t"+variant.getEnd()+"\n");
			}
		}
		
		iter.close();*/
		writer.close();
		vcf.close();
	}
	
	private double getFractionSNPs(Gene gene, VCFFileReader vcf, String sp1, String sp2, int readLength) {
		Map<SingleInterval, Integer> map=new TreeMap<SingleInterval, Integer>();
		IntervalTree<SingleInterval> windows=new IntervalTree<SingleInterval>();
		for(int i=0; i<gene.getGenomicLength(); i++){
			int start=gene.getReferenceStartPosition()+i;
			int end=start+readLength;
			SingleInterval window=new SingleInterval(gene.getReferenceName(), start, end);
			windows.put(start, end, window);
			map.put(window, 0);
		}
		
		
		CloseableIterator<VariantContext> iter=vcf.query(gene.getReferenceName().replace("chr", ""), gene.getReferenceStartPosition(), gene.getReferenceEndPosition());
		
		while(iter.hasNext()){
			VariantContext variant=iter.next();
			boolean isHet=isHetero(variant, sp1, sp2);
			if(isHet){
				Iterator<SingleInterval> witer=windows.overlappingValueIterator(variant.getStart(), variant.getEnd());
				while(witer.hasNext()){
					SingleInterval window=witer.next();
					map.put(window, 1);
				}
			}
		}
		
		iter.close();
		
		double sum=0;
		double count=0;
		for(SingleInterval region: map.keySet()){
			sum+=map.get(region);
			count++;
		}
		return sum/count;
	}

	private boolean isHetero(VariantContext variant, String sp1, String sp2) {
		return !variant.getGenotype(sp1).sameGenotype(variant.getGenotype(sp2));
	}

	public static void main(String[] args) throws IOException{
		if(args.length>5){
			
			VCFFileReader vcf=new VCFFileReader(new File(args[0]));
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[1]);
			String species1=args[2];
			String species2=args[3];
			String save=args[4];
			int readLength=new Integer(args[5]);
			new GetSNPs(vcf, genes, species1, species2, save, readLength);
		}
		else{System.err.println(usage);}
		
		
		
	}
	
	static String usage=" args[0]=VCF file \n args[1]=genes \n args[2]=species name 1 \n args[3]=species name 2 \n args[4]=save \n args[5]=read length";

	
	
}

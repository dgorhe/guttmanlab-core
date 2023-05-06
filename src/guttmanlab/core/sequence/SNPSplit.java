package guttmanlab.core.sequence;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import guttmanlab.core.annotation.SingleInterval;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class SNPSplit {

	public SNPSplit(VCFFileReader snps, String sp1, String sp2, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		CloseableIterator<VariantContext> iter=snps.iterator();
		
		//Map<SingleInterval, String> variants=new TreeMap<SingleInterval, String>();
		
		int counter=0;
		int hets=0;
		while(iter.hasNext()){
			VariantContext variant=iter.next();
			if(isHetero(variant, sp1, sp2)){
				SingleInterval region=new SingleInterval(variant.getChr(), variant.getStart(), variant.getEnd());
				String g1=variant.getGenotype(sp1).getGenotypeString();
				String g2=variant.getGenotype(sp2).getGenotypeString();
				//variants.put(region, g1+","+g2);
				writer.write("chr"+region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+g1+","+g2+"\n");
				hets++;
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter+" "+hets);}
		}
		
		writer.close();
		iter.close();
	}
	
	private boolean isHetero(VariantContext variant, String sp1, String sp2) {
		if(!variant.hasID()){return false;}
		if(!variant.hasGenotypes()){return false;}
		Genotype g1=variant.getGenotype(sp1);
		Genotype g2=variant.getGenotype(sp2);
		if(g1 ==null || g2==null){return false;}
		if(g1.getGenotypeString().equals("\\.") || g2.getGenotypeString().equals("\\.")){return false;}
		return !g1.sameGenotype(g2);
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>4){
			VCFFileReader vcf=new VCFFileReader(new File(args[0]));
			String save=args[1];
			String species1=args[2];
			String species2=args[3];
			new SNPSplit(vcf, species1, species2, save);
			vcf.close();
		}
		else{System.err.println(usage);}
		
		
	}
	
	static String usage=" args[0]=VCF file \n args[1]=save file name \n args[2]=species name 1 \n args[3]=species name 2";


	
	
}

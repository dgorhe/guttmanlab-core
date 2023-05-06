package guttmanlab.core.primer3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class GenerateSNPPrimers {
	String primer3core="/central/groups/guttman/software/arraydesign/primer3_core";
	int extension=40;
	
	public GenerateSNPPrimers(VCFFileReader vcf, Collection<Annotation> genes, Map<String, Sequence> sequenceByChr, String sp1, String sp2, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Primer3Configuration primer3config=Primer3ConfigurationFactory.getSNPPrimerConfiguration();
		
		Collection<SingleInterval> probes=new TreeSet<SingleInterval>();
		for(Annotation gene: genes){
			Collection<SingleInterval> temp=getProbes(gene, vcf, sp1, sp2);
			probes.addAll(temp);
			//System.err.println(gene.getName()+" "+temp.size());
			
		}
		
		for(SingleInterval probe: probes){
			Sequence probeSeq=sequenceByChr.get(probe.getReferenceName()).getSubsequence(probe);
			//String primer=getPrimers(probeSeq);
			PrimerPair primer=PcrPrimerDesigner.designPrimerPairAroundSNP(primer3config, probeSeq.getSequenceBases(), primer3core, extension-2, extension+2, 70);
			//PrimerPair primer=PcrPrimerDesigner.designBestPrimer(primer3config, probeSeq.getSequenceBases(), primer3core);
			if(primer!=null){writer.write(probe.getName()+"\t"+probe.toUCSC()+"\t"+primer.getLeftPrimer()+"\t"+primer.getRightPrimer()+"\t"+primer.getLeftPrimerTM()+"\t"+primer.getRightPrimerTM()+"\n");}
			//else{System.err.println("skipped "+probe.toUCSC());}
		}
		
		writer.close();
		vcf.close();
	}
	
	
	
	
	private String getPrimers(Sequence probeSeq) {
		String bases=probeSeq.getSequenceBases();
		
		String bestPrimer=null;
		double minDiff=Double.MAX_VALUE;
		for(int i=0; i<bases.length()-20; i++){
			String primer=bases.substring(i, i+20);
			double Tm=PrimerUtils.computeTM(primer);
			double diff=Math.abs(60-Tm);
			if(diff<minDiff){bestPrimer=primer; minDiff=diff;}
		}
		return bestPrimer;
		
	}




	private Collection<SingleInterval> getProbes(Annotation gene, VCFFileReader vcf, String sp1, String sp2) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		CloseableIterator<VariantContext> iter=vcf.query(gene.getReferenceName().replace("chr", ""), gene.getReferenceStartPosition(), gene.getReferenceEndPosition());
		
		while(iter.hasNext()){
			VariantContext variant=iter.next();
			boolean isHet=isHetero(variant, sp1, sp2);
			SingleInterval variantRegion=new SingleInterval(gene.getReferenceName(), variant.getStart(), variant.getEnd(), gene.getOrientation());
			if(isHet){
				//&& gene.overlaps(variantRegion)
				//generate probe
				String chr=gene.getReferenceName();
				int start=variant.getStart()-extension;
				int end=variant.getEnd()+extension;
				/*if(gene.getOrientation().equals(Strand.POSITIVE)){
					start=variant.getStart()-35;
					end=variant.getEnd();
				}*/
					
				SingleInterval region=new SingleInterval(chr, start, end, gene.getOrientation());
				region.setName(gene.getName());
				rtrn.add(region);
				
			}
		}
		
		iter.close();
		
		return rtrn;
	}
	
	private boolean isHetero(VariantContext variant, String sp1, String sp2) {
		return !variant.getGenotype(sp1).sameGenotype(variant.getGenotype(sp2));
	}

	
	private static Collection<Annotation> getExons(Collection<Gene> genes) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		for(Gene gene: genes){
			Iterator<SingleInterval> exon=gene.getBlocks();
			int index=1;
			while(exon.hasNext()){
				SingleInterval e=exon.next();
				e.setName(gene.getName()+"_Exon"+index);
				rtrn.add(e);
				index++;
			}
			//System.err.println(gene.getName()+" "+gene.getNumberOfBlocks());
			//rtrn.addAll(gene.getExonSet());
			//System.err.println(gene.getName()+" "+gene.getExonSet().size());
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>5){
			VCFFileReader vcf=new VCFFileReader(new File(args[0]));
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[1]);
			Map<String, Sequence> sequenceByChr=FastaFileIOImpl.readFromFileByName(args[2]);
			String species1=args[3];
			String species2=args[4];
			String save=args[5];
			
			
			Collection<Annotation> exons=getExons(genes);
			new GenerateSNPPrimers(vcf, exons, sequenceByChr, species1, species2, save);
		}
		else{System.err.println(usage);}
		
		
		
	}
	
	

	static String usage=" args[0]=VCF file \n args[1]=genes \n args[2]=genome fasta file \n args[3]=species name 1 \n args[4]=species name 2 \n args[5]=save";

	
	
}

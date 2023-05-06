package guttmanlab.core.probegeneration;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
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

public class GenerateSNPProbes {

	
	public GenerateSNPProbes(VCFFileReader vcf, Collection<Gene> genes, Map<String, Sequence> sequenceByChr, String sp1, String sp2, String save, int probeLength) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Annotation> probes=new TreeSet<Annotation>();
		for(Gene gene: genes){
			Collection<Annotation> temp=getProbes(gene, vcf, sp1, sp2, probeLength);
			probes.addAll(temp);
			System.err.println(gene.getName()+" "+temp.size());
		}
		
		for(Annotation probe: probes){
			Sequence probeSeq=sequenceByChr.get(probe.getReferenceName()).getSubsequence(probe);
			boolean rejectProbe=rejectProbe(probeSeq.getSequenceBases());
			if(!rejectProbe){
				writer.write(probe.getReferenceName()+"\t"+probe.getReferenceStartPosition()+"\t"+probe.getReferenceEndPosition()+"\t"+probeSeq+"\n");
			}
		}
		
		writer.close();
		vcf.close();
	}
	
	private boolean rejectProbe(String probe) {
		LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 15, 12);
		RepeatFilter rf=new RepeatFilter(0.07, true, true);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe) || rf.rejectSequence(probe);
	}

	
	
	private Collection<Annotation> getProbes(Gene gene, VCFFileReader vcf, String sp1, String sp2, int readLength) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		CloseableIterator<VariantContext> iter=vcf.query(gene.getReferenceName().replace("chr", ""), gene.getReferenceStartPosition(), gene.getReferenceEndPosition());
		
		//System.err.println(gene.toBED());
		
		while(iter.hasNext()){
			VariantContext variant=iter.next();
			boolean isHet=isHetero(variant, sp1, sp2);
			if(isHet){
				SingleInterval region=new SingleInterval(gene.getReferenceName(), variant.getStart(), variant.getEnd());
				region.setOrientation(gene.getOrientation());
				//System.err.println(region.toUCSC()+" "+gene.toUCSC()+" "+gene.overlaps(region) +" "+region.overlaps(gene)+" "+gene.getNumberOfBlocks());
				if(gene.overlapsExon(region)){
					
					//generate probe
					int relativePos=gene.getRelativePositionFrom5PrimeOfFeature(variant.getStart());
					
					//System.err.println(relativePos);
					int start=Math.max(0, relativePos-(readLength/2));
					int end=start+readLength;
					Annotation window=gene.trimByRelativePositions(start, end);
					rtrn.add(window);
					
					/*String chr=gene.getReferenceName();
					int start=variant.getStart()-(readLength/2);
					int end=variant.getStart()+(readLength/2);
					SingleInterval region=new SingleInterval(chr, start, end, Strand.POSITIVE);
					rtrn.add(region);*/
				}
			}
		}
		
		iter.close();
		
		return rtrn;
	}
	
	private boolean isHetero(VariantContext variant, String sp1, String sp2) {
		return !variant.getGenotype(sp1).sameGenotype(variant.getGenotype(sp2));
	}

	public static void main(String[] args) throws IOException{
		if(args.length>6){
			VCFFileReader vcf=new VCFFileReader(new File(args[0]));
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[1]);
			Map<String, Sequence> sequenceByChr=FastaFileIOImpl.readFromFileByName(args[2]);
			String species1=args[3];
			String species2=args[4];
			String save=args[5];
			int readLength=new Integer(args[6]);
			new GenerateSNPProbes(vcf, genes, sequenceByChr, species1, species2, save, readLength);
		}
		else{System.err.println(usage);}
		
		
		
	}
	
	static String usage=" args[0]=VCF file \n args[1]=genes \n args[2]=genome fasta file \n args[3]=species name 1 \n args[4]=species name 2 \n args[5]=save \n args[6]=probe length";

	
	
}

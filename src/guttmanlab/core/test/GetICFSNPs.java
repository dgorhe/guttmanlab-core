package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.probegeneration.GenerateSNPProbes;
import guttmanlab.core.probegeneration.LowComplexityFilter;
import guttmanlab.core.probegeneration.PolyBaseFilter;
import guttmanlab.core.probegeneration.RepeatFilter;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class GetICFSNPs {

	public GetICFSNPs(Collection<SingleInterval> snps, Map<String, Sequence> sequenceByChr, String save, int probeLength) throws IOException {
			
		Collection<SingleInterval> probes=getProbes(snps, probeLength);
		
		FileWriter writer=new FileWriter(save);
		for(SingleInterval probe: probes){
			Sequence probeSeq=sequenceByChr.get(probe.getReferenceName()).getSubsequence(probe);
			boolean rejectProbe=rejectProbe(probeSeq.getSequenceBases());
			if(!rejectProbe){
				writer.write(probe.getReferenceName()+"\t"+probe.getReferenceStartPosition()+"\t"+probe.getReferenceEndPosition()+"\t"+probeSeq+"\n");
			}
		}
		
		writer.close();
	}
	
	private Collection<SingleInterval> getProbes(Collection<SingleInterval> snps, int readLength) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval snp: snps){
			
				String chr=snp.getReferenceName();
				int start=snp.getReferenceStartPosition()-(readLength/2);
				int end=snp.getReferenceStartPosition()+(readLength/2);
				SingleInterval region=new SingleInterval(chr, start, end, Strand.POSITIVE);
				rtrn.add(region);
			
			
		}
		
		
		return rtrn;
	}

	private boolean rejectProbe(String probe) {
		LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 15, 12);
		RepeatFilter rf=new RepeatFilter(0.07, true, true);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe) || rf.rejectSequence(probe);
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			//VCFFileReader vcf=new VCFFileReader(new File(args[0]));
			Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[0]);
			Map<String, Sequence> sequenceByChr=FastaFileIOImpl.readFromFileByName(args[1]);
			String save=args[2];
			int readLength=new Integer(args[3]);
			new GetICFSNPs(regions, sequenceByChr, save, readLength);
		}
		else{System.err.println(usage);}
		
		
		
	}
	
	static String usage=" args[0]=snp list \n args[1]=genome fasta file \n args[2]=save \n args[3]=probe length";

	
	
}

package guttmanlab.core.test;

import java.util.Collection;
import guttmanlab.core.primer3.PcrPrimerDesigner;
import guttmanlab.core.primer3.PrimerPair;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class DesignCloningPrimer {

	int fragmentLength=1000;
	String primer3core="/central/groups/guttman/software/arraydesign/primer3_core";
	
	
	public DesignCloningPrimer(Sequence geneSeq) throws Exception{
		Collection<PrimerPair> primers=PcrPrimerDesigner.designTilingPrimers(geneSeq, primer3core);
		
		for(PrimerPair p: primers){
			System.out.println(p.getLeftPrimer()+" "+p.getRightPrimer()+" "+p.getLeftPrimerPosition()+" "+ p.getRightPrimerPosition()+" "+p.getProductSize());
		}
		
		
	}

	
	public static void main(String[] args) throws Exception{
		Collection<Sequence> seqs=FastaFileIOImpl.readFromFile(args[0]);
		
		for(Sequence seq: seqs){
			System.err.println(seq.getSequenceBases().length());
			new DesignCloningPrimer(seq);
		}
	}
	
}

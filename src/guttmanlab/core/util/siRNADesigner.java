package guttmanlab.core.util;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.probegeneration.PolyBaseFilter;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class siRNADesigner {

	static int length=21;
	static String common="CCTGTCTC";
	static double minGC=0.3;
	static double maxGC=0.5;
	/*
	 * "Template oligonucleotides (DNA) are 29 nt in length. The 8 nucleotides at the 3' end of
		each oligonucleotide must be 5'-CCTGTCTC-3'. See section “siRNA design” on page 8
		for complete instructions on designing template oligonucleotides."
	 */
	
	
	public static Collection<String> designPossibleSiRNAs(Sequence cDNA, FileWriter writer) throws IOException {
		PolyBaseFilter filter=new PolyBaseFilter("ACGT", 4,4);
		//FileWriter writer=new FileWriter(save);
		Collection<String> rtrn=new TreeSet<String>();
		
		char[] chars=cDNA.getSequenceBases().toCharArray();
		for(int i=0; i<cDNA.getLength()-1; i++) {
			String dinuc=chars[i]+""+chars[i+1];
			if(dinuc.equals("AA")) {
				String seq=get(cDNA, i, length);
				if(seq!=null) {
					double percentGC=Sequence.computeGCContent(seq);
					if(filter.rejectSequence(seq)) {}
					else {
						String sense=getSense(seq);
						String antisense=getAntisense(seq);
						if(percentGC>=minGC && percentGC<=maxGC) {
							writer.write(cDNA.getName()+"\t"+i+"\t"+seq+"\t"+percentGC+"\t"+sense+"\t"+antisense+"\n");
						}
						
						//System.err.println(seq+common);
						rtrn.add(seq);
					}
				}
			}
		}
		//writer.close();
		return rtrn;
	}
	
	

	private static String getAntisense(String target) {
		String antisense="AA"+Sequence.reverseComplement(target.substring(2, target.length()))+common;
		return antisense;
	}

	private static String getSense(String target) {
		return target+common;
	}
	private static String get(Sequence cDNA, int i, int length2) {
		if(i+length2<cDNA.getLength()) {
			return cDNA.getSubSequence(i, i+length2);
		}
		return null;
	}
	
	public static void main(String[] args) throws IOException{
		Collection<Sequence> sequences=FastaFileIOImpl.readFromFile(args[0]);
		String save=args[1];
		FileWriter writer=new FileWriter(save);
		
		for(Sequence seq: sequences) {
			System.err.println(seq.getName());
			designPossibleSiRNAs(seq, writer);
		}
		
		writer.close();
		/*String target="AACGATTGACAGCGGATTGCC";
		
		String sense=target+common;
		String antisense="AA"+Sequence.reverseComplement(target.substring(2, target.length()))+common;
		System.err.println(sense+"\t"+antisense);*/
		
	}
	
	
}

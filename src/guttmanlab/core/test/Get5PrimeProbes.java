package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.probegeneration.LowComplexityFilter;
import guttmanlab.core.probegeneration.PolyBaseFilter;
import guttmanlab.core.probegeneration.RepeatFilter;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class Get5PrimeProbes {

	String overhang="CAAGTCA";
	int probeLength=50;
	//String fastaFile="/groups/guttman/genomes/mm10/mm10withchr.fa";
	
	public Get5PrimeProbes(Collection<Gene> genes, String fastaFile, String save) throws IOException{
		Map<String, Sequence> chrSequences=FastaFileIOImpl.readFromFileByName(fastaFile);
		
		//Get probes
		//Add 5' overhang
		Map<String, String> probes=getProbes(genes, chrSequences);
		
		write(save+".probes", probes);
		
		
		
	}

	private void write(String save, Map<String, String> probes) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String probe: probes.keySet()){
			writer.write(probe+"\t"+probes.get(probe)+"\n");
		}
		
		writer.close();
	}

	private Collection<String> getNames(SingleInterval region, Iterator<Gene> iterator) {
		Collection<String> rtrn=new TreeSet<String>();
		while(iterator.hasNext()){
			Gene gene=iterator.next();
			if(gene.overlaps(region)){rtrn.add(gene.getName());}
		}
		return rtrn;
	}

	private Collection<SingleInterval> getTopN(Map<String, Integer> chrCountsPos, Map<String, Integer> chrCountsNeg, int n) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		Map<Integer, Collection<SingleInterval>> temp=new TreeMap<Integer, Collection<SingleInterval>>();
		
		Map<String, Integer> merge=new TreeMap<String, Integer>();
		merge.putAll(chrCountsPos);
		merge.putAll(chrCountsNeg);
		
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		for(String region: merge.keySet()){
			SingleInterval r=new SingleInterval(region);
			int score;
			if(chrCountsNeg.containsKey(region) && chrCountsPos.containsKey(region)){
				int posScore=chrCountsPos.get(region);
				int negScore=chrCountsNeg.get(region);
				score=-Math.max(posScore, negScore);
				if(posScore>negScore){
					r.setOrientation(Strand.POSITIVE);
				}
				else{r.setOrientation(Strand.NEGATIVE);}
				
			}
			else if(chrCountsNeg.containsKey(region)){
				r.setOrientation(Strand.NEGATIVE);
				score=-chrCountsNeg.get(region);
			}
			else{
				r.setOrientation(Strand.POSITIVE);
				score=-chrCountsPos.get(region);
			}
			scores.put(r, score);
			Collection<SingleInterval> list=new ArrayList<SingleInterval>();
			if(temp.containsKey(score)){list=temp.get(score);}
			list.add(r);
			temp.put(score, list);
		}
		
		
		for(Integer score: temp.keySet()){
			rtrn.addAll(temp.get(score));
			if(rtrn.size()>=n){return rtrn;}
		}
		
		
		return rtrn;
	}

	private Map<String, String> getProbes(Collection<Gene> genes, Map<String, Sequence> chrSequences2) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		int counter=0;
		Collection<SingleInterval> exons=new TreeSet<SingleInterval>();
		for(Gene gene: genes){
			if(gene.hasCodingRegion() && chrSequences2.containsKey(gene.getReferenceName())){exons.add(gene.get5PrimeExon());}
		}
		
		for(SingleInterval exon: exons){
			Sequence regionSeq=chrSequences2.get(exon.getReferenceName()).getSubsequence(exon);
			String probe=getProbe(regionSeq, probeLength);
			if(probe!=null){
				String fullProbe=probe;
				rtrn.put(exon.getName()+"_"+exon.toUCSC(exon.getOrientation()), fullProbe);
			}
			counter++;
			if(counter%1000 ==0){System.err.println(counter+" "+exons.size()+" "+exon.getReferenceName());}
			
		}
		return rtrn;
	}
	
	
	
	
	private Map<String, String> getProbes(Map<String, Integer> chrCounts, Map<String, Sequence> chrSequences2) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(String string: chrCounts.keySet()){
			SingleInterval region=parse(string);
			//System.out.println(region.toBED());
			Sequence regionSeq=chrSequences2.get(region.getReferenceName()).getSubsequence(region);
			
			String probe=getProbe(regionSeq, probeLength);
			if(probe!=null){
				String fullProbe=this.overhang+probe;
				rtrn.put(string, fullProbe);
			}
		}
		return rtrn;
	}

	private SingleInterval parse(String string) {
		SingleInterval rtrn=new SingleInterval(string.split("_")[1]);
		rtrn.setOrientation(Strand.fromString(string.split("_")[2]));
		rtrn.setName(string.split("_")[0]);
		return rtrn;
	}

	private String getProbe(Sequence regionSeq, int probeLength2) {
		for(int i=0; i<regionSeq.getLength()-probeLength2; i++){
			String probe=Sequence.reverseComplement(regionSeq.getSubSequence(i, i+probeLength2));
			if(!rejectProbe(probe)){return probe;}
		}
		
		return null;
	}

	private boolean rejectProbe(String probe) {
		LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 15, 12);
		RepeatFilter rf=new RepeatFilter(0.07, true, true);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe) || rf.rejectSequence(probe);
	}
	
	private Map<String, Integer> getByStrand(Map<String, Integer> chrCountsPos, Map<String, Integer> chrCountsNeg, boolean firstOfPairFlag, boolean readNegativeStrandFlag) {
		boolean pos=false;
		if(firstOfPairFlag){
			pos=!readNegativeStrandFlag;
		}
		else{pos=readNegativeStrandFlag;}
		
		if(pos){return chrCountsPos;}
		return chrCountsNeg;
	}

	/*private void write(String save, Map<String, Integer> chrCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String chr: chrCounts.keySet()){
			SingleInterval region=new SingleInterval(chr);
			int count=chrCounts.get(chr);
			writer.write(region.toBedgraph(count)+"\n");
		}
		
		writer.close();
	}*/
	
	public static void main(String[] args) throws IOException{
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		String save=args[1];
		String fastaFile=args[2];
		new Get5PrimeProbes(genes, fastaFile, save);
	}
	
}

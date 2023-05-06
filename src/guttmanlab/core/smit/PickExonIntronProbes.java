package guttmanlab.core.smit;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.probegeneration.LowComplexityFilter;
import guttmanlab.core.probegeneration.PolyBaseFilter;
import guttmanlab.core.probegeneration.RepeatFilter;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class PickExonIntronProbes {

	int probeLength;
	int numberProbes=5;
	
	public PickExonIntronProbes(Collection<Gene> genes, ArrayList<String> geneOrder, Map<String, Sequence> sequenceByChr, String save, int totalNumProbes, int probeLength) throws IOException{
		this.probeLength=probeLength;
		Collection<String> probes=new TreeSet<String>();
		Collection<String> sequences=new TreeSet<String>();
		
		Map<String, Collection<Annotation>> exons=getExonsByGene(genes);
		Map<String, Collection<Annotation>> introns=getIntronsByGene(genes);
		
		Iterator<String> geneIter=geneOrder.iterator();
		while(geneIter.hasNext() && probes.size()<totalNumProbes){
			String gene=geneIter.next();
			System.out.println(gene+" "+probes.size());
			Collection<Annotation> exonList=exons.get(gene);
			Collection<Annotation> intronList=introns.get(gene);
			
			for(Annotation exon: exonList){
				Collection<Sequence> temp=pickProbes(exon, sequenceByChr, numberProbes, probeLength, sequences);
				for(Sequence seq: temp){probes.add(seq.getSequenceBases());}
			}
			
			for(Annotation intron: intronList){
				Collection<Sequence> temp=pickProbes(intron, sequenceByChr, numberProbes, probeLength, sequences);
				for(Sequence seq: temp){probes.add(seq.getSequenceBases());}
			}
			
		}
		
		
		
		
		write(save, probes);
	}
	
	private Map<String, Collection<Annotation>> getIntronsByGene(Collection<Gene> genes) {
		Map<String, Collection<Annotation>> rtrn=new TreeMap<String, Collection<Annotation>>();
		
		for(Gene gene: genes){
			rtrn.put(gene.getName(), gene.getIntrons());
		}
		
		return rtrn;
	}

	private Map<String, Collection<Annotation>> getExonsByGene(Collection<Gene> genes) {
		Map<String, Collection<Annotation>> rtrn=new TreeMap<String, Collection<Annotation>>();
		
		for(Gene gene: genes){
			rtrn.put(gene.getName(), gene.getBlockSet());
		}
		
		return rtrn;
 	}

	private Collection<Sequence> pickProbes(Annotation intron, Map<String, Sequence> sequenceByChr, int numberProbes, int probeLength, Collection<String> sequences) {
		List<Sequence> rtrn=new ArrayList<Sequence>();
		Sequence geneSeq=sequenceByChr.get(intron.getReferenceName()).getSubsequence(intron);
		
		//enumerate all probes
		for(int i=0; i<geneSeq.getSequenceBases().length(); i+=probeLength){
			int start=i;
			int end=i+probeLength;
			String seq=geneSeq.getSubSequence(start, end);
			//seq=Sequence.reverseComplement(seq);
			
			boolean rejectProbe=rejectProbe(seq, sequences);
			if(!rejectProbe){
				Sequence probe=new Sequence(intron.getName()+"_"+i, seq);
				sequences.add(seq);
				//System.out.println(intron.getReferenceName()+"\t"+(intron.getReferenceStartPosition()+start)+"\t"+(intron.getReferenceStartPosition()+end));
				rtrn.add(probe);
			}
		}
		
		//TODO Get numberOfProbes
		rtrn=get(rtrn, numberProbes);
		
		return rtrn;
	}

	private List<Sequence> get(List<Sequence> rtrn, int numberProbes2) {
		List<Sequence> list=new ArrayList<Sequence>();
		if(rtrn.size()<=numberProbes2){return rtrn;}
		for(int i=0; i<numberProbes2; i++){
			int index=new Double(Math.random()*(i*((double)rtrn.size()/(double)numberProbes2))).intValue();
			list.add(rtrn.remove(index));
		}
		return list;
	}

	private boolean rejectProbe(String probe, Collection<String> sequences) {
		LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 15, 12);
		RepeatFilter rf=new RepeatFilter(0.07, true, true);
		boolean seqLength=probe.length()<probeLength;
		boolean hasSeq=sequences.contains(probe);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe) || rf.rejectSequence(probe) || seqLength || hasSeq;
	}

	private void write(String save, Collection<String> probes) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String probe: probes){
			writer.write(probe+"\n");
		}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		if(args.length>5){
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
			Map<String, Sequence> sequenceByChr=FastaFileIOImpl.readFromFileByName(args[1]);
			ArrayList<String> geneOrder=parseList(args[2]);
			String save=args[3];
			int totalNumProbes=new Integer(args[4]);
			int probeLength=new Integer(args[5]);
			new PickExonIntronProbes(genes, geneOrder, sequenceByChr, save, totalNumProbes, probeLength);
		}
		else{System.err.println(usage);}
	}
	
	private static ArrayList<String> parseList(String file) throws IOException {
		ArrayList<String> rtrn=new ArrayList<String>();
		List<String> lines=BEDFileIO.loadLines(file);
		for(String line: lines){
			String geneName=line.split("\t")[0];
			rtrn.add(geneName);
		}
		return rtrn;
	}

	static String usage=" args[0]=gene BED file \n args[1]=genome fasta file \n args[2]=genes in order \n args[3]=save \n args[4]=total number of probes \n args[5]=probe length";
}

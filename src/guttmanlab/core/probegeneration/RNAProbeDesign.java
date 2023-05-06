package guttmanlab.core.probegeneration;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.ParseGTF;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;



public class RNAProbeDesign {
	
	private static String genomeFasta="/central/groups/guttman/mguttman/ProbeGeneration/mm10MaskedGenome/";
	private static String geneFile="/central/groups/guttman/mguttman/ProbeGeneration/gencode.vM25.annotation.gtf";

	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length>3) {
			
			String geneName=args[0];
			String save=args[1];
			int probeLength=Integer.parseInt(args[2]);
			boolean includeIntrons=Boolean.parseBoolean(args[3]);
			
			
			Collection<Gene> genes=parseGenes(geneFile, geneName);
			
			Map<String, Sequence> genomeSeq=FastaFileIOImpl.readFromFilesByName(new File(genomeFasta).listFiles(), genes);
			
			
			FileWriter writer=new FileWriter(save);
			String sequenceOutput=save+".raw";
			FileWriter writerRaw=new FileWriter(sequenceOutput);
			
			
			Collection<Sequence> sequences=new ArrayList<Sequence>();
			Map<SingleInterval, Sequence> intronSequences=new TreeMap<SingleInterval, Sequence>();
			
			for(Gene g: genes) {
				System.out.println(g.toBED());
				System.err.println(g.getName());
				
				Sequence geneSeq=g.getSequence(genomeSeq.get(g.getReferenceName()));
			
				if(includeIntrons) {
					for(Annotation intron: g.getIntrons()) {
						System.out.println(intron.toBED());
						Sequence intronSeq=new Gene(intron).getSequence(genomeSeq.get(g.getReferenceName()));
						intronSequences.put(intron.getSingleInterval(), intronSeq);
					}
				}
				
				geneSeq.setName(g.getName());
				sequences.add(geneSeq);
			}
			
			
			Map<String, Collection<SingleInterval>> map=new TreeMap<String, Collection<SingleInterval>>();
			
			for(Sequence seq: sequences) {
				System.err.println(seq.getName());
				Map<Integer, String> kmers=seq.enumerateKmer(probeLength);
				for(Integer pos: kmers.keySet()) {
					String kmer= kmers.get(pos);
					kmer=Sequence.reverseComplement(kmer);
					boolean isGood=!rejectProbe(kmer);
					if(isGood) {
						Collection<SingleInterval> list=new TreeSet<SingleInterval>();
						SingleInterval interval=new SingleInterval(seq.getName(), pos, pos+probeLength);
						if(map.containsKey(kmer)) {
							list=map.get(kmer);
						}
						list.add(interval);
						map.put(kmer, list);
					}
				}
			}
			
			if(includeIntrons) {
				for(SingleInterval intron: intronSequences.keySet()) {
					String name=intron.getName();
					System.err.println(name);
					Sequence seq=intronSequences.get(intron);
					Map<Integer, String> kmers=seq.enumerateKmer(probeLength); 
					for(Integer pos: kmers.keySet()) {
						String kmer= kmers.get(pos);
						kmer=Sequence.reverseComplement(kmer);
						boolean isGood=!rejectProbe(kmer);
						if(isGood) {
							Collection<SingleInterval> list=new TreeSet<SingleInterval>();
							SingleInterval interval=new SingleInterval(name, pos, pos+probeLength);
							if(map.containsKey(kmer)) {
								list=map.get(kmer);
							}
							list.add(interval);
							//System.err.println(kmer+" "+interval.toUCSC());
							map.put(kmer, list);
						}
					}
				}
			}
			
			
			for(String kmer: map.keySet()) {
				writerRaw.write(">"+kmer+"\n"+kmer+"\n");
			}
			writerRaw.close();
			
			
			Collection<String> seqsToFilter1=alignToRepeats(sequenceOutput);
			Collection<String> seqsToFilter2=alignToGenome(sequenceOutput); //TODO need to get the sequence of the read, not alignment
			Collection<String> seqsToFilter3=alignToGenomeLocal(sequenceOutput);
			Collection<String> seqsToFilter4=alignToRepeatsLocal(sequenceOutput);
			
			Collection<String> toFilter=new TreeSet<String>();
			toFilter.addAll(seqsToFilter1);
			toFilter.addAll(seqsToFilter2);
			toFilter.addAll(seqsToFilter3);
			toFilter.addAll(seqsToFilter4);
			
			System.err.println(seqsToFilter1.size()+" "+seqsToFilter2.size()+" "+seqsToFilter3.size()+" "+seqsToFilter4.size()+" "+toFilter.size());
			
			for(String kmer: map.keySet()) {
				if(!toFilter.contains(kmer)) {
					writer.write(kmer);
					Collection<SingleInterval> list=map.get(kmer);
					String gene=getGene(list);
					String type=getType(list);
					writer.write("\t"+gene+"\t"+type);
					
					for(SingleInterval pos: list) {writer.write("\t"+pos.toUCSC());}
					writer.write("\n");
				}
				//else {System.out.println(kmer);}
			}
			
			writer.close();
			
		}else {System.err.println(usage);}
		
		
		
		///central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex
		//bowtie2 -r --un probes.noRepeat.seq -x /central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/allRepeats.mouse -U probes.seq -S repeatAligned.sam
		//bowtie2 -r -k 5 -x /groups/guttman/genomes/mus_musculus/GRCm38/bowtie2/GRCm38.primary_assembly.genome -U /groups/guttman/mguttman/pou5f1.probes.raw -S /groups/guttman/mguttman/pou5f1.probes.raw.sam
		//TODO Also consider a local alignment search
		//CountSAM to filter multi aligners
	}
	
	
	
	private static Collection<Gene> parseGenes(String geneFile, String geneName) throws IOException {
		String save=geneFile+"."+geneName+".gtf";
		ParseGTF.writeGTFByGene(new File(geneFile), save, geneName);
		
		Collection<Gene> genes=ParseGTF.getHighQualityProteinCodingGenes(new File(save));
		if(genes.isEmpty()) {
			genes=ParseGTF.getHighQualityGenes(new File(save));
		}
		if(genes.isEmpty()) {
			genes=ParseGTF.getAllGenes(new File(save));
		}
		
		return genes;
	}



	private static Collection<String> alignToGenome(String sequenceOutput) throws IOException, InterruptedException {
		String output=sequenceOutput+".genomeAligned.sam";
		
		//bowtie2 -r -k 5 -x /groups/guttman/genomes/mus_musculus/GRCm38/bowtie2/GRCm38.primary_assembly.genome -U /groups/guttman/mguttman/pou5f1.probes.raw -S /groups/guttman/mguttman/pou5f1.probes.raw.sam
		String cmd="bowtie2 --no-unal --no-head -f -k 5 -x /groups/guttman/genomes/mus_musculus/GRCm38/bowtie2/GRCm38.primary_assembly.genome -U "+sequenceOutput +" -S "+output;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		System.err.println(p.exitValue());
		
		
		Collection<String> rtrn=new TreeSet<String>();
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(output)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens=nextLine.split("\t");
			if(tokens.length>8) {
				String key=tokens[0];
				int count=0;
				if(counts.containsKey(key)) {
					count=counts.get(key);
				}
				count++;
				if(count>1) {
					//System.err.println(key+" "+count);
					rtrn.add(key);
				}
				counts.put(key, count);
			}
		}
		reader.close();
		return rtrn;
	}
	
	private static Collection<String> alignToGenomeLocal(String sequenceOutput) throws IOException, InterruptedException {
		String output=sequenceOutput+".genomeAligned.sam";
		
		//bowtie2 -r -k 5 -x /groups/guttman/genomes/mus_musculus/GRCm38/bowtie2/GRCm38.primary_assembly.genome -U /groups/guttman/mguttman/pou5f1.probes.raw -S /groups/guttman/mguttman/pou5f1.probes.raw.sam
		String cmd="bowtie2 --local --no-unal --no-head -f -k 5 -x /groups/guttman/genomes/mus_musculus/GRCm38/bowtie2/GRCm38.primary_assembly.genome -U "+sequenceOutput +" -S "+output;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		System.err.println(p.exitValue());
		
		
		Collection<String> rtrn=new TreeSet<String>();
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(output)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens=nextLine.split("\t");
			if(tokens.length>8) {
				String key=tokens[0];
				int count=0;
				if(counts.containsKey(key)) {
					count=counts.get(key);
				}
				count++;
				if(count>1) {
					//System.err.println(key+" "+count);
					rtrn.add(key);
				}
				counts.put(key, count);
			}
		}
		reader.close();
		return rtrn;
	}



	private static Collection<String> alignToRepeats(String sequenceOutput) throws IOException, InterruptedException {
		//Process p=Runtime.getRuntime().exec("module add bowtie/2.3.4.1");
		//p.waitFor();
		String output=sequenceOutput+".repeatAligned.sam";
		//bowtie2 -r --un probes.noRepeat.seq -x /central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/allRepeats.mouse -U probes.seq -S repeatAligned.sam
		String cmd="bowtie2 -f  --no-unal --no-head -x /central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/allRepeats.mouse -U "+sequenceOutput+" -S "+output;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		//get sequences that align to repeats
		Collection<String> rtrn=new TreeSet<String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(output)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens=nextLine.split("\t");
			if(tokens.length>8) {
				String key=tokens[0];
				System.err.println(key);
				rtrn.add(key);
			}
		}
		reader.close();
		return rtrn;
	}
	
	private static Collection<String> alignToRepeatsLocal(String sequenceOutput) throws IOException, InterruptedException {
		//Process p=Runtime.getRuntime().exec("module add bowtie/2.3.4.1");
		//p.waitFor();
		
		String output=sequenceOutput+".repeatAligned.sam";
		//bowtie2 -r --un probes.noRepeat.seq -x /central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/allRepeats.mouse -U probes.seq -S repeatAligned.sam
		String cmd="bowtie2 -f  --no-unal --no-head --local -x /central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex/allRepeats.mouse -U "+sequenceOutput+" -S "+output;
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		//get sequences that align to repeats
		Collection<String> rtrn=new TreeSet<String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(output)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens=nextLine.split("\t");
			if(tokens.length>8) {
				String key=tokens[0];
				System.err.println(key);
				rtrn.add(key);
			}
		}
		reader.close();
		return rtrn;
	}



	private static String getType(Collection<SingleInterval> list) {
		Collection<String> set=new TreeSet<String>();
		
		for(SingleInterval r: list) {
			if(r.getReferenceName().contains("intron")) {set.add("intron");}
			else {set.add("exon");}
		}
		
		if(set.size()==1) {return set.iterator().next();}
		return "amb";
	}


	private static String getGene(Collection<SingleInterval> list) {
		Collection<String> set=new TreeSet<String>();
		
		for(SingleInterval r: list) {
			String[] tokens=r.getReferenceName().split("_");
			set.add(tokens[0]);
		}
		
		if(set.size()==1) {return set.iterator().next();}
		return "amb";
	}



	private static boolean rejectProbe(String probe) {
		if(probe.contains("N")) {return true;}
		LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 5, 5); //8,6
		PolyBaseFilter pbf2=new PolyBaseFilter("ACGT", 15, 12);
		PolyBaseFilter pbf3=new PolyBaseFilter("ACGT", 8, 6);
		//RepeatFilter rf=new RepeatFilter(0.07, true, true);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe) || pbf2.rejectSequence(probe) || pbf3.rejectSequence(probe);
	}
	
	static String usage="module add bowtie/2.3.4.1 \n args[0]=gene name \n args[1]=save \n args[2]=probe length \n args[3]=include introns (true/false)";
}




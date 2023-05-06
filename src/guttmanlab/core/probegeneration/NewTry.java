package guttmanlab.core.probegeneration;

	import java.io.File;
	import java.io.FileWriter;
	import java.io.IOException;
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

public class NewTry {


		public static void main(String[] args) throws IOException {
			if(args.length>3) {
			Collection<Gene> genes=ParseGTF.getHighQualityGenes(new File(args[0]));
			Map<String, Sequence> genomeSeq=FastaFileIOImpl.readFromFilesByName(new File(args[1]).listFiles());
			FileWriter writer=new FileWriter(args[2]);
			
			FileWriter writerRaw=new FileWriter(args[2]+".raw");
			
			int probeLength=Integer.parseInt(args[3]);
			
			
			
			Collection<Sequence> sequences=new ArrayList<Sequence>();
			Map<SingleInterval, Sequence> intronSequences=new TreeMap<SingleInterval, Sequence>();
			
			int counter=0;
			for(Gene g: genes) {
				String name=g.getName()+":"+counter;
				g.setName(name);
				System.out.println(g.toBED());
				System.err.println(g.getName());
				
				Sequence geneSeq=g.getSequence(genomeSeq.get(g.getReferenceName()));
			
				for(Annotation intron: g.getIntrons()) {
					System.out.println(intron.toBED());
					Sequence intronSeq=new Gene(intron).getSequence(genomeSeq.get(g.getReferenceName()));
					intronSequences.put(intron.getSingleInterval(), intronSeq);
				}
				
				geneSeq.setName(name);
				sequences.add(geneSeq);
				counter++;
				
				
				
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
			
			for(SingleInterval intron: intronSequences.keySet()) {
				String name=intron.getName()+"_"+intron.toUCSC();
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
			
			
			for(String kmer: map.keySet()) {
				writer.write(kmer);
				writerRaw.write(kmer+"\n");
				Collection<SingleInterval> list=map.get(kmer);
				for(SingleInterval pos: list) {writer.write("\t"+pos.toUCSC());}
				writer.write("\n");
			}
			
			writer.close();
			writerRaw.close();
			}else {System.err.println(usage);}
			
			
			
			///central/groups/guttman/mguttman/bowtieIndex/mouse/bowtieIndex
			//bowtie2 -r -k 5 -x /groups/guttman/genomes/mus_musculus/GRCm38/bowtie2/GRCm38.primary_assembly.genome -U /groups/guttman/mguttman/pou5f1.probes.raw -S /groups/guttman/mguttman/pou5f1.probes.raw.sam
		}
		
		
		
		private static boolean rejectProbe(String probe) {
			if(probe.contains("N")) {return true;}
			LowComplexityFilter lcf=new LowComplexityFilter();
			PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 8, 6);
			//RepeatFilter rf=new RepeatFilter(0.07, true, true);
			return lcf.rejectSequence(probe) || pbf.rejectSequence(probe);
		}
		static String usage=" args[0]=GTF \n args[1]=genome sequence folder \n args[2]=save \n args[3]=probe length";
	}




	


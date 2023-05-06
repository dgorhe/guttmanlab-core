package guttmanlab.core.primer3;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class DesignRTPrimers {
	//TODO Introns and exons
	//TODO add repeat mispriming library
	//TODO Must end with C
	
	static String random=Math.random()+""+System.currentTimeMillis();
	//static String primer3core="/Users/mguttman/primer3/src/primer3_core";
	static String primer3core="/central/groups/guttman/software/arraydesign/primer3_core";
	static String overhang="TGACTTG";
	static int layoutDistance=100;
	static int primerSize=25;
	static int[] primerSizes= {25,24,23,22,21,20};
	
	/*public static Map<String, Map<Integer, String>> design(Collection<Sequence> geneSeqs) throws IOException, InterruptedException {
		Map<String, Map<Integer, String>> rtrn=new TreeMap<String, Map<Integer, String>>();
		for(Sequence geneSeq: geneSeqs) {
			System.err.println(geneSeq.getName());
			Sequence rc=geneSeq.reverseComplement();
			Map<Integer, String> possiblePrimers=rc.enumerateKmer(primerSize);
			possiblePrimers=endsWithC(possiblePrimers);
			Map<Integer, String> passFilter=filter(possiblePrimers, layoutDistance);
			rtrn.put(geneSeq.getName(), passFilter);
		}
		return rtrn;
	}*/
	
	
	public static Map<Integer, String> design(Sequence geneSeq) throws IOException, InterruptedException {
		Sequence rc=geneSeq.reverseComplement();
		Map<Integer, List<String>> possiblePrimers=rc.enumerateKmer(primerSizes);
		possiblePrimers=endsWithC(possiblePrimers);
		Map<Integer, String> passFilter=filter(possiblePrimers, layoutDistance);
		System.err.println(possiblePrimers.size()+" "+passFilter.size());
		return passFilter;
	}
	
	
	private static Map<Integer, List<String>> endsWithC(Map<Integer, List<String>> possiblePrimers) {
		Map<Integer, List<String>> rtrn=new TreeMap<Integer, List<String>>();
		
		for(Integer pos: possiblePrimers.keySet()) {
			List<String> primers=possiblePrimers.get(pos);
			List<String> filtered=new ArrayList<String>();
			for(String primer: primers) {
				if(endsWithC(primer)) {filtered.add(primer);}
			}
			if(!filtered.isEmpty()) {rtrn.put(pos, filtered);}
		}
		
		return rtrn;
	}


	private static boolean endsWithC(String primer) {
		if(primer.endsWith("C") || primer.endsWith("c")) {return true;}
		return false;
	}

	
	private static Map<Integer, String> filter(Map<Integer, List<String>> possiblePrimers, int layoutDistance) throws IOException, InterruptedException {
		Map<Integer, String> rtrn=new TreeMap<Integer, String>();
		
		int currentPos=-1;
		for(Integer pos: possiblePrimers.keySet()) {
			if(checkPos(pos, currentPos, layoutDistance)) {
				String probe=getProbe(possiblePrimers.get(pos));
				if(probe!=null) {
				//if(!rejectProbe(possiblePrimers.get(pos))) {
					rtrn.put(pos, probe);
					currentPos=pos;
				}
			}
		}
		
		System.err.println(possiblePrimers.size()+" "+rtrn.size());
		return rtrn;
	}
	

	/*private static Map<Integer, String> filter(Map<Integer, String> possiblePrimers, int layoutDistance) throws IOException, InterruptedException {
		Map<Integer, String> rtrn=new TreeMap<Integer, String>();
		
		int currentPos=-1;
		for(Integer pos: possiblePrimers.keySet()) {
			if(checkPos(pos, currentPos, layoutDistance)) {
				if(!rejectProbe(possiblePrimers.get(pos))) {
					rtrn.put(pos, possiblePrimers.get(pos));
					currentPos=pos;
				}
			}
		}
		
		System.err.println(possiblePrimers.size()+" "+rtrn.size());
		return rtrn;
	}*/
	
	private static boolean checkPos(Integer pos, int currentPos, int layoutDistance) {
		if(currentPos<0) {return true;}
		if((pos-currentPos)>layoutDistance) {return true;}
		return false;
	}


	private static Collection<String> filter(ArrayList<String> possiblePrimers) throws IOException, InterruptedException {
		Collection<String> rtrn=new ArrayList<String>();
		for(String primer: possiblePrimers) {
			if(!rejectProbe(primer)) {rtrn.add(primer);}
		}
		System.err.println(possiblePrimers.size()+" "+rtrn.size());
		return rtrn;
	}
	
	private static boolean rejectProbe(String probe) throws IOException, InterruptedException {
		/*LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 5, 5);
		//RepeatFilter rf=new RepeatFilter(0.07, true, true);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe);*/
		
		return !passPrimer3(probe);
	}
	
	
	private static String getProbe(List<String> probes) throws IOException, InterruptedException {
		/*LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 5, 5);
		//RepeatFilter rf=new RepeatFilter(0.07, true, true);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe);*/
		
		for(String probe: probes) {
			if(passPrimer3(probe)) {return probe;}
		}
		return null;
	}


	
	private static boolean passPrimer3(String primer) throws IOException, InterruptedException {
		String file=makeTempFile(primer);
		String primerOut=random+"p3Out.txt";
		String command=primer3core+" --output "+primerOut+" "+file;
		Process p=Runtime.getRuntime().exec(command);
		p.waitFor();
		p.getOutputStream();
		return passFilter(primerOut);
	}
	
	/*public static Collection<String> design(Sequence geneSeq, int k) throws IOException {
		Sequence rc=geneSeq.reverseComplement();
		ArrayList<String> possiblePrimers=rc.enumerateKmer(k);
		return possiblePrimers;
	}*/
	
	private static Collection<String> runPrimer3(Collection<String> list) throws IOException, InterruptedException {
		Collection<String> passFilter=new ArrayList<String>();
		for(String primer: list) {
			String file=makeTempFile(primer);
			String primerOut=random+"p3Out.txt";
			String command="/Users/mguttman/primer3/src/primer3_core --output "+primerOut+" "+file;
			Process p=Runtime.getRuntime().exec(command);
			p.waitFor();
			p.getOutputStream();
			if(passFilter(primerOut)) {
				passFilter.add(primer);
				//System.err.println(primer);
			}
		}
		return passFilter;
	}

	private static boolean passFilter(String primerOut) throws IOException {
		List<String> lines=BEDFileIO.loadLines(primerOut);
		
		for(String line: lines) {
			if(line.startsWith("PRIMER_LEFT_NUM_RETURNED")) {
				int num=Integer.parseInt(line.split("=")[1]);
				if(num>0) {return true;}
				else {return false;}
			}
		}
		return false;
	}

	private static String makeTempFile(String primer) throws IOException {
		String out=random+"temp.txt";
		FileWriter writer=new FileWriter(out);
		
		writer.write("SEQUENCE_ID=example3\n"+"SEQUENCE_PRIMER="+primer+"\n"
				+ "PRIMER_TASK=check_primers\n"
				+ "PRIMER_OPT_SIZE=22\n"
				+ "PRIMER_MIN_SIZE=20\n"
				+ "PRIMER_MAX_SIZE=25\n"
				+ "PRIMER_EXPLAIN_FLAG=1\n"
				+ "PRIMER_GC_CLAMP=1\n"
				+ "SEQUENCE_OVERHANG_LEFT="+overhang+"\n"
				//+ "PRIMER_MISPRIMING_LIBRARY==/central/groups/guttman/mguttman/RTPrimers/cat_rodrep_and_simple.fa\n"
				+ "=");
		
		writer.close();
		return out;
	}

	private static void write(String save, Collection<PrimerPair> primers) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(PrimerPair p: primers) {
			writer.write(p.getLeftPrimer()+"\t"+p.getRightPrimer()+"\t"+p.getLeftPrimerPosition()+"\t"+p.getRightPrimerPosition()+"\t"+p.getLeftPrimerPenalty()+"\t"+p.getRightPrimerPenalty()+"\n");
		}
		
		writer.close();
	}
	
	private static void write(Map<Integer, String> primers) throws IOException {
		for(Integer p: primers.keySet()) {
			//System.err.println(p.getLeftPrimer()+"\t"+p.getRightPrimer()+"\t"+p.getLeftPrimerPosition()+"\t"+p.getRightPrimerPosition()+"\t"+p.getLeftPrimerPenalty()+"\t"+p.getRightPrimerPenalty());
			System.out.println(p+"\t"+primers.get(p)+"\t"+(overhang+primers.get(p)));
		}
	}
	
	private static void write(Map<String, Map<Integer, String>> primers, Map<String, String> nameToChr, String save, String type) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(String gene: primers.keySet()) {
			for(Integer pos: primers.get(gene).keySet()) {
				String name1=gene.split(" ")[0];
				String chr=getChr(nameToChr, name1);
				writer.write(chr+"\t"+name1+"\t"+type+"\t"+pos+"\t"+primers.get(gene).get(pos)+"\n");
			}
		}
		writer.close();
	}
	
	
	private static void write(FileWriter writer, Map<Integer, String> primers, Map<String, String> nameToChr, String type, String gene) throws IOException {
		for(Integer pos: primers.keySet()) {
			String name1=gene.split(" ")[0];
			String chr=getChr(nameToChr, name1);
			String primer=primers.get(pos);
			writer.write(chr+"\t"+name1+"\t"+type+"\t"+pos+"\t"+primer+"\t"+(overhang+primer)+"\n");
		}
		
	}
	
	private static String getChr(Map<String, String> nameToChr, String name1) {
		if(nameToChr.containsKey(name1)) {
			return nameToChr.get(name1);
		}
		return "unknown";
	}


	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length>1) {
			Collection<Sequence> sequences=FastaFileIOImpl.readFromFile(args[0]);
			Map<String, String> nameToChr=BEDFileIO.loadNameToChromosome(args[1]);
			String save=args[2];
			String type=args[3];
			
			FileWriter writer=new FileWriter(save);
			
			int counter=0;
			for(Sequence seq: sequences) {
				double percent=100*((double)counter/(double)sequences.size());
				System.err.println(seq.getName()+"\t"+counter+"\t"+sequences.size()+"\t"+percent+"\t"+seq.getLength());
				Map<Integer, String> primers=design(seq);
				write(writer, primers, nameToChr, type, seq.getName());
				counter++;
			}
			writer.close();
			
			
		}
		else {System.err.println(usage);}
	}

	

	static String usage=" args[0]=gene sequences \n args[1]=bed file \n args[2]=save \n args[3]=type (intron or exon)";
	
}

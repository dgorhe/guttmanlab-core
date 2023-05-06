package guttmanlab.core.orfcloning;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotationcollection.GeneCollection;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class DesignORFPrimers {

	//TODO Add all primer at end
	//TODO Add name for each gene
	//TODO Add percentiles by expression
	//TODO Make sure primers are compatible across sets
	
	int primerLength=25;
	String padding="AAAAA";
	
	public DesignORFPrimers(GeneCollection genes, Map<String, String> geneIDToName, Map<String, String> geneTiers, String fastaDir, Collection<Pair<Sequence>> permanentPrimers, String save) throws IOException{
		Sequence currentSeq=null;
		String currentChr=null;
		
		Map<String, Pair<Sequence>> possiblePrimers=getPossiblePrimers(permanentPrimers);
		
		
		Map<Gene, Pair<String>> genePrimers=new TreeMap<Gene, Pair<String>>();
		
		for(Gene gene: genes){
			if(gene.hasCodingRegion()){
				Annotation reverse=get5Primer(gene);
				Annotation forward=get3Primer(gene);
				if(!gene.getReferenceName().equalsIgnoreCase(currentChr)){
					String fasta=fastaDir+"/"+gene.getReferenceName()+".fa";
					Sequence seq=FastaFileIOImpl.readFromFile(fasta).iterator().next();
					currentChr=gene.getReferenceName();
					currentSeq=seq;
				}
				
				String reverseSeq=currentSeq.getSubsequence(reverse).getSequenceBases();
				String forwardSeq=currentSeq.getSubsequence(forward).reverseComplement().getSequenceBases();
				
				Pair<String> primers=new Pair<String>(reverseSeq, forwardSeq);
				
				genePrimers.put(gene, primers);
			}
		}
		
		
		
		write(save, genePrimers, possiblePrimers, geneIDToName, geneTiers);
		
		
		
	}
	
	private Map<String, Pair<Sequence>> getPossiblePrimers(Collection<Pair<Sequence>> permanentPrimers) {
		Map<String, Pair<Sequence>> rtrn=new TreeMap<String, Pair<Sequence>>();
		
		//5' ends need to be mod 3, no stop codons (ideally on either strand), 3' ends can be anything
		for(Pair<Sequence> primer: permanentPrimers){
			String name=primer.getValue1().getName();
			
			Collection<Sequence> primer1Extensions=extend(primer.getValue1());
			Collection<Sequence> primer2Extensions=extend(primer.getValue2());
			
			System.err.println(name+" "+primer1Extensions.size()+" "+primer2Extensions.size());
			
			for(Sequence primer1: primer1Extensions){
				for(Sequence primer2: primer2Extensions){
					if(!hasStopCodon(primer1) && !hasStopCodon(primer2)){
						Pair<Sequence> newPrimer=new Pair<Sequence>(primer1, primer2);
						rtrn.put(name, newPrimer);
					}
				}
			}	
		}
			
		System.err.println(rtrn.size());
		return rtrn;
	}

	private Collection<Sequence> additionalPrimers(Sequence primer, int add) {
		Collection<Sequence> rtrn=new ArrayList<Sequence>();
		rtrn.add(primer);
		for(int i=0; i<add; i++){
			Collection<Sequence> temp=new ArrayList<Sequence>();
			for(Sequence primer1: rtrn){
				temp.addAll(extend(primer1));
			}
			rtrn=temp;
		}
		return rtrn;
	}

	private Collection<Sequence> extend(Sequence primer1) {
		String name=primer1.getName();
		
		Collection<Sequence> rtrn=new ArrayList<Sequence>();
		
		int add=3-primer1.getSequenceBases().length()%3;
		//System.err.println(add);
		
		if(add==3){
			rtrn.add(primer1);
		}
		
		if(add==1){
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"A"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"T"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"C"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"G"));
		}
		
		if(add==2){
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"AA"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"TA"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"CA"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"GA"));
			
			rtrn.add(new Sequence(name, primer1.getSequenceBases()+"AC"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"TC"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"CC"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"GC"));
			
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"AG"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"TG"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"CG"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"GG"));
			
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"AT"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"TT"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"CT"));
			rtrn.add(new Sequence(name,primer1.getSequenceBases()+"GT"));
		}
		
		return rtrn;
	}

	private boolean hasStopCodon(Sequence primer) {
		Collection<String> stopCodons=new ArrayList<String>();
		stopCodons.add("TAG");
		stopCodons.add("TAA");
		stopCodons.add("TGA");
		stopCodons.add("CTA");
		stopCodons.add("TTA");
		stopCodons.add("TCA");
		
		//char[] bases=primer.getSequenceBases().toCharArray();
		for(int i=0; i<primer.getSequenceBases().length(); i+=3){
			//StringBuilder codon=new StringBuilder(bases[i]+bases[i+1]+bases[i+2]);
			String codon=primer.getSequenceBases().substring(i, i+3);
			if(stopCodons.contains(codon.toString())){return true;}
		}
		return false;
	}

	private void write(String save, Map<Gene, Pair<String>> genePrimers, Map<String, Pair<Sequence>> possiblePrimers, Map<String, String> geneIDToName, Map<String, String> tiers) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Iterator<String> primerIter=possiblePrimers.keySet().iterator();
		String name1=primerIter.next();
		String name2=primerIter.next();
		String name3=primerIter.next();
		Pair<Sequence> allPP=possiblePrimers.remove(name1);
		Pair<Sequence> forwardPP=possiblePrimers.remove(name2);
		Pair<Sequence> reversePP=possiblePrimers.remove(name3);
		
		Map<String, Pair<Sequence>> forwardPrimerTiers=getPrimersByTier(possiblePrimers, tiers);
		Map<String, Pair<Sequence>> reversePrimerTiers=getPrimersByTier(possiblePrimers, tiers);
		
		writeRemainingPrimers(save+".remaining", possiblePrimers, allPP, forwardPP, reversePP);
		
		Collection<Pair<String>> primerList=new HashSet<Pair<String>>();
		
		for(Gene gene: genePrimers.keySet()){
			Pair<String> oligos=genePrimers.get(gene);
			if(!primerList.contains(oligos)){
				String geneID=gene.getName();
				String geneName=geneIDToName.get(geneID);
				String tier=tiers.get(geneID);
				if(tier==null){tier="12";}
				if(forwardPrimerTiers.containsKey(tier)){
					Pair<Sequence> forPrimersTier=forwardPrimerTiers.get(tier);
					Pair<Sequence> revPrimersTier=reversePrimerTiers.get(tier);
					
					String fullOligoForward=""+allPP.getValue1()+forPrimersTier.getValue1()+"CACC"+forwardPP.getValue1()+oligos.getValue1()+forwardPP.getValue2().reverseComplement()+forPrimersTier.getValue2().reverseComplement()+allPP.getValue2().reverseComplement();		
					String fullOligoReverse=""+allPP.getValue1()+revPrimersTier.getValue1()+reversePP.getValue1()+oligos.getValue2()+reversePP.getValue2().reverseComplement()+revPrimersTier.getValue2().reverseComplement()+allPP.getValue2().reverseComplement();
					
					String fullOligoForwardRC=new Sequence(fullOligoForward).reverseComplement().getSequenceBases();
					String fullOligoReverseRC=new Sequence(fullOligoReverse).reverseComplement().getSequenceBases();
					
					String plateOligoForward="CACC"+forwardPP.getValue1()+oligos.getValue1();		
					String plateOligoReverse=reversePP.getValue1()+oligos.getValue2();
					
					String plateOligoForward2="CACC"+forwardPP.getValue1()+oligos.getValue1()+forwardPP.getValue2();		
					String plateOligoReverse2=reversePP.getValue1()+oligos.getValue2()+reversePP.getValue2();
					
					writer.write(geneID+"\t"+geneName+"\t"+gene.getCodingRegion().size()+"\tForward\t"+tier+"\t"+oligos.getValue1()+"\t"+plateOligoForward+"\t"+plateOligoForward2+"\t"+(fullOligoForward+padding)+"\t"+(fullOligoForwardRC+padding)+"\n");
					writer.write(geneID+"\t"+geneName+"\t"+gene.getCodingRegion().size()+"\tReverse\t"+tier+"\t"+oligos.getValue2()+"\t"+plateOligoReverse+"\t"+plateOligoReverse2+"\t"+(fullOligoReverse+padding)+"\t"+(fullOligoReverseRC+padding)+"\n");	
				}
				else{System.err.println(geneID+" "+geneName+" "+tier);}
			}
			//else{System.err.println(gene.getName());}
			primerList.add(oligos);
		}
		
		writer.close();
		
		FileWriter writerPrimers=new FileWriter(save+".primerList");
		writerPrimers.write("All\t"+name1+"\t"+allPP.getValue1()+"\t"+allPP.getValue2()+"\n");
		writerPrimers.write("allForward\t"+name2+"\t"+forwardPP.getValue1()+"\t"+forwardPP.getValue2()+"\n");
		writerPrimers.write("allReverse\t"+name3+"\t"+reversePP.getValue1()+"\t"+reversePP.getValue2()+"\n");
		
		for(String tier: forwardPrimerTiers.keySet()){
			writerPrimers.write("Tier "+tier+" Forward\t"+forwardPrimerTiers.get(tier).getValue1().getName()+"\t"+forwardPrimerTiers.get(tier).getValue1()+"\t"+forwardPrimerTiers.get(tier).getValue2()+"\n");
			writerPrimers.write("Tier "+tier+" Reverse\t"+reversePrimerTiers.get(tier).getValue1().getName()+"\t"+reversePrimerTiers.get(tier).getValue1()+"\t"+reversePrimerTiers.get(tier).getValue2()+"\n");
		}
		
		writerPrimers.close();
		
	}

	private void writeRemainingPrimers(String save, Map<String, Pair<Sequence>> possiblePrimers, Pair<Sequence> allPP, Pair<Sequence> forwardPP, Pair<Sequence> reversePP) throws IOException {
		
		FileWriter writer=new FileWriter(save);
		
		writer.write(allPP.getValue1().getName()+"\t"+allPP.getValue1().getSequenceBases()+"\t"+allPP.getValue2().getSequenceBases()+"\n");
		writer.write(forwardPP.getValue1().getName()+"\t"+forwardPP.getValue1().getSequenceBases()+"\t"+forwardPP.getValue2().getSequenceBases()+"\n");
		writer.write(reversePP.getValue1().getName()+"\t"+reversePP.getValue1().getSequenceBases()+"\t"+reversePP.getValue2().getSequenceBases()+"\n");
		
		for(String name: possiblePrimers.keySet()){
			writer.write(name+"\t"+possiblePrimers.get(name).getValue1().getSequenceBases()+"\t"+possiblePrimers.get(name).getValue2()+"\n");
		}
		
		writer.close();
	}

	private Map<String, Pair<Sequence>> getPrimersByTier(Map<String, Pair<Sequence>> possiblePrimers, Map<String, String> tiers) {
		Map<String, Pair<Sequence>> rtrn=new TreeMap<String, Pair<Sequence>>();
		
		Collection<String> allTiers=new TreeSet<String>();
		allTiers.addAll(tiers.values());
		
		Iterator<String> names=possiblePrimers.keySet().iterator();
		
		Collection<String> toRemove=new TreeSet<String>();
		for(String tier: allTiers){
			String name=names.next();
			toRemove.add(name);
			Pair<Sequence> primer=possiblePrimers.get(name);
			rtrn.put(tier, primer);
		}
		
		for(String name: toRemove){
			possiblePrimers.remove(name);
		}
		
		return rtrn;
	}

	private Annotation get3Primer(Gene gene) {
		Gene cds=new Gene(gene.getCodingRegion());
		return cds.get3Prime(3, primerLength);
	}

	private Annotation get5Primer(Gene gene) {
		Gene cds=new Gene(gene.getCodingRegion());
		return cds.get5Prime(3, primerLength);
	}
	
	private static Collection<Pair<Sequence>> parse(String fileName) throws IOException {
		Collection<Pair<Sequence>> rtrn=new ArrayList<Pair<Sequence>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String name=tokens[0];
			Sequence primer1=new Sequence(name, tokens[1]);
			Sequence primer2=new Sequence(name, tokens[2]);
			Pair<Sequence> primers=new Pair<Sequence>(primer1, primer2);
			rtrn.add(primers);
		}
		reader.close();
		return rtrn;
	}
	
	private static Map<String, String> getIDToName(String fileName) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String id=tokens[0];
			String name=tokens[1];
			rtrn.put(id, name);
		}
		reader.close();
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		//Global Mouse Tiers
		/*GeneCollection genes=new GeneCollection("/Users/mguttman/Data/Ref.mm9.bed");
		Map<String, String> geneIDToName=getIDToName("/Users/mguttman/Data/mm9.names");
		Map<String, String> geneTiers=getIDToName("/Users/mguttman/Data/mm9.tiers.txt");
		String fasta="/Users/mguttman/Data/mm9Genome/";
		Collection<Pair<Sequence>> primers=parse("/Users/mguttman/Data/Permanent_Primers_Array_Design.txt");
		String save="/Users/mguttman/Data/Primers/mm9.globalTiers.primers";
		new DesignORFPrimers(genes, geneIDToName, geneTiers, fasta, primers, save);*/
		
		
		//ES Mouse Tiers
		/*GeneCollection genes=new GeneCollection("/Users/mguttman/Data/Ref.mm9.bed");
		Map<String, String> geneIDToName=getIDToName("/Users/mguttman/Data/mm9.names");
		Map<String, String> geneTiers=getIDToName("/Users/mguttman/Data/ESTiers.txt");
		String fasta="/Users/mguttman/Data/mm9Genome/";
		Collection<Pair<Sequence>> primers=parse("/Users/mguttman/Data/mm9.globalTiers.primers.remaining");
		String save="/Users/mguttman/Data/mm9.ESTiers.primers";
		new DesignORFPrimers(genes, geneIDToName, geneTiers, fasta, primers, save);*/
		
		//Human Tiers
		GeneCollection genes=new GeneCollection("/Users/mguttman/Data/Ref.hg19.bed");
		Map<String, String> geneIDToName=getIDToName("/Users/mguttman/Data/hg19.refseqToNames.bed");
		Map<String, String> geneTiers=getIDToName("/Users/mguttman/Data/hg19.tiers.txt");
		String fasta="/Users/mguttman/Data/hg19Genome/";
		Collection<Pair<Sequence>> primers=parse("/Users/mguttman/Data/Primers/mm9.globalTiers.primers.remaining");
		String save="/Users/mguttman/Data/Primers/hg19.globalTiers.primers";
		new DesignORFPrimers(genes, geneIDToName, geneTiers, fasta, primers, save);
		
		
		
		/*if(args.length>5){
			GeneCollection genes=new GeneCollection(args[0]);
			Map<String, String> geneIDToName=getIDToName(args[1]);
			Map<String, String> geneTiers=getIDToName(args[2]);
			
			String fasta=args[3];
			Collection<Pair<Sequence>> primers=parse(args[4]);
			String save=args[5];
			new DesignORFPrimers(genes, geneIDToName, geneTiers, fasta, primers, save);
		}
		else{System.err.println(usage);}*/
	}

	

	static String usage=" args[0]=gene BED file \n args[1]=ID to name \n args[2]=gene tiers \n args[3]=fasta dir \n args[4]=permanent primers \n args[5]=save";
	
}

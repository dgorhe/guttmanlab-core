package guttmanlab.core.annotation.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;

public class ParseGTF {

	
	public static Collection<Gene> getHighQualityGenes(File gtfFile) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
		
		Map<String, Map<String, Collection<SingleInterval>>> map=new TreeMap<String, Map<String, Collection<SingleInterval>>>();
		
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
				String[] tokens=nextLine.split("\t");
				if(tokens[2].equals("exon")) {
					SingleInterval exon=new SingleInterval(tokens[0], Integer.parseInt(tokens[3])-1, Integer.parseInt(tokens[4]));
					exon.setOrientation(Strand.fromString(tokens[6]));
					String geneName=getGeneName(tokens[8]);
					String transcript=getTranscriptID(tokens[8]);
					String type=getTranscriptType(tokens[8]);
					String supportLevel=getTranscriptSupportLevel(tokens[8]);
					String transcriptStatus=getTranscriptStatus(tokens[8]); //transcript_status "PUTATIVE"
					if(transcript!=null) {// && type.equals("protein_coding")
						//System.err.println(supportLevel);
						if(supportLevel!=null && supportLevel.equals("1")) {
							Collection<SingleInterval> set=set(map, exon, geneName, transcript);
						}
					}
				}
			}
		}
		reader.close();
		return getJunctions(map);
	}
	
	public static Collection<Gene> getHighQualityProteinCodingGenes(File gtfFile) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
		
		Map<String, Map<String, Collection<SingleInterval>>> map=new TreeMap<String, Map<String, Collection<SingleInterval>>>();
		
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
				String[] tokens=nextLine.split("\t");
				if(tokens[2].equals("exon")) {
					SingleInterval exon=new SingleInterval(tokens[0], Integer.parseInt(tokens[3])-1, Integer.parseInt(tokens[4]));
					exon.setOrientation(Strand.fromString(tokens[6]));
					String geneName=getGeneName(tokens[8]);
					String transcript=getTranscriptID(tokens[8]);
					String type=getTranscriptType(tokens[8]);
					String supportLevel=getTranscriptSupportLevel(tokens[8]);
					String transcriptStatus=getTranscriptStatus(tokens[8]); //transcript_status "PUTATIVE"
					if(transcript!=null && type.equals("protein_coding")) {// 
						//System.err.println(supportLevel);
						if(supportLevel!=null && supportLevel.equals("1")) {
							Collection<SingleInterval> set=set(map, exon, geneName, transcript);
						}
					}
				}
			}
		}
		reader.close();
		return getJunctions(map);
	}
	
	
	public static Collection<Gene> getAllGenes(File gtfFile) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
		
		Map<String, Map<String, Collection<SingleInterval>>> map=new TreeMap<String, Map<String, Collection<SingleInterval>>>();
		
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
				String[] tokens=nextLine.split("\t");
				if(tokens[2].equals("exon")) {
					SingleInterval exon=new SingleInterval(tokens[0], Integer.parseInt(tokens[3])-1, Integer.parseInt(tokens[4]));
					exon.setOrientation(Strand.fromString(tokens[6]));
					String geneName=getGeneName(tokens[8]);
					String transcript=getTranscriptID(tokens[8]);
					String type=getTranscriptType(tokens[8]);
					String supportLevel=getTranscriptSupportLevel(tokens[8]);
					String transcriptStatus=getTranscriptStatus(tokens[8]); //transcript_status "PUTATIVE"
					if(transcript!=null) {// && type.equals("protein_coding")
						//System.err.println(supportLevel);
						Collection<SingleInterval> set=set(map, exon, geneName, transcript);
					}
				}
			}
		}
		reader.close();
		return getJunctions(map);
	}
	
	
	private static Collection<Gene> getJunctions(Map<String, Map<String, Collection<SingleInterval>>> map) {
		Collection<Gene> set=new TreeSet<Gene>();
		for(String gene: map.keySet()) {
			for(String transcript: map.get(gene).keySet()) {
				Collection<SingleInterval> exons=map.get(gene).get(transcript);
				Gene g=new Gene(exons, gene);
				g.setName(gene+"_"+transcript);
				set.add(g);
			}
		}
		return set;
	}
	
	
	private static Collection<SingleInterval> set(Map<String, Map<String, Collection<SingleInterval>>> map, SingleInterval exon, String geneName, String transcript) {
		if(!map.containsKey(geneName)) {
			map.put(geneName, new TreeMap<String, Collection<SingleInterval>>());
		}
		
		Map<String, Collection<SingleInterval>> rtrn=map.get(geneName);
		if(!rtrn.containsKey(transcript)) {
			rtrn.put(transcript, new TreeSet<SingleInterval>());
		}
		
		Collection<SingleInterval> set=rtrn.get(transcript);
		set.add(exon);
		return set;
	}
	
	private static String getGeneName(String string) {
		return getTag(string, "gene_name");
	}
	
	private static String getTag(String string, String tag) {
		String[] tokens=string.split(";");
		
		for(int i=0; i<tokens.length; i++) {
			String token=tokens[i].trim();
			String key=token.split(" ")[0];
			String val=token.split(" ")[1];
			val=val.replaceAll("\"","");
			if(key.equals(tag)) {return val;}	
		}
		return null;
	}
	
	private static String getTranscriptID(String string) {
		return getTag(string, "transcript_id");
	}
	
	private static String getTranscriptType(String string) {
		return getTag(string, "transcript_type");
	}
	
	
	private static String getTranscriptSupportLevel(String string) {
		return getTag(string, "transcript_support_level");
	}
	
	
	private static String getTranscriptStatus(String string) {
		return getTag(string, "transcript_status");
	}

	public static void writeGTFByGene(File gtfFile, String save, String geneNameToUse) throws NumberFormatException, IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gtfFile)));
		FileWriter writer=new FileWriter(save);
		
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			if(!nextLine.startsWith("#")) {
				String[] tokens=nextLine.split("\t");
				String geneName=getGeneName(tokens[8]);
				if(geneName.equalsIgnoreCase(geneNameToUse)) {writer.write(nextLine+"\n");}
				
			}
		}
		reader.close();
		writer.close();
		
	}
	
	
}

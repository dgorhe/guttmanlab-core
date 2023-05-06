package guttmanlab.core.test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class CompressProteins {
	
	Collection<Character> aminoacids;
	
	public CompressProteins(String fasta, String save, Collection<Character> aminoacids) throws IOException{
		this.aminoacids=aminoacids;
		FileWriter writer=new FileWriter(save);
		Collection<Sequence> proteins=FastaFileIOImpl.readFromFile(fasta);
		//Map<String, Collection<String>> compressedToName=new TreeMap<String, Collection<String>>();
		
		//Map<String, Collection<String>> peptidePatterns=new TreeMap<String, Collection<String>>();
		for(Sequence protein: proteins){
			String name=protein.getName();
			String seq=protein.getSequenceBases();
			
			String totalCompress=compress(protein.getSequenceBases());
			
			writer.write(name);
			
			Collection<String> list=new TreeSet<String>();
			Collection<String> peptides=cleave(seq, "R");
			//System.err.println(seq+"\t"+peptides.size());
			for(String peptide: peptides){
				String compressedPeptide=compress(peptide);
				if(!compressedPeptide.isEmpty()){
				Collection<String> sub=cleave(peptide, "P");
				//System.err.println(peptide+" "+sub.size());
				for(String s: sub){
					String cs=compress(s);
					if(!cs.isEmpty()){
					String temp=compressedPeptide+":"+cs;
					//System.err.println(compressedPeptide+"\t"+name);
					//if(!list.contains(temp)){writer.write("\t"+temp);}
					list.add(temp);
					}
				}
				}
			}
			
			writer.write("\t");
			for(String l: list){writer.write("\t"+totalCompress+":"+l);}
			
			writer.write("\n");
			
			
			
		}
		
		
		writer.close();
		
		/*writer=new FileWriter(save+".counts");
		
		for(String c: compressedToName.keySet()){
			writer.write(c+"\t"+compressedToName.get(c).size());
			Collection<String> names=compressedToName.get(c);
			for(String name: names){
				writer.write("\t"+name);
			}
			writer.write("\n");
		}
		
		writer.close();*/
		
	}

	private Collection<String> cleave(Collection<String> seqs, String pattern) {
		Collection<String> rtrn=new TreeSet<String>();
		for(String seq: seqs){
			rtrn.addAll(cleave(seq, pattern));
		}
		return rtrn;
	}
	
	private Collection<String> cleave(String seq, String pattern) {
		Collection<String> rtrn=new TreeSet<String>();
		String[] tokens=seq.split(pattern);
		
		for(int i=0; i<tokens.length; i++){rtrn.add(tokens[i]+""+pattern);}
		
		return rtrn;
	}

	private Collection<String> compress(Collection<String> peptides) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String p: peptides){
			String c=compress(p);
			rtrn.add(c);
		}
		
		return rtrn;
	}

	private String compress(String seq) {
		Map<Character, Integer> counter=new TreeMap<Character, Integer>();
		
		for(int i=0; i<seq.toCharArray().length; i++){
			char aa=seq.toCharArray()[i];
			int count=0;
			if(counter.containsKey(aa)){count=counter.get(aa);}
			count++;
			counter.put(aa, count);
		}
		
		String rtrn="";
		for(Character c: counter.keySet()){
			if(aminoacids.contains(c)){
				int count=counter.get(c);
				rtrn+=""+c+""+count;
			}
		}
		
		//rtrn+="L"+seq.length();
		
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException{
		Collection<Character> aminoacids=new TreeSet<Character>();
		aminoacids.add('K');
		aminoacids.add('C');
		aminoacids.add('R');
		aminoacids.add('H');
		aminoacids.add('D');
		aminoacids.add('E');
		aminoacids.add('Y');
		
		String in=args[0];
		String out=args[1];
		new CompressProteins(in, out, aminoacids);
	}
	
}

package guttmanlab.core.sars.leader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class ComputedNdS {

	public ComputedNdS(File fasta, String save) throws IOException {
		FastaFileIOImpl f=new FastaFileIOImpl(fasta);
		f.hasNext();
		Sequence reference=f.next();
		String refSeq=Sequence.trim(reference.getSequenceBases());
		
		Map<String, Pair<Integer>> possibleNandS=getPossibleNandS(Sequence.getCodonTable());
		
		Map<Integer, Integer> dN=new TreeMap<Integer, Integer>();
		Map<Integer, Integer> dS=new TreeMap<Integer, Integer>();
		
		int total=0;
		int counter=0;
		while(f.hasNext()) {
			Sequence other=f.next();
			String otherSeq=Sequence.trim(other.getSequenceBases());
			if(otherSeq.length()==refSeq.length()) {
				Map<Integer, Pair<Integer>> counts=countNandS(refSeq, otherSeq);
				add(dN, dS, counts);
			}
			else {counter++;}
			total++;
			if(total%10000==0) {System.err.println(total);}
		}
		System.err.println(counter);
		
		write(save, dN, dS, refSeq, reference.translate(), possibleNandS);
		
	}
	
	
	
	
	private Map<String, Pair<Integer>> getPossibleNandS(Map<String, String> codonTable) {
		Map<String, Pair<Integer>> rtrn=new TreeMap<String, Pair<Integer>>();
		for(String codon: codonTable.keySet()) {
			Pair<Integer> n=getPossibleNandS(codon, codonTable);
			System.err.println(codon+" "+n.getValue1()+" "+n.getValue2());
			rtrn.put(codon, n);
		}
		return rtrn;
	}




	private void write(String save, Map<Integer, Integer> dN, Map<Integer, Integer> dS, String dnaSeq, String protein, Map<String, Pair<Integer>> possibleNandS2) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Integer> allPos=new TreeSet<Integer>();
		allPos.addAll(dN.keySet());
		allPos.addAll(dS.keySet());
		
		char[] dnaChars=dnaSeq.toCharArray();
		char[] chars=protein.toCharArray();
		for(Integer pos: allPos) {
			int n=get(dN, pos);
			int s=get(dS, pos);
			String codon=getCodon(dnaChars, pos);
			Pair<Integer> possibleNandS=possibleNandS2.get(codon);
			writer.write(chars[pos]+"\t"+codon+"\t"+(pos+1)+"\t"+n+"\t"+s+"\t"+possibleNandS.getValue1()+"\t"+possibleNandS.getValue2()+"\n");
		}
		
		writer.close();
	}




	private Pair<Integer> getPossibleNandS(String codon, Map<String, String> codonTable) {
		Pair<Integer> rtrn=new Pair<Integer>();
		
		Collection<String> allPossibleMutants=mutate(codon);
		
		String refAA=codonTable.get(codon);
		
		int s=0;
		int n=0;
		
		for(String other: allPossibleMutants) {
			if(!other.equals(codon)) {
				String otherAA=codonTable.get(other);
				if(otherAA==null) {System.err.println(codon+" "+other);}
				if(otherAA.equals(refAA)) {s++;}
				else {n++;}
			}
		}
		
		rtrn.setValue1(n);
		rtrn.setValue2(s);
		
		return rtrn;
	}







	private Collection<String> mutate(String codon) {
		Collection<String> rtrn=new TreeSet<String>();
		
		char[] chars=codon.toCharArray();
		
		rtrn.add("A"+chars[1]+""+chars[2]);
		rtrn.add("C"+chars[1]+""+chars[2]);
		rtrn.add("T"+chars[1]+""+chars[2]);
		rtrn.add("G"+chars[1]+""+chars[2]);
		
		rtrn.add(chars[0]+"A"+chars[2]);
		rtrn.add(chars[0]+"C"+chars[2]);
		rtrn.add(chars[0]+"T"+chars[2]);
		rtrn.add(chars[0]+"G"+chars[2]);
		
		rtrn.add(chars[0]+""+chars[1]+"A");
		rtrn.add(chars[0]+""+chars[1]+"C");
		rtrn.add(chars[0]+""+chars[1]+"T");
		rtrn.add(chars[0]+""+chars[1]+"G");
		
		return rtrn;
	}




	private Collection<String> swapPosition(char[] chars, int i) {
		Collection<String> rtrn=new TreeSet<String>();
	
		chars[i]='A';
		makeString(chars, rtrn);
		chars[i]='C';
		makeString(chars, rtrn);
		chars[i]='G';
		makeString(chars, rtrn);
		chars[i]='T';
		makeString(chars, rtrn);
		
		return rtrn;
	}




	private void makeString(char[] chars, Collection<String> rtrn) {
		String temp=chars[0]+""+chars[1]+""+chars[2];
		rtrn.add(temp);
	}




	private String getCodon(char[] dnaChars, Integer pos) {
		int start=pos*3;
		return dnaChars[start]+""+dnaChars[start+1]+""+dnaChars[start+2];
	}




	private int get(Map<Integer, Integer> dN, Integer pos) {
		int rtrn=0;
		if(dN.containsKey(pos)) {rtrn=dN.get(pos);}
		return rtrn;
	}




	private Map<Integer, Pair<Integer>> countNandS(String reference, String other) {
		Map<Integer, Pair<Integer>> rtrn=new TreeMap<Integer, Pair<Integer>>();
		char[] chars1=reference.toCharArray();
		char[] chars2=other.toCharArray();
		
		Map<String, String> codonTable=Sequence.getCodonTable();
		boolean hitEnd=false;
		for(int i=0; i<chars1.length-3; i+=3) {
			if(!hitEnd) {
				String codon1=chars1[i]+""+chars1[i+1]+""+chars1[i+2];
				String codon2=chars2[i]+""+chars2[i+1]+""+chars2[i+2];
				String AA1=codonTable.get(codon1);
				String AA2=codonTable.get(codon2);
				
				/*if(AA1==null) {AA1="?";}
				if(AA2==null) {AA2="?";}*/
				
				if(AA1!=null && AA2!=null) {
					if(!AA1.equals(AA2)) {
						rtrn.put(i/3, new Pair<Integer>(1,0));
					}
					else if(!codon1.equalsIgnoreCase(codon2)){
						//System.err.println(codon1+ " "+codon2+" "+AA1+" "+AA2);
						rtrn.put(i/3, new Pair<Integer>(0,1));
					}
				}
				
				if(AA1.equals("*") && AA2.equals("*")) {hitEnd=true;}
			}
			
			
		}
		return rtrn;
	}




	private void add(Map<Integer, Integer> dN, Map<Integer, Integer> dS, Map<Integer, Pair<Integer>> counts) {
		for(Integer pos: counts.keySet()) {
			Pair<Integer> pair=counts.get(pos);
			int nCount=0;
			int sCount=0;
			if(dN.containsKey(pos)) {nCount=dN.get(pos);}
			if(dS.containsKey(pos)) {sCount=dS.get(pos);}
			nCount+=pair.getValue1();
			sCount+=pair.getValue2();
			dS.put(pos, sCount);
			dN.put(pos, nCount);
		}
	}


	public static void main(String[] args) throws IOException {
		File f=new File(args[0]);
		String save=args[1];
		new ComputedNdS(f, save);
	}
	
}

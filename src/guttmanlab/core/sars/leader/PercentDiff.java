package guttmanlab.core.sars.leader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class PercentDiff {

	public PercentDiff(File fasta, String save) throws IOException {
		FastaFileIOImpl f=new FastaFileIOImpl(fasta);
		
		FileWriter writer=new FileWriter(save);
		Sequence reference=null;
		if(f.hasNext()) {
			reference=f.next();
		}
		writer.write(">"+reference.getName()+"\n");
		writer.write(reference.getSequenceBases()+"\n");
		
		
		
		int counter=0;
		while(f.hasNext()) {
			Sequence other=f.next();
			if(!other.getSequenceBases().startsWith("-")) {
				writer.write(">"+other.getName()+"\n");
				writer.write(other.getSequenceBases()+"\n");
			}
			/*
			String trim=Sequence.trim(other.getSequenceBases());
			if(!isEmpty(other.getSequenceBases())) {
				if(!reference.getSequenceBases().equals(other.getSequenceBases())) {
					double percentDiff=difference(reference, other);
					if(percentDiff>0 && percentDiff<0.2) {
						double gapFreq=(double)trim.length()/(double)reference.getLength();
						if(gapFreq>0.75) {
							String diffString=differenceString(reference, other);
							//writer.write(diffString+"\t"+percentDiff+"\t"+trim.length()+"\t"+gapFreq+"\n");
							if(!diffString.contains("N")) {
							writer.write(">"+other.getName()+"\n");
							writer.write(other.getSequenceBases()+"\n");
							}
						}
					}
				}
			}*/
			counter++;
			if(counter%10000==0) {System.err.println(counter);}
		}
	
		
		writer.close();
	}

	private boolean isEmpty(String diff) {
		int other=0;
		char[] chars=diff.toCharArray();
		
		for(int i=0; i<chars.length; i++) {
			if(chars[i]=='-' || chars[i]=='.' || chars[i]=='N') {}
			else {other++;}
		}
		if(other>0) {return false;}
		return true;
	}

	
	private void write(String save, Map<Integer, Map<Character, Integer>> variantCounter, Sequence reference) throws IOException {
		FileWriter writer=new FileWriter(save);
		char[] chars=reference.getSequenceBases().toCharArray();
		for(Integer pos: variantCounter.keySet()) {
			
			writer.write(chars[pos]+"\t"+(pos+1));
			Map<Character, Integer> temp=variantCounter.get(pos);
			for(Character c: temp.keySet()) {writer.write("\t"+c+"="+temp.get(c));}
			
			writer.write("\n");
			
		}
		
		writer.close();
	}
	
	private void write(String string, Map<Integer, Integer> diffs, int counter, Sequence reference) throws IOException {
		FileWriter writer=new FileWriter(string);
		char[] chars=reference.getSequenceBases().toCharArray();
		for(Integer pos: diffs.keySet()) {
			int count=diffs.get(pos);
			double freq=(double)count/(double)counter;
			writer.write(chars[pos]+"\t"+(pos+1)+"\t"+count+"\t"+freq+"\n");
		}
		
		writer.close();
	}
	
	
	private double difference(Sequence reference, Sequence other) {
		double counter=0;
		double total=0;
		char[] chars1=reference.getSequenceBases().toCharArray();
		char[] chars2=other.getSequenceBases().toCharArray();
		
		for(int i=0; i<chars1.length; i++) {
			if(chars1[i]==chars2[i]) {total++;}
			else {
				//counter++;
				if(chars2[i]!='?' && chars2[i]!='-' && chars2[i]!='N') {
					counter++;
					total++;
				}
				
			}
			//total++;
		}
		return counter/total;
	}

	private String differenceString(Sequence reference, Sequence other) {
		String rtrn="";
		char[] chars1=reference.getSequenceBases().toCharArray();
		char[] chars2=other.getSequenceBases().toCharArray();
		
		for(int i=0; i<chars1.length; i++) {
			if(chars1[i]==chars2[i]) {rtrn+=".";}
			else {
				rtrn+=chars2[i];
			}
		}
		return rtrn;
	}
	
	
	private void add(int i, Map<Integer, Map<Character, Integer>> variantCount, char c) {
		Map<Character, Integer> temp=new TreeMap<Character, Integer>();
		if(variantCount.containsKey(i)) {temp=variantCount.get(i);}
		int counter=0;
		if(temp.containsKey(c)) {counter=temp.get(c);}
		counter++;
		temp.put(c, counter);
		variantCount.put(i, temp);
	}

	private void add(int i, Map<Integer, Integer> diffs) {
		int count=0;
		if(diffs.containsKey(i)) {count=diffs.get(i);}
		count++;
		diffs.put(i, count);
	}

	public static void main(String[] args) throws IOException {
		File f=new File(args[0]);
		String save=args[1];
		new PercentDiff(f, save);
	}
	
}

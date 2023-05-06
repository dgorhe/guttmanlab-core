package guttmanlab.core.barcodeidentification;

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class EnumerateMismatches {

	
	
	private static Map<String, String> parse(String fastaFile) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		Collection<Sequence> seqs= FastaFileIOImpl.readFromFile(fastaFile);
		Map<String, String> mismatch1=new TreeMap<String, String>();
		
		for(Sequence seq: seqs){
			Map<String, String> mismatch=mismatch(seq);
			mismatch1.putAll(mismatch);
			rtrn.putAll(mismatch);
			rtrn.put(seq.getSequenceBases(), seq.getName());
		}
		
		
		//mismatch 2
		for(String seq: mismatch1.keySet()){
			Map<String, String> mismatch=mismatch(seq, mismatch1.get(seq));
			rtrn.putAll(mismatch);
		}
		
		return rtrn;
	}
	
	private static Map<String, String> mismatch(Sequence seq) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(int i=0; i<seq.getSequenceBases().length(); i++){
			String first=seq.getSequenceBases().substring(0, i);
			String second=seq.getSequenceBases().substring(i+1, seq.getSequenceBases().length());
			rtrn.put(first+"A"+second, seq.getName());
			rtrn.put(first+"C"+second, seq.getName());
			rtrn.put(first+"G"+second, seq.getName());
			rtrn.put(first+"T"+second, seq.getName());
			rtrn.put(first+"N"+second, seq.getName());
			//rtrn.put(first+second, seq.getName()); //Deletion
			//System.err.println(seq.getSequenceBases()+" "+i+" "+first+" "+second);
		}
		return rtrn;
	}
	
	private static Map<String, String> mismatch(String seq, String name) {
		Sequence s=new Sequence(name, seq);
		return mismatch(s);
	}
	
	
	public static void main(String[] args){
		String seq="TTCAACGTCCATGTCG";
		String name="pid";
		Sequence s=new Sequence(name, seq);
		Map<String, String> map=mismatch(s);
		for(String m: map.keySet()){
			Map<String, String> sub=mismatch(m, map.get(m));
			System.err.print(m);
			for(String m2: sub.keySet()){System.err.print("\t"+m2);}
			System.err.println();
			
		}
	}

	
}

package guttmanlab.core.sars.leader;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class GetMSACoordinates {

	public static Pair<Integer> getCoordinates(int referenceStart, int referenceEnd, Sequence referenceFasta) {
		
		Pair<Integer> rtrn=new Pair<Integer>();
		
		//go through reference sequence and count by excluding indels
		
		int counter=0;
		char[] chars=referenceFasta.getSequenceBases().toCharArray();
		for(int i=0; i<chars.length; i++) {
			if(chars[i]!='-') {counter++;}
			if(counter==referenceStart) {rtrn.setValue1(i);}
			if(counter==referenceEnd) {rtrn.setValue2(i);}
		}
		
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException {
		Sequence seq=FastaFileIOImpl.readFromFile(args[0]).iterator().next();
		int start=Integer.parseInt(args[1]);
		int end=Integer.parseInt(args[2]);
		Pair<Integer> coordinates=getCoordinates(start, end, seq);
		
		String sub=seq.getSubSequence(coordinates.getValue1(), coordinates.getValue2());
		System.err.println(coordinates.getValue1()+" "+coordinates.getValue2()+"\n"+sub);
	}
	
}

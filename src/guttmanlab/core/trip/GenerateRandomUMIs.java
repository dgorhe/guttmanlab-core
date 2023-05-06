package guttmanlab.core.trip;

import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.sequence.Sequence;

public class GenerateRandomUMIs {

	private static Collection<String> randomUMIs(int number, int length){
		Collection<String> rtrn=new TreeSet<String>();
		
		for(int i=0; i<number; i++) {
			String seq=Sequence.randomSequence(length);
			rtrn.add(seq);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) {
		
		int[] number= {1,10,100,1000, 5000, 10000, 25000, 50000, 100000,1000000};
		int length=10;
		int[] distance= {0};
		
		for(int i=0; i<number.length; i++) {
			Collection<String> list= randomUMIs(number[i], length);
			
			for(int j=0; j<distance.length; j++) {
				Collection<String> collapsed=CollapseCDNA.collapseUMI(list, distance[j]);
				System.out.println(number[i]+" "+distance[j]+"\t"+list.size()+" "+collapsed.size());
			}
		}
		
		
	}
	
}

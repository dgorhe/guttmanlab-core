package guttmanlab.core.test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

public class EnumeratePossibleSequences {

	public EnumeratePossibleSequences(List<String> dpm, List<String> round1, List<String> round2){
		for(String d1: dpm){
			for(String r1: round1){
				for(String r2: round2){
					String product=d1+r1+r2;
					System.out.println(product);
				}
			}
		}
	}
	
	
	public static void main(String[] args){
		Collection<String> kmers=enumerateKmer();
		String[] startingList={"CAAGTCA","CTCTCGT","GCTGTTG", "AGCAGAT","GATGTCA"};
		kmers=filterGC(kmers, 3);
		kmers=filterDistance(kmers, startingList, 4);
		kmers=excludeHomo(kmers);
		
		String kmer1=randomlySelect(kmers);
		kmers=filterDistance(kmers, kmer1, 4);
		
		String kmer2=randomlySelect(kmers);
		kmers=filterDistance(kmers, kmer2, 4);
		
		String kmer3=randomlySelect(kmers);
		kmers=filterDistance(kmers, kmer3, 4);
		
		String kmer4=randomlySelect(kmers);
		kmers=filterDistance(kmers, kmer4, 4);
		
		
		System.err.println(kmer1+"\n"+kmer2+"\n"+kmer3+"\n"+kmer4);
		
		
		System.err.println(kmers.size());
		
		
		/*List<String> dpm=new ArrayList<String>();
				dpm.add("AAGACCACCAGATCGGAAGAGCGTCGTGTA");
				dpm.add("ACAAGAGGCAGATCGGAAGAGCGTCGTGTA"); 
				dpm.add("ATCGAGGGTAGATCGGAAGAGCGTCGTGTA");
				dpm.add("AGAGGAGAAAGATCGGAAGAGCGTCGTGTA"); 
				dpm.add("AATACCTGGAGATCGGAAGAGCGTCGTGTA");
				dpm.add("ACTCTCTTAAGATCGGAAGAGCGTCGTGTA");
				
		List<String> round1=new ArrayList<String>();
			round1.add("CAAGTCAGCAGCCAC");
			round1.add("CAAGTCACGCAGCTG");
			round1.add("CAAGTCAAGCGCGTC"); 
			round1.add("CAAGTCACCGTGGGT");
			round1.add("CAAGTCAGGAGGCCA");
			round1.add("CAAGTCAGCGTCCCT");
			round1.add("CAAGTCACAGTCCGC");
			round1.add("CAAGTCACCGCCTCG");
			round1.add("CAAGTCATCGGCGGG");
			round1.add("CAAGTCACGGCTGCC");
			round1.add("CAAGTCAACCGCCCG"); 
			round1.add("CAAGTCAGGGCCGGT");
		
		List<String> round2=new ArrayList<String>();
			round2.add("AGTTGTCGGGTGCAG");
			round2.add("AGTTGTCCCTCGGCA");
			round2.add("AGTTGTCCGACACGC");
			round2.add("AGTTGTCGTACGGGC");
			round2.add("AGTTGTCAGCCCAGG");
			round2.add("AGTTGTCCCGACGAC");
			round2.add("AGTTGTCCCGGTCAG");
			round2.add("AGTTGTCGACGCAGC");
			round2.add("AGTTGTCCACACGCG");
			round2.add("AGTTGTCGTCCCGAG");
			round2.add("AGTTGTCACGGACGC");
			round2.add("AGTTGTCCTCGTCCC");
			
			new EnumeratePossibleSequences(dpm, round1, round2);*/
	}


	private static String randomlySelect(Collection<String> kmers) {
		int index=new Double(Math.random()*kmers.size()).intValue();
		return (String)kmers.toArray()[index];
	}


	private static Collection<String> excludeHomo(Collection<String> kmers) {
		Collection<String> rtrn=new TreeSet<String>();
		for(String kmer: kmers){
			if(kmer.contains("AAA") || kmer.contains("GGG")|| kmer.contains("CCC") || kmer.contains("TTT")){
				
			}
			else{
				//System.err.println(kmer);
				rtrn.add(kmer);
			}
		}
		return rtrn;
	}


	private static Collection<String> filterDistance(Collection<String> kmers, String[] startingList, int minDistance) {
		Collection<String> rtrn=new TreeSet<String>();
		for(String kmer: kmers){
			int distance=distance(kmer, startingList);
			if(distance>minDistance){
				rtrn.add(kmer);
			}
		}
		return rtrn;
	}
	
	private static Collection<String> filterDistance(Collection<String> kmers, String startingList, int minDistance) {
		Collection<String> rtrn=new TreeSet<String>();
		for(String kmer: kmers){
			int distance=distance(kmer, startingList);
			if(distance>minDistance){
				rtrn.add(kmer);
			}
		}
		return rtrn;
	}


	private static int distance(String kmer, String[] startingList) {
		int minDistance=7;
		for(int i=0; i<startingList.length; i++){
			minDistance=Math.min(minDistance, distance(kmer, startingList[i]));
		}
		return minDistance;
	}


	private static int distance(String kmer, String string) {
		int distance=0;
		char[] list1=kmer.toCharArray();
		char[] list2=string.toCharArray();
		
		for(int i=0; i<list1.length; i++){
			if(list1[i]!=list2[i]){distance++;}
		}
		return distance;
	}


	private static Collection<String> filterGC(Collection<String> kmers, int minGC) {
		Collection<String> rtrn=new TreeSet<String>();
		for(String kmer: kmers){
			if(countGC(kmer)>=minGC){rtrn.add(kmer);}
		}
		return rtrn;
	}


	private static int countGC(String kmer) {
		int count=0;
		char[] charList=kmer.toCharArray();
		for(int i=0; i<charList.length; i++){
			if(charList[i]=='C' || charList[i]=='G'){count++;}
		}
		return count;
	}


	private static Collection<String> enumerateKmer() {
		String[] list={"A", "C", "G", "T"};
		Collection<String> rtrn=new TreeSet<String>();
		
		for(int i1=0; i1<list.length; i1++){
			for(int i2=0; i2<list.length; i2++){
				for(int i3=0; i3<list.length; i3++){
					for(int i4=0; i4<list.length; i4++){
						for(int i5=0; i5<list.length; i5++){
							for(int i6=0; i6<list.length; i6++){
								for(int i7=0; i7<list.length; i7++){
									String kmer=list[i1]+list[i2]+list[i3]+list[i4]+list[i5]+list[i6]+list[i7];
									rtrn.add(kmer);
								}
						}
					}
				}
			}
		}
		}
		return rtrn;
	}
	
	
}

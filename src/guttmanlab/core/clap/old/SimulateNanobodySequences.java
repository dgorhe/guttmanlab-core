package guttmanlab.core.clap.old;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import guttmanlab.core.sequence.Sequence;

public class SimulateNanobodySequences {

	//CDR1 = NNKx7 = 21 nucleotides total
	//CDR3 = NNKx10 = 30 nucleotides total
	static char[] N= {'A', 'C', 'G', 'T'};
	static char[] K= {'G', 'T'}; //K = G,T
	
	public static List<String> randomSequences(int numberOfCodons, int uniqueSequences){
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<uniqueSequences; i++) {
			String seq=randomSequence(numberOfCodons);
			rtrn.add(seq);
		}
		return rtrn;
	}
	
	
	private static String randomSequence(int numberOfCodons) {
		String seq="";
		for(int i=0; i<numberOfCodons; i++) {
			String triple=sampleTriple();
			seq+=triple;
		}
		return seq;
	}

	private static String sampleTriple() {
		char char1=sample(N);
		char char2=sample(N);
		char char3=sample(K);
		String rtrn=char1+""+char2+""+char3;
		return rtrn;
	}

	private static char sample(char[] n2) {
		int rand=new Double(Math.random()*n2.length).intValue();
		return n2[rand];
	}
	
	private static void writeFastq(List<String> list, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		int counter=0;
		for(String seq: list) {
			String name="read"+counter;
			Sequence sequence=new Sequence(name, seq);
			writer.write(sequence.toFastq()+"\n");
			counter++;
		}
		writer.close();
	}
	
	private static int[] simulateSequencingDepth(List<String> sequences, int sequencingDepth, int numPerm) {
		//List<String> rtrn=new ArrayList<String>();
		
		Set<String>[] sets=new Set[numPerm];
		for(int i=0; i<numPerm; i++) {sets[i]=new TreeSet<String>();}
		
		//sample from sequences
		for(int i=0; i<sequencingDepth; i++) {
			String[] hits=sample(sequences, numPerm);
			for(int j=0; j<hits.length; j++) {
				sets[j].add(hits[j]);
			}
			//rtrn=sum(rtrn, hits);
		}
		
		int[] rtrn=new int[numPerm];
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=sets[i].size();
		}
		
		return rtrn;
	}
	
	private static int[] sum(int[] val1, int[] hits) {
		int[] rtrn=new int[val1.length];
		
		for(int i=0; i<rtrn.length; i++) {
			rtrn[i]=val1[i]+hits[i];
		}
		
		return rtrn;
	}


	private static String[] sample(List<String> sequences, int numPerm) {
		String[] rtrn=new String[numPerm];
		for(int i=0; i<numPerm; i++) {
			int randPos=new Double(Math.random()*sequences.size()).intValue();
			rtrn[i]=sequences.get(randPos);
		}
		return rtrn;
	}
	
	private static String sample(List<String> sequences) {
		int randPos=new Double(Math.random()*sequences.size()).intValue();
		return sequences.get(randPos);
	}

	
	private static int uniqueCount(List<String> list) {
		Set<String> rtrn=new TreeSet<String>();
		
		rtrn.addAll(list);
		
		return rtrn.size();
	}
	
	private static void write(double d, int sequencingDepth, int[] uniqueCounts) {
		String out=sequencingDepth+"\t"+d;
		for(int i=0; i<uniqueCounts.length; i++) {
			double fraction=(double)uniqueCounts[i]/(double)sequencingDepth;
			out+="\t"+fraction;
		}
		System.out.println(out);
	}

	public static void main(String[] args) throws IOException {
		if(args.length>0) {
			int codonLength=17;
			//int complexity=Integer.parseInt(args[0]);
			int sequencingDepth=Integer.parseInt(args[0]);
			//String save=args[2];
			int numPerm=10;
			
			double[] complexityFactors= {0.001, 0.01, 0.1, 1, 10, 100};
			
			for(int i=0; i<complexityFactors.length; i++) {
				System.err.println(complexityFactors[i]);
				int complexity=new Double(complexityFactors[i]*sequencingDepth).intValue();
				List<String> sequences=randomSequences(codonLength, complexity);
				int[] uniqueCounts=simulateSequencingDepth(sequences, sequencingDepth, numPerm);
				write(complexityFactors[i], sequencingDepth, uniqueCounts);
				
			}
			
			
			
			
				
		}
			
			
		else {System.err.println(usage);}
	}
	
	

	



	







	static String usage="args[0]=sequencing depth";
}

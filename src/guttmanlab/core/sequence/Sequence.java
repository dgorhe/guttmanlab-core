package guttmanlab.core.sequence;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.nanobody.Cluster;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * A nucleotide sequence
 * @author prussell
 *
 */
public class Sequence implements Comparable<Sequence>{
	
	private String sequence;
	private String name;
	
	/**
	 * @param name Sequence name
	 * @param seq Nucleotide sequence
	 */
	public Sequence(String name, String seq){
		this.name=name;
		this.sequence=seq;
	}
	
	/**
	 * @param seq Nucleotide sequence
	 */
	public Sequence(String seq){
		this.sequence=seq;
	}
	
	/**
	 * @return A new sequence that is the reverse complement of this sequence
	 */
	public Sequence reverseComplement() {
		Sequence tmpSeq = new Sequence(this.name, this.sequence);
		tmpSeq.reverse();
		return tmpSeq;
	}
	
	public static String reverseComplement(String sequence) {
		Sequence tmpSeq = new Sequence(sequence);
		tmpSeq.reverse();
		return tmpSeq.sequence;
	}
	
	@Override
	public String toString(){
		return this.sequence;
	}
	
	public String toFastq() {
		return "@"+name+"\n"+this.sequence+"\n+\n"+quality(this.sequence);
	}
	
	private String quality(String sequence2) {
		String rtrn="";
		for(int i=0; i<sequence2.length(); i++) {
			rtrn+="I";
		}
		return rtrn;
	}

	public String toString(boolean mask){
		if(!mask){return toString();}
		else{
			String seq=getSequenceBases();
			seq.replaceAll("a", "N");
			seq.replaceAll("c", "N");
			seq.replaceAll("g", "N");
			seq.replaceAll("t", "N");
			return seq;
		}
		
	}
	
	private void reverse() {
		String seqBases = getSequenceBases();
		String reversedSeq = "";
		for(int j = seqBases.length() - 1; j >= 0 ; j--) {
			char c = seqBases.charAt(j);
			if('c' == c) {
				reversedSeq+='g';
			}else if ('C' == c) {
				reversedSeq+='G';
			}else if ('G' == c) {
				reversedSeq+='C';
			}else if ('g' == c) {
				reversedSeq+=('c');
			}else if ('a' == c) {
				reversedSeq+=('t');
			}else if ('A' == c) {
				reversedSeq+=('T');
			}else if('t' == c) {
				reversedSeq+=('a');
			}else if('T' == c) {
				reversedSeq+=('A');
			}else if('N'==c){
				reversedSeq+=('N');
			}else if('n'==c){
				reversedSeq+=('n');
			}else {
				reversedSeq+=(c);
			}
		}
		
		this.sequence = reversedSeq;
	}
	
	/**
	 * @return Return the sequence bases
	 */
	public String getSequenceBases() {
		return this.sequence;
	}
	
	/**
	 * @return Sequence name
	 */
	public String getName(){
		return this.name;
	}
	
	/**
	 * @return Sequence length
	 */
	public int getLength() {
		return sequence.length();
	}
	
	/**
	 * Get subsequence
	 * @param name Name of new sequence to return
	 * @param start Start position of subsequence
	 * @param end Position after last position to include
	 * @return The subsequence
	 */
	public Sequence getSubSequence(String name, int start, int end) {
		String subSeq = sequence.substring(Math.max(start, 0), Math.min(end, sequence.length()));
		Sequence seq = new Sequence(name, subSeq);
		return seq;
	}
	
	public String getSubSequence(int start, int end) {
		String subSeq = sequence.substring(Math.max(start, 0), Math.min(end, sequence.length()));
		return subSeq;
	}
	
	/**
	 * Get the spliced transcribed sequence of an annotation
	 * Bases are reported in 5' to 3' direction
	 * @param annot The annotation
	 * @return Sequence with same name as annotation containing the transcribed sequence
	 */
	public Sequence getSubsequence(Annotation annot) {
		if(!annot.getOrientation().equals(Strand.POSITIVE) && !annot.getOrientation().equals(Strand.NEGATIVE)) {
			throw new IllegalArgumentException("Strand must be known");
		}
		
		
		if(annot.getOrientation().equals(Strand.NEGATIVE)) {return getNegativeSubsequence(annot);}
		
		return getPositiveSubseqence(annot);
		
		
	}
	
	private Sequence getNegativeSubsequence(Annotation annot) {
		String seq = "";
		Iterator<SingleInterval> blockIter = annot.getBlocks();
		
		SingleInterval[] list=new SingleInterval[annot.getNumberOfBlocks()];
		
		int i=list.length-1;
		while(blockIter.hasNext()) {
			list[i]=blockIter.next();
			i--;
		}
		
		for(SingleInterval block: list) {
			Sequence blockSequence = getSubSequence("", block.getReferenceStartPosition(), block.getReferenceEndPosition());
			Sequence rc=blockSequence.reverseComplement();
			seq += rc;
		}
		
		
		return new Sequence(annot.getName(), seq);
	}

	private Sequence getPositiveSubseqence(Annotation annot) {
		String seq = "";
		Iterator<SingleInterval> blockIter = annot.getBlocks();
		while(blockIter.hasNext()) {
			SingleInterval block = blockIter.next();
			Sequence blockSequence = getSubSequence("", block.getReferenceStartPosition(), block.getReferenceEndPosition());
			
			String forwardBases = blockSequence.getSequenceBases();
			seq += forwardBases;
		}
		
		return new Sequence(annot.getName(), seq);
	}

	/**
	 * Soft mask the specified regions
	 * Throw an exception if the annotations do not refer to this sequence
	 * @param regions The regions to mask with respect to this sequence
	 * @return A new sequence with the regions soft masked
	 */
	public Sequence softMask(Collection<Annotation> regions) {
		// TODO
		throw new UnsupportedOperationException();
	}
	
	/**
	 * Hard mask the specified regions
	 * Throw an exception if the annotations do not refer to this sequence
	 * @param regions The regions to mask with respect to this sequence
	 * @return A new sequence with the regions hard masked
	 */
	public Sequence hardMask(Collection<Annotation> regions) {
		// TODO
		throw new UnsupportedOperationException("TODO");
	}
	
	@Override
	public boolean equals(Object o) {
		if(!o.getClass().equals(Sequence.class)) {
			return false;
		}
		Sequence otherSeq = (Sequence)o;
		if(!getName().equals(otherSeq.getName()))	{
			return false;
		}
		if(!getSequenceBases().equals(otherSeq.getSequenceBases()))	{
			return false;
		}
		return true;
	}
	
	@Override
	public int hashCode() {
		String s = getName() + ":" + getSequenceBases();
		return s.hashCode();
	}

	public Collection<SingleInterval> find(Sequence reSeq) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		int length=reSeq.getSequenceBases().length();
		for(int i=0; i<getSequenceBases().length(); i++){
			String sub=getSubSequence(i, i+length);
			if(sub.equalsIgnoreCase(reSeq.getSequenceBases()) || sub.equalsIgnoreCase(reSeq.reverseComplement().getSequenceBases())){
				SingleInterval interval=new SingleInterval(this.getName(), i, i+length);
				//System.err.println(reSeq.getName()+" "+interval.toUCSC()+" "+sub);
				rtrn.add(interval);
			}
		}
		return rtrn;
		
	}

	@Override
	public int compareTo(Sequence o) {
		if(!this.getName().equals(o.getName())){return this.getName().compareTo(o.getName());}
		if(!this.getSequenceBases().equals(o.getSequenceBases())){return this.getSequenceBases().compareTo(o.getSequenceBases());}
		return 0;
	}

	public static double computeGCContent(String productSequence) {
		int count=0;
		int total=0;
		for(int i=0; i< productSequence.toCharArray().length; i++){
			if(productSequence.toCharArray()[i]=='C' || productSequence.toCharArray()[i]=='G'){count++;}
			total++;
		}
		return (double)count/(double)total;
	}

	public int getLength() {
		return this.sequence.length();
	}

	public String toFasta() {
		return ">"+this.getName()+"\n"+this.sequence;
	}

	public void setName(String ucsc) {
		this.name=ucsc;
		
	}

	public static String trim(String seq) {
		String rtrn="";
		char[] chars=seq.toCharArray();
		
		for(int i=0; i<chars.length; i++) {
			if(chars[i]!='-') {rtrn+=chars[i];}
		}
		return rtrn;
	}
	
	public String translate() {
		String rtrn="";
		//go through nucleotide sequence and get triples.
		//compare to genetic code
		boolean hitEnd=false;
		
		String trim=trim(this.getSequenceBases());
		Map<String, String> codonTable=getCodonTable();
		
		char[] chars=trim.toCharArray();
		for(int i=0; i<chars.length-3; i+=3) {
			if(!hitEnd) {
				String codon=chars[i]+""+chars[i+1]+""+chars[i+2];
				String AA=codonTable.get(codon);
				if(AA==null) {AA="?";}
				if(AA.equals("*")) {hitEnd=true;}
				rtrn+=AA;
			}
		}
		
		return rtrn;
	}
	
	public String getCodingSequence() {
		//go through nucleotide sequence and get triples.
		//compare to genetic code
		boolean hitEnd=false;
		boolean hitStart=false;
		int startPos=0;
		int endPos=0;
		
		String trim=trim(this.getSequenceBases());
		Map<String, String> codonTable=getCodonTable();
		
		char[] chars=trim.toCharArray();
		for(int i=0; i<chars.length-3; i++) {
			String codon=chars[i]+""+chars[i+1]+""+chars[i+2];
			if(!hitStart) {
				if(codon.equals("ATG")) {
					hitStart=true;
					startPos=i;
				}
			}
		}
		
		for(int i=startPos; i<chars.length-3; i+=3) {
			String codon=chars[i]+""+chars[i+1]+""+chars[i+2];
			if(!hitEnd) {
				String AA=codonTable.get(codon);
				System.err.println(i+" "+codon+" "+AA);
				if(AA.equals("*")) {
					hitEnd=true;
					endPos=i;
				}
			}
		}
		
			
		
		
		System.err.println(startPos+" "+endPos);
		return this.getSequenceBases().substring(startPos, endPos);
	}

	public static Map<String, String> getCodonTable() {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		rtrn.put("TTT","F");
		rtrn.put("TTC","F");
		rtrn.put("TTA","L");
		rtrn.put("TTG","L");
		rtrn.put("CTT","L");
		rtrn.put("CTC","L");
		rtrn.put("CTA","L");
		rtrn.put("CTG","L");
		rtrn.put("ATT","I");
		rtrn.put("ATC","I");
		rtrn.put("ATA","I");
		rtrn.put("ATG","M");
		rtrn.put("GTT","V");
		rtrn.put("GTC","V");
		rtrn.put("GTA","V");
		rtrn.put("GTG","V");
		rtrn.put("TCT","S");
		rtrn.put("TCC","S");
		rtrn.put("TCA","S");
		rtrn.put("TCG","S");
		rtrn.put("CCT","P");
		rtrn.put("CCC","P");
		rtrn.put("CCA","P");
		rtrn.put("CCG","P");
		rtrn.put("ACT","T");
		rtrn.put("ACC","T");
		rtrn.put("ACA","T");
		rtrn.put("ACG","T");
		rtrn.put("GCT","A");
		rtrn.put("GCC","A");
		rtrn.put("GCA","A");
		rtrn.put("GCG","A");
		rtrn.put("TAT","Y");
		rtrn.put("TAC","Y");
		rtrn.put("TAA","*");
		rtrn.put("TAG","*");
		rtrn.put("CAT","H");
		rtrn.put("CAC","H");
		rtrn.put("CAA","Q");
		rtrn.put("CAG","Q");
		rtrn.put("AAT","N");
		rtrn.put("AAC","N");
		rtrn.put("AAA","K");
		rtrn.put("AAG","K");
		rtrn.put("GAT","D");
		rtrn.put("GAC","D");
		rtrn.put("GAA","E");
		rtrn.put("GAG","E");
		rtrn.put("TGT","C");
		rtrn.put("TGC","C");
		rtrn.put("TGA","*");
		rtrn.put("TGG","W");
		rtrn.put("CGT","R");
		rtrn.put("CGC","R");
		rtrn.put("CGA","R");
		rtrn.put("CGG","R");
		rtrn.put("AGT","S");
		rtrn.put("AGC","S");
		rtrn.put("AGA","R");
		rtrn.put("AGG","R");
		rtrn.put("GGT","G");
		rtrn.put("GGC","G");
		rtrn.put("GGA","G");
		rtrn.put("GGG","G");
		
		
		return rtrn;
	}

	/*public ArrayList<String> enumerateKmer(int k) {
		ArrayList<String> rtrn=new ArrayList<String>();
		
		char[] chars=this.sequence.toCharArray();
		
		for(int i=0; i<chars.length-k; i++) {
			String val=sub(chars, i, i+k);
			//System.err.println(i+" "+val);
			rtrn.add(val);
		}
		
		return rtrn;
	}*/
	
	public Map<Integer, String> enumerateKmer(int k) {
		Map<Integer, String> rtrn=new TreeMap<Integer, String>();
		
		char[] chars=this.sequence.toCharArray();
		
		for(int i=0; i<chars.length-k; i++) {
			String val=sub(chars, i, i+k);
			//System.err.println(i+" "+val);
			rtrn.put(i, val);
		}
		
		return rtrn;
	}
	
	public static Map<Integer, String> enumerateKmer(String seq, int k) {
		Map<Integer, String> rtrn=new TreeMap<Integer, String>();
		
		char[] chars=seq.toCharArray();
		
		for(int i=0; i<=chars.length-k; i++) {
			String val=sub(chars, i, i+k);
			rtrn.put(i, val);
		}
		
		return rtrn;
	}

	private static String sub(char[] chars, int i, int j) {
		
		StringBuilder rtrn=new StringBuilder();
		for(int index=i; index<j; index++) {rtrn.append(chars[index]);}
		
		return rtrn.toString();
	}

	public Map<Integer, List<String>> enumerateKmer(int[] primerSizes) {
		Map<Integer, List<String>> rtrn=new TreeMap<Integer, List<String>>();
		for(int i=0; i< primerSizes.length ; i++) {
			Map<Integer, String> kmer=enumerateKmer(primerSizes[i]);
			add(kmer, rtrn);
		}
		return rtrn;
	}

	private void add(Map<Integer, String> kmer, Map<Integer, List<String>> rtrn) {
		for(Integer pos: kmer.keySet()) {
			if(!rtrn.containsKey(pos)) {
				List<String> newList=new ArrayList<String>();
				rtrn.put(pos, newList);
			}
			rtrn.get(pos).add(kmer.get(pos));
		}
		
	}

	public static String translate(String rCDR2) {
		char[] seqs=rCDR2.toCharArray();
		
		Map<String, String> codonTable=getCodonTable();
		
		String rtrn="";
		for(int i=0; i<seqs.length-2; i+=3) {
			String codon=seqs[i]+""+seqs[i+1]+""+seqs[i+2];
			rtrn+=codonTable.get(codon);
		}
		
		return rtrn;
	}
	
	private static boolean isClose(String seq1, String seq2, int minDistance) {
		char[] char1=seq1.toCharArray();
		char[] char2=seq2.toCharArray();
		
		int distance=0;
		for(int i=0; i<char1.length; i++) {
			if(char1[i]!=char2[i]) {
				distance++;
				if(distance>minDistance) {return false;}
			}
		}
		return true;
	}

	public static Collection<Cluster> cluster(Collection<String> list, int minDistance) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		Map<String, Cluster> toMerge=new TreeMap<String, Cluster>();
		
		//iterate through each pair
		for(String str1: list) {
			Cluster c1=new Cluster(str1);
			for(String str2: list) {
				if(isClose(str1, str2, minDistance)) {
					c1.addProtein(str2);
				}
			}
			if(c1.getSize()==1) {rtrn.add(c1);}
			else {toMerge.put(str1, c1);}
		}
		
		rtrn.addAll(merge(toMerge));
		return rtrn;
	}

	private static Collection<Cluster> merge(Map<String, Cluster> toMerge) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(String b: toMerge.keySet()) {
			Cluster c=toMerge.get(b);
			Cluster newC=new Cluster(c.getProteins());
			for(String p: c.getProteins()) {
				newC.addProteins(toMerge.get(p).getProteins());
			}
			rtrn.add(newC);
		}
		
		return merge(rtrn);
	}

	private static Collection<Cluster> merge(Collection<Cluster> toMerge) {
		
		Collection<Cluster> subset=new TreeSet<Cluster>();
		int maxSize=getMaxSize(toMerge);
		
		for(int i=2; i<maxSize; i++) {
			Collection<Cluster> clusters1=getClusters(toMerge, i);
			Collection<Cluster> clusters2=getClustersAbove(toMerge, i);
			
			for(Cluster c1: clusters1) {
				for(Cluster c2: clusters2) {
					if(c1.isSubset(c2)) {subset.add(c1);}
				}
			}
		}
		
		toMerge.removeAll(subset);
		return toMerge;
	}

	private static Collection<Cluster> getClustersAbove(Collection<Cluster> toMerge, int i) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Cluster c: toMerge) {
			if(c.getSize()>i) {rtrn.add(c);}
		}
		
		return rtrn;
	}

	private static Collection<Cluster> getClusters(Collection<Cluster> toMerge, int i) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		
		for(Cluster c: toMerge) {
			if(c.getSize()==i) {rtrn.add(c);}
		}
		
		return rtrn;
	}

	private static int getMaxSize(Collection<Cluster> toMerge) {
		int max=0;
		
		for(Cluster c: toMerge) {
			max=Math.max(c.getSize(), max);
		}
		
		return max;
	}

	public static int distance(String seq1, String seq2) {
		char[] char1=seq1.toCharArray();
		char[] char2=seq2.toCharArray();
			
		int distance=0;
		for(int i=0; i<char1.length; i++) {
			if(char1[i]!=char2[i]) {distance++;}
		}
		return distance;
	}

	public static String randomSequence(int length) {
		String rtrn="";
		
		for(int i=0; i<length; i++) {
			rtrn+=randomNucleotide();
		}
		
		return rtrn;
	}

	private static String randomNucleotide() {
		double rand=Math.random();
		if(rand<0.25) {return "A";}
		if(rand>0.25 && rand<0.5) {return "C";}
		if(rand>0.5 && rand<0.75) {return "G";}
		return "T";
	}

	

	
	
}

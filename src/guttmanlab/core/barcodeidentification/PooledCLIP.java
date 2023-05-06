package guttmanlab.core.barcodeidentification;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;

public class PooledCLIP {
	
	String proteinTag="TTCAACGTCCATGTCG";
	double threshold=0.9;
	Map<String, String> b8;
	Map<String, String> b7;
	Map<String, String> b6;
	Map<String, String> b5;
	Map<String, String> b4;
	Map<String, String> b3;
	Map<String, String> b2;
	Map<String, String> pid;
	
	Map<String, String> b8m1;
	Map<String, String> b7m1;
	Map<String, String> b6m1;
	Map<String, String> b5m1;
	Map<String, String> b4m1;
	Map<String, String> b3m1;
	Map<String, String> b2m1;
	Map<String, String> pidm1;
	
	Map<String, String> b8m2;
	Map<String, String> b7m2;
	Map<String, String> b6m2;
	Map<String, String> b5m2;
	Map<String, String> b4m2;
	Map<String, String> b3m2;
	Map<String, String> b2m2;
	Map<String, String> pidm2;
	
	int tagLength=16;

	public PooledCLIP(File file1, File file2, Map<String, String> b8, Map<String, String> b7, Map<String, String> b6, Map<String, String> b5, Map<String, String> b4, Map<String, String> b3, Map<String, String> b2, Map<String, String> proteinIDs, String save) throws IOException{
		this.b8=b8;
		this.b7=b7;
		this.b6=b6;
		this.b5=b5;
		this.b4=b4;
		this.b3=b3;
		this.b2=b2;
		this.pid=proteinIDs;
		mismatch1();
		mismatch2();
		
		
		FastqReader f1=new FastqReader(file1);
		FastqReader f2=new FastqReader(file2);
		
		Iterator<FastqRecord> iter1=f1.iterator();
		Iterator<FastqRecord> iter2=f2.iterator();
		
		
		FastqWriter writer1=new FastqWriterFactory().newWriter(new File(save+"_R1.fq"));
		FastqWriter writer2=new FastqWriterFactory().newWriter(new File(save+"_R2.fq"));
		
		int counter=0;
		int missed=0;
		
		long start=System.currentTimeMillis();
		while(iter1.hasNext()){
			FastqRecord record1=iter1.next();
			FastqRecord record2=iter2.next();
			
			String barcode=parseBarcodes(record2, b8, b7, b6, b5, b4, b3, b2); //TODO Cleanup this implementation
			
			String pid=findProteinTag(record1.getReadString(), proteinIDs, proteinTag);
			barcode="["+pid+"]"+barcode;
			
			String header1=removeSpace(record1.getReadHeader())+"::"+barcode;
			String header2=removeSpace(record2.getReadHeader())+"::"+barcode;
			
			//String seqHeaderPrefix, java.lang.String seqLine, java.lang.String qualHeaderPrefix, java.lang.String qualLine
			FastqRecord newRecord1=new FastqRecord(header1, record1.getReadString(), record1.getBaseQualityHeader(), record1.getBaseQualityString());
			writer1.write(newRecord1);
			
			
			FastqRecord newRecord2=new FastqRecord(header2, record2.getReadString(), record2.getBaseQualityHeader(), record2.getBaseQualityString());
			writer2.write(newRecord2);
			
			if(barcode.contains("NF")) {
				missed++;
				
				/*String header1=removeSpace(record1.getReadHeader())+"::"+barcode;
				String header2=removeSpace(record2.getReadHeader())+"::"+barcode;
				
				//String seqHeaderPrefix, java.lang.String seqLine, java.lang.String qualHeaderPrefix, java.lang.String qualLine
				FastqRecord newRecord1=new FastqRecord(header1, record1.getReadString(), record1.getBaseQualityHeader(), record1.getBaseQualityString());
				writer1.write(newRecord1);
				
				
				FastqRecord newRecord2=new FastqRecord(header2, record2.getReadString(), record2.getBaseQualityHeader(), record2.getBaseQualityString());
				writer2.write(newRecord2);*/
				
				
			}
			
			counter++;
			if(counter%10000==0){
				long end=System.currentTimeMillis();
				double sec=(end-start)/1000.0;
				start=end;
				System.err.println(missed+" "+counter+" "+sec);
			}
			
		}
		
		System.err.println(missed+" "+counter);
		f1.close();
		f2.close();
		writer1.close();
		writer2.close();
	}

	private void mismatch1() {
		b8m1=new TreeMap<String, String>();
		b7m1=new TreeMap<String, String>();
		b6m1=new TreeMap<String, String>();
		b5m1=new TreeMap<String, String>();
		b4m1=new TreeMap<String, String>();
		b3m1=new TreeMap<String, String>();
		b2m1=new TreeMap<String, String>();
		pidm1=new TreeMap<String, String>();
		
		b8m1.putAll(mismatch(b8, true));
		b7m1.putAll(mismatch(b7, false));
		b6m1.putAll(mismatch(b6, false));
		b5m1.putAll(mismatch(b5, false));
		b4m1.putAll(mismatch(b4, false));
		b3m1.putAll(mismatch(b3, false));
		b2m1.putAll(mismatch(b2, false));
		pidm1.putAll(mismatch(pid, false));
		
	}
	
	private void mismatch2() {
		b8m2=new TreeMap<String, String>();
		b7m2=new TreeMap<String, String>();
		b6m2=new TreeMap<String, String>();
		b5m2=new TreeMap<String, String>();
		b4m2=new TreeMap<String, String>();
		b3m2=new TreeMap<String, String>();
		b2m2=new TreeMap<String, String>();
		pidm2=new TreeMap<String, String>();
		
		b8m2.putAll(mismatch(b8m1, true));
		b7m2.putAll(mismatch(b7m1, false));
		b6m2.putAll(mismatch(b6m1, false));
		b5m2.putAll(mismatch(b5m1, false));
		b4m2.putAll(mismatch(b4m1, false));
		b3m2.putAll(mismatch(b3m1, false));
		b2m2.putAll(mismatch(b2m1, false));
		pidm2.putAll(mismatch(pidm1, false));
		
	}

	private Map<String, String> mismatch(Map<String, String> b82, boolean withN) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String seq: b82.keySet()) {
			String name=b82.get(seq);
			Map<String, String> mismatch=mismatch(seq, name, withN);
			rtrn.putAll(mismatch);
		}
		
		return rtrn;
	}

	private String findProteinTag(String readString, Map<String, String> proteinIDs, String proteinTag2) {
		//String rtrn=searchHarder(readString, proteinIDs, pidm1, pidm2, "NF");
		
		String rtrn=getProteinID(readString, proteinIDs, pidm1, pidm2, "NF");
		
		if(rtrn.equals("NF")){
			if(matches(readString, proteinTag2)){rtrn="PROT_unassigned";}
			else{rtrn="RNA";}
		}
		else{rtrn="PROT_"+rtrn;}
		
		return rtrn;
	}

	private String getProteinID(String readString, Map<String, String> proteinIDs, Map<String, String> pidm1, Map<String, String> pidm2, String NF) {
		int startIndex=38;
		int endIndex=startIndex+tagLength;
		String subseq=readString.substring(startIndex, endIndex);
		if(proteinIDs.containsKey(subseq)) {return proteinIDs.get(subseq);}
		if(pidm1.containsKey(subseq)) {return pidm1.get(subseq);}
		if(pidm2.containsKey(subseq)) {return pidm2.get(subseq);}
		
		int fudge=3;
		//TODO Search up and down
		for(int i=1; i<fudge; i++) {
			if(startIndex-i>=0) {
				String kmer1=readString.substring(startIndex-i, endIndex-i);
				if(proteinIDs.containsKey(kmer1)) {return "["+proteinIDs.get(kmer1)+"]";}
				//if(pidm1.containsKey(kmer1)) {return "["+pidm1.get(kmer1)+"]";}
				//if(pidm2.containsKey(kmer1)) {return "["+pidm2.get(kmer1)+"]";}
			}
			if(endIndex+i<readString.length()) {
				String kmer2=readString.substring(startIndex+i, endIndex+i);
				if(proteinIDs.containsKey(kmer2)) {return "["+proteinIDs.get(kmer2)+"]";}
				//if(pidm1.containsKey(kmer2)) {return "["+pidm1.get(kmer2)+"]";}
				//if(pidm2.containsKey(kmer2)) {return "["+pidm2.get(kmer2)+"]";}
			}
		}
		
		return NF;
		
	}

	private boolean matches(String readString, String proteinTag2) {
		int size=proteinTag2.length();
		for(int i=0; i<readString.length()-size; i++){
			String kmer=readString.substring(i, i+proteinTag2.length());
			if((kmer.equalsIgnoreCase(proteinTag2))){
				return true;
			}
		}
		return false;
	}
	

	private String removeSpace(String readHeader) {
		return readHeader.replaceAll(" ", ":");
	}

	private String assign(Map<String, Integer> counts) {
		double total=0;
		for(String name: counts.keySet()){
			total+=counts.get(name);
		}
		
		for(String name: counts.keySet()){
			double score=counts.get(name);
			double ratio=score/total;
			if(ratio>threshold){return name;}
		}
		return "NF";
	}

	private Map<String, Integer> count(Collection<String> collection) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(String val: collection){
			String protein=val.split("_")[0];
			int count=0;
			if(rtrn.containsKey(protein)){count=rtrn.get(protein);}
			count++;
			rtrn.put(protein, count);
		}
		
		return rtrn;
	}

	private String getProteinID(FastqRecord record1, Map<String, String> proteinIDs) {
		String seq=record1.getReadString();
		String rtrn="NF";
		if(seq.startsWith(proteinTag)){
			rtrn="["+searchHarder(seq, proteinIDs, "NF");
			String UMI=seq.substring(22, 30)+"]";
			rtrn+="_"+UMI;
			
		}
		return rtrn;
	}

	private String parseBarcodes(FastqRecord record2, Map<String, String> b8, Map<String, String> b7, Map<String, String> b6, Map<String, String> b5, Map<String, String> b4, Map<String, String> b3, Map<String, String> b2) {
		String seq=record2.getReadString();
		
		int index=getR8Index(record2, b8, b8m1, b8m2);
		
		//if(index==-1){index=getR7Index(seq, b7);}
		
		String rtrn="NF_R8";
		if(index==-1){index=8;}
		//if(index!=-1){
		
			String R8=get(seq,0, index, b8, b8m1, b8m2, "NF_R8");
			
			//String R8=seq.substring(0, index);
			
			String R7=get(seq, index+7, index+23, b7, b7m1, b7m2, "NF_R7");
			
			//System.err.println(record2.getReadString()+" "+R7);
			
			String R6=get(seq, index+30, index+46, b6, b6m1, b6m2, "NF_R6");
			
			
			String R5=get(seq,index+53, index+69, b5, b5m1, b5m2, "NF_R5");
			
			String R4=get(seq, index+76, index+92, b4, b4m1, b4m2, "NF_R4");
			
			String R3=get(seq, index+99, index+115, b3, b3m1, b3m2, "NF_R3");
			
			String R2=get(seq, index+122, index+138, b2, b2m1, b2m2, "NF_R2");
			
			
			rtrn=R8+R7+R6+R5+R4+R3+R2;
		
			//System.err.println("Round 1: "+rtrn);
			
			/*if(rtrn.contains("NF")) {
				
				for(int i=0; i<seq.length()-tagLength; i++){
					String kmer=seq.substring(i, i+tagLength);
					if(b2.containsKey(kmer)){R2="["+b2.get(kmer)+"]";}
					if(b3.containsKey(kmer)) {R3="["+b3.get(kmer)+"]";}
					if(b4.containsKey(kmer)) {R4="["+b4.get(kmer)+"]";}
					if(b5.containsKey(kmer)) {R5="["+b5.get(kmer)+"]";}
					if(b6.containsKey(kmer)) {R6="["+b6.get(kmer)+"]";}
					if(b7.containsKey(kmer)) {R7="["+b7.get(kmer)+"]";}
				}
			}
		
		rtrn=R8+R7+R6+R5+R4+R3+R2;*/
			
		//System.err.println("Round 2: "+rtrn);
		
		return rtrn;
	}

	private int getR7Index(String seq, Map<String, String> barcodeNames) {
		if(barcodeNames.containsKey(seq.substring(8+7, 8+23))){
			return 8;
		}
		if(barcodeNames.containsKey(seq.substring(9+7, 9+23))){
			return 9;
		}
		if(barcodeNames.containsKey(seq.substring(10+7, 10+23))){
			return 10;
		}
		if(barcodeNames.containsKey(seq.substring(11+7, 11+23))){
			return 11;
		}
		if(barcodeNames.containsKey(seq.substring(12+7, 12+23))){
			return 12;
		}
		return -1;
	}

	private String get(String seq, int start, int end, Map<String, String> b2) {
		String name="[NF_"+seq.substring(start, end)+"]";
		
		String substring=seq.substring(start, end);
		if(b2.containsKey(substring)){
			return "["+b2.get(substring)+"]";
		}
		
		name=searchHarder(seq, b2, name);
		
		return name;
	
	}
	
	
	private String get(String seq, int start, int end, Map<String, String> b2, Map<String, String> b2m1, Map<String, String> b2m2, String NF) {
		String name=NF+"_"+seq.substring(start, end);
		
		String substring=seq.substring(start, end);
		if(b2.containsKey(substring)){
			return "["+b2.get(substring)+"]";
		}
		
		if(b2m1.containsKey(substring)){
			return "["+b2m1.get(substring)+"]";
		}
		
		if(b2m2.containsKey(substring)){
			return "["+b2m2.get(substring)+"]";
		}
		
		int fudge=5;
		//TODO Search up and down
		for(int i=1; i<fudge; i++) {
			if(start-i>=0) {
				String kmer1=seq.substring(start-i, end-i);
				if(b2.containsKey(kmer1)) {return "["+b2.get(kmer1)+"]";}
				if(b2m1.containsKey(kmer1)) {return "["+b2m1.get(kmer1)+"]";}
				//if(b2m2.containsKey(kmer1)) {return "["+b2m2.get(kmer1)+"]";}
			}
			if(end+i<seq.length()) {
				String kmer2=seq.substring(start+i, end+i);
				if(b2.containsKey(kmer2)) {return "["+b2.get(kmer2)+"]";}
				if(b2m1.containsKey(kmer2)) {return "["+b2m1.get(kmer2)+"]";}
				//if(b2m2.containsKey(kmer2)) {return "["+b2m2.get(kmer2)+"]";}
			}
		}
		
		
		//name=searchHarder(seq, b2, b2m1, b2m2, name);
		name="["+name+"]";
		
		return name;
	
	}
		
		
		

	private String searchHarder(String seq, Map<String, String> b2, String NF) {
		Collection<Integer> sizes=new TreeSet<Integer>();
			for(String b: b2.keySet()){sizes.add(b.length());}
			
			for(int i=0; i<seq.length(); i++){
				Collection<String> kmers=getKmer(seq, i, sizes);
				for(String kmer: kmers){
					//System.err.println(kmer+" "+kmer.length());
					if(b2.containsKey(kmer)){
						String name=b2.get(kmer);
						//System.err.println(i+" "+name);
						return name;
					}
				}
			}
			return NF;
		}
		
	
	private String searchHarder(String seq, Map<String, String> b2, Map<String, String> b2m1, Map<String, String> b2m2, String NF) {
		Collection<Integer> sizes=new TreeSet<Integer>();
			for(String b: b2.keySet()){sizes.add(b.length());}
			
			for(int i=0; i<seq.length(); i++){
				Collection<String> kmers=getKmer(seq, i, sizes);
				for(String kmer: kmers){
					//System.err.println(kmer+" "+kmer.length());
					if(b2.containsKey(kmer)){
						String name=b2.get(kmer);
						System.err.println(i+" "+name+" "+kmer.length());
						return name;
					}
					/*else if(b2m1.containsKey(kmer)) {
						String name=b2m1.get(kmer);
						//System.err.println(i+" "+name);
						return name;
					}
					else if(b2m2.containsKey(kmer)) {
						String name=b2m2.get(kmer);
						//System.err.println(i+" "+name);
						return name;
					}*/
				}
			}
			return NF;
		}
	

	private Collection<String> getKmer(String seq, int i, Collection<Integer> sizes) {
		Collection<String> rtrn=new TreeSet<String>();
		for(Integer size: sizes){
			if(i+size<seq.length()){
				rtrn.add(seq.substring(i, i+size));
			}
		}
		return rtrn;
	}
	
	

	private int getR8Index(FastqRecord record2, Map<String, String> barcodeNames) {
		String seq=record2.getReadString();
		
		String R8=seq.substring(0,  12);
		if(barcodeNames.containsKey(R8)){return 12;}
		R8=seq.substring(0,  11);
		if(barcodeNames.containsKey(R8)){return 11;}
		R8=seq.substring(0,  10);
		if(barcodeNames.containsKey(R8)){return 10;}
		R8=seq.substring(0,  9);
		if(barcodeNames.containsKey(R8)){return 9;}
		R8=seq.substring(0, 8);
		if(barcodeNames.containsKey(R8)){return 8;}
		R8=seq.substring(0, 5);
		if(barcodeNames.containsKey(R8)){return 5;}
		
		return -1;
	}

	
	private int getR8Index(FastqRecord record2, Map<String, String> b8, Map<String, String> b8m1, Map<String, String> b8m2) {
		String seq=record2.getReadString();
		
		String R8=seq.substring(0,  12);
		if(b8.containsKey(R8)){return 12;}
		R8=seq.substring(0,  11);
		if(b8.containsKey(R8)){return 11;}
		R8=seq.substring(0,  10);
		if(b8.containsKey(R8)){return 10;}
		R8=seq.substring(0,  9);
		if(b8.containsKey(R8)){return 9;}
		R8=seq.substring(0, 8);
		if(b8.containsKey(R8)){return 8;}
		R8=seq.substring(0, 5);
		if(b8.containsKey(R8)){return 5;}
		
		
		
		R8=seq.substring(0,  12);
		if(b8m1.containsKey(R8)){return 12;}
		R8=seq.substring(0,  11);
		if(b8m1.containsKey(R8)){return 11;}
		R8=seq.substring(0,  10);
		if(b8m1.containsKey(R8)){return 10;}
		R8=seq.substring(0,  9);
		if(b8m1.containsKey(R8)){return 9;}
		R8=seq.substring(0, 8);
		if(b8m1.containsKey(R8)){return 8;}
		R8=seq.substring(0, 5);
		if(b8m1.containsKey(R8)){return 5;}
		
		R8=seq.substring(0,  12);
		if(b8m2.containsKey(R8)){return 12;}
		R8=seq.substring(0,  11);
		if(b8m2.containsKey(R8)){return 11;}
		R8=seq.substring(0,  10);
		if(b8m2.containsKey(R8)){return 10;}
		R8=seq.substring(0,  9);
		if(b8m2.containsKey(R8)){return 9;}
		R8=seq.substring(0, 8);
		if(b8m2.containsKey(R8)){return 8;}
		R8=seq.substring(0, 5);
		if(b8m2.containsKey(R8)){return 5;}
		
		return -1;
	}
	
	
	private String get(Map<String, String> barcodeNames, String r8) {
		String rtrn="NF_"+r8;
		
		if(barcodeNames.containsKey(r8)){rtrn=barcodeNames.get(r8);}
		
		return "["+rtrn+"]";
	}
	
	private static Map<String, String> parse(String string) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		Collection<Sequence> seqs= FastaFileIOImpl.readFromFile(string);
		//Map<String, String> mismatch1=new TreeMap<String, String>();
		
		for(Sequence seq: seqs){
			//Map<String, String> mismatch=mismatch(seq);
			//mismatch1.putAll(mismatch);
			//rtrn.putAll(mismatch);
			rtrn.put(seq.getSequenceBases(), seq.getName());
		}
		
		
		//mismatch 2
		/*for(String seq: mismatch1.keySet()){
			Map<String, String> mismatch=mismatch(seq, mismatch1.get(seq));
			rtrn.putAll(mismatch);
		}*/
		
		return rtrn;
	}
	
	private static Map<String, String> mismatch(String seq, String name, boolean withN) {
		Sequence s=new Sequence(name, seq);
		return mismatch(s, withN);
	}
	
	private static Map<String, String> mismatch(Sequence seq, boolean withN) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(int i=0; i<seq.getSequenceBases().length(); i++){
			Collection<String> newStrings=getStrings(seq, i, withN);
			for(String b: newStrings) {
				if(!b.equalsIgnoreCase(seq.getSequenceBases())) {
					rtrn.put(b, seq.getName());
				}
			}
		}
		return rtrn;
	}

	private static Collection<String> getStrings(Sequence seq, int i, boolean withN) {
		Collection<String> rtrn=new TreeSet<String>();
		String first=seq.getSequenceBases().substring(0, i);
		String second=seq.getSequenceBases().substring(i+1, seq.getSequenceBases().length());
		rtrn.add(first+"A"+second);
		rtrn.add(first+"C"+second);
		rtrn.add(first+"G"+second);
		rtrn.add(first+"T"+second);
		if(withN) {rtrn.add(first+"N"+second);}
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		if(args.length>10){
		File file1=new File(args[0]);
		File file2=new File(args[1]);
		//Map<String, String> barcodeNames=parse(args[2]);
		
		Map<String, String> b8=parse(args[2]);
		Map<String, String> b7=parse(args[3]);
		Map<String, String> b6=parse(args[4]);
		Map<String, String> b5=parse(args[5]);
		Map<String, String> b4=parse(args[6]);
		Map<String, String> b3=parse(args[7]);
		Map<String, String> b2=parse(args[8]);
		
		Map<String, String> pids=parse(args[9]);
		String save=args[10];
		
		new PooledCLIP(file1, file2, b8, b7, b6, b5, b4, b3, b2, pids, save);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=fastq (read 1) \n args[1]=fastq (read 2) \n args[2]=round 8 barcodes (fasta) \n args[3]=round 7 \n args[4]=round 6 \n args[5]=round 5 \n args[6]=round 4 \n args[7]=round 3 \n args[8]=round 3 \n args[9]=protein id barcodes \n args[10]=save";
	
}

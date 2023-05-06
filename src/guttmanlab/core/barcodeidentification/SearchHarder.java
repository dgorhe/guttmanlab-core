package guttmanlab.core.barcodeidentification;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class SearchHarder {

	
	
	public SearchHarder(File file1, File file2, Map<String, String> b8, String out1, String out2, int numRounds) throws IOException{
		FileWriter writer1=new FileWriter(out1);
		FileWriter writer2=new FileWriter(out2);
		
		FastqReader f1=new FastqReader(file1);
		FastqReader f2=new FastqReader(file2);
		
		Iterator<FastqRecord> iter1=f1.iterator();
		Iterator<FastqRecord> iter2=f2.iterator();
		
		Map<String, Map<String, String>> mismatchesByRound=mismatch(b8);
		Map<String, Map<String, String>> mismatchesByRound2=mismatch2(mismatchesByRound);
		
		System.err.println("mismatches made "+mismatchesByRound.size()+" "+mismatchesByRound.keySet());
		System.err.println("mismatches 2 made "+mismatchesByRound2.size()+" "+mismatchesByRound2.keySet()+" "+mismatchesByRound2.get("ROUND1").size());
		
		
		int counter=0;
		int fixed=0;
		while(iter2.hasNext()){
			FastqRecord record1=iter1.next();
			FastqRecord record2=iter2.next();
			
			Collection<String> barcodes=searchHarder(record2.getReadString(), b8); //TODO Cleanup this implementation
			
			
			Collection<String> list=new TreeSet<String>();
			for(String b: barcodes) {list.add(b.split("_")[0]);}
			
			
			//Mismatch1
			if(list.size()<numRounds) {
				Map<String, String> subset=filter(list, mismatchesByRound);
				Collection<String> temp=searchHarder(record2.getReadString(), subset);
				barcodes.addAll(temp);
				for(String b: temp) {
					list.add(b.split("_")[0]);
				}
				//System.err.println(counter+" "+temp.size()+" "+barcodes.size()+" "+barcodes);
			}
			
			//Mismatch2
			if(list.size()<numRounds) {
				Map<String, String> subset=filter(list, mismatchesByRound2);
				Collection<String> temp=searchHarder(record2.getReadString(), subset);
				barcodes.addAll(temp);
				for(String b: temp) {
					list.add(b.split("_")[0]);
				}
			}
			
			
			
			if(list.size()==numRounds) {
				fixed++;
			}
			
			addMissing(list, mismatchesByRound2, barcodes);
			
			counter++;
			if(counter%10000==0){
				System.err.println(counter+" "+fixed);
			}
			
			String name="@"+record2.getReadHeader();
			//int numPrevious=countPrev(record2.getReadHeader());
			name+="::"+updateName(record2.getReadHeader(), barcodes);
			
			write(record1, name, writer1);
			write(record2, name, writer2);
			
		}
		
		f1.close();
		f2.close();
		
		writer1.close();
		writer2.close();
	}
	
	
	private int countPrev(String readHeader) {
		int counter=0;
		String[] originalBarcodes=parse2(readHeader.split("::")[1]);
		for(int i=0; i<originalBarcodes.length; i++) {
			if(!originalBarcodes[i].equals("NOT_FOUND")) {counter++;}
		}
		return counter;
	}


	private String updateName(String readHeader, Collection<String> barcodes) {
		//Parse read 1
		// replace read name for both
		//Illuminaflowcell:: DPM/Bead:: NY
		
		/*String runName=readHeader.split("::")[0];
		String[] originalBarcodes=parse2(readHeader.split("::")[1]);
		String name="@"+runName+"::"+"["+originalBarcodes[0]+"]"+"["+originalBarcodes[1]+"]";*/
		
		String name="";
		
		
		//TODO Add new
		
		String[] list=new String[barcodes.size()];
		
		int pos=barcodes.size()-1;
		for(String b: barcodes) {
			list[pos]=b;
			pos--;
		}
		
		for(int i=0; i<list.length; i++) {
			String add=list[i];
			if(list[i].contains("_NF")) {add="NOT_FOUND";}
			name+="["+add+"]";
		}
		
		//System.err.println(readHeader+" "+name);		
	
		
		return name;
	}

	private String[] parse2(String list) {
		String[] barcode=list.split("\\[");
		String[] rtrn=new String[barcode.length-1];
		
		for(int i=1; i<barcode.length; i++) {
			rtrn[i-1]=barcode[i].replaceAll("\\]", "");
		}
		
		return rtrn;
		
	}

	private void write(FastqRecord record, String name, FileWriter writer) throws IOException {
		writer.write(name+"\n");
		writer.write(record.getReadString()+"\n");
		writer.write("+\n");
		writer.write(record.getBaseQualityString()+"\n");
	}

	
	private void addMissing(Collection<String> list, Map<String, Map<String, String>> b8m1, Collection<String> barcodes) {
		
		for(String round: b8m1.keySet()) {
			if(!list.contains(round)) {
				barcodes.add(round+"_NF");
			}
		}
	}
	

	private Map<String, String> filter(Collection<String> list, Map<String, Map<String, String>> b8m1) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String round: b8m1.keySet()) {
			if(!list.contains(round)) {
				rtrn.putAll(b8m1.get(round));
			}
		}
		
		return rtrn;
	}


	
	private Map<String, Map<String, String>> mismatch2(Map<String, Map<String, String>> map) {
		Map<String, Map<String, String>> rtrn=new TreeMap<String, Map<String, String>>();
		
		for(String round: map.keySet()) {
			Map<String, String> seqs=map.get(round);
			Map<String, String> temp=new TreeMap<String, String>();
			for(String seq: seqs.keySet()) {
				//String name=seqs.get(seq).replaceAll("_M1", "_M2");
				String name=seqs.get(seq);
				Map<String, String> mismatch=mismatch(seq, name);
				temp.putAll(mismatch);
			}
			rtrn.put(round, temp);	
		}
		
		return rtrn;
	}
	
	private Map<String, Map<String, String>> mismatch3(Map<String, Map<String, String>> map) {
		Map<String, Map<String, String>> rtrn=new TreeMap<String, Map<String, String>>();
		
		for(String round: map.keySet()) {
			Map<String, String> seqs=map.get(round);
			Map<String, String> temp=new TreeMap<String, String>();
			for(String seq: seqs.keySet()) {
				String name=seqs.get(seq).replaceAll("_M2", "_M3");
				Map<String, String> mismatch=mismatch(seq, name);
				temp.putAll(mismatch);
			}
			rtrn.put(round, temp);	
		}
		
		return rtrn;
	}


	private Map<String, Map<String, String>> mismatch(Map<String, String> b82) {
		Map<String, Map<String, String>> rtrn=new TreeMap<String, Map<String, String>>();
		
		for(String seq: b82.keySet()) {
			String round=b82.get(seq).split("_")[0];
			if(!rtrn.containsKey(round)) {rtrn.put(round, new TreeMap<String, String>());}
			Map<String, String> map=rtrn.get(round);
			//String name=b82.get(seq)+"_M1";
			
			String name=b82.get(seq);
			Map<String, String> mismatch=mismatch(seq, name);
			map.putAll(mismatch);
			rtrn.put(round, map);
		}
		
		return rtrn;
	}
	
	private static Map<String, String> mismatch(String seq, String name) {
		Sequence s=new Sequence(name, seq);
		return mismatch(s);
	}
	
	private static Map<String, String> mismatch(Sequence seq) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(int i=0; i<seq.getSequenceBases().length(); i++){
			Collection<String> newStrings=getStrings(seq, i);
			for(String b: newStrings) {
				if(!b.equalsIgnoreCase(seq.getSequenceBases())) {
					rtrn.put(b, seq.getName());
				}
			}
		}
		return rtrn;
	}
	
	private static Collection<String> getStrings(Sequence seq, int i) {
		Collection<String> rtrn=new TreeSet<String>();
		String first=seq.getSequenceBases().substring(0, i);
		String second=seq.getSequenceBases().substring(i+1, seq.getSequenceBases().length());
		rtrn.add(first+"A"+second);
		rtrn.add(first+"C"+second);
		rtrn.add(first+"G"+second);
		rtrn.add(first+"T"+second);
		//rtrn.add(first+"N"+second);
		return rtrn;
	}
	
	/*private String searchHarder(String seq, Map<String, String> b2) {
		Collection<Integer> sizes=new TreeSet<Integer>();
			for(String b: b2.keySet()){sizes.add(b.length());}
			
			for(int i=0; i<seq.length(); i++){
				Collection<String> kmers=getKmer(seq, i, sizes);
				for(String kmer: kmers){
					if(b2.containsKey(kmer)){
						String name=b2.get(kmer);
						return name;
					}
				}
			}
			return "NF";
		}*/
	
	
	private Collection<String> searchHarder(String seq, Map<String, String> b2) {
		Collection<String> rtrn=new TreeSet<String>();
		Collection<Integer> sizes=new TreeSet<Integer>();
			for(String b: b2.keySet()){sizes.add(b.length());}
			
			for(int i=0; i<seq.length(); i++){
				Collection<String> kmers=getKmer(seq, i, sizes);
				for(String kmer: kmers){
					if(b2.containsKey(kmer)){
						String name=b2.get(kmer);
						rtrn.add(name);
					}
				}
			}
			return rtrn;
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
	
	private static Map<String, String> parse(String string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		//Collection<Sequence> seqs= FastaFileIOImpl.readFromFile(string);
		
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			rtrn.put(tokens[2], tokens[1]);
			//rtrn.put(Sequence.reverseComplement(tokens[2]), tokens[1]+"_RC");
		}
		
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException {
		if(args.length>5) {
			File file1=new File(args[0]);
			File file2=new File(args[1]);
			Map<String, String> b8=parse(args[2]);
			String out1=args[3];
			String out2=args[4];
			int numBarcodes=Integer.parseInt(args[5]);
			new SearchHarder(file1, file2, b8, out1, out2, numBarcodes);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=read1 (fq) \n args[1]=read2 (fq) \n args[2]=barcodes \n args[3]=out1 \n args[4]=out2 \n args[5]=number of barcode rounds";
}

package guttmanlab.core.sars.leader;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

public class ParseVariants {
	
	private String proteinName="NSP1_";

	public ParseVariants (String file, String save, String proteinName) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file);
		this.proteinName=proteinName;
		
		List<String> positions=getPositions(lines);
		List<String> dates=getDates();
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(positions, dates);
			
		for(String date: dates) {
			computeByDate(lines, date, mwh);
		}
			
		mwh.write(save);
	}
	
	
	private Collection<String> getLocations(List<String> lines) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String line: lines) {
			String loc=line.split("\t")[8];
			rtrn.add(loc);
		}
		
		for(String loc: rtrn) {System.out.println(loc);}
		
		return rtrn;
	}






	private List<String> getDates() {
		List<String> rtrn=new ArrayList<String>();
		
		rtrn.add("2020-03");
		rtrn.add("2020-04");
		rtrn.add("2020-05");
		rtrn.add("2020-06");
		rtrn.add("2020-07");
		rtrn.add("2020-08");
		rtrn.add("2020-09");
		rtrn.add("2020-10");
		rtrn.add("2020-11");
		rtrn.add("2020-12");
		
		rtrn.add("2021-01");
		rtrn.add("2021-02");
		rtrn.add("2021-03");
		rtrn.add("2021-04");
		rtrn.add("2021-05");
		rtrn.add("2021-06");
		rtrn.add("2021-07");
		rtrn.add("2021-08");
		rtrn.add("all");
		
		return rtrn;
	}






	private List<String> getPositions(List<String> lines) {
		List<String> rtrn=new ArrayList<String>();
		rtrn.add("total");
		rtrn.add("allNSP1");
		for(String line: lines) {
			String[] tokens=line.split("\t");
			Collection<String> nsp1Mutants=parseNSP1(tokens[5]);
			for(String mut: nsp1Mutants) {
				if(!rtrn.contains(mut)) {rtrn.add(mut);}
			}
		}
		return rtrn;		
	}






	private Map<String, Integer> computeByDate(List<String> lines, String dateStart, MatrixWithHeaders mwh) {
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		
		int total=0;
		int allNSP1=0;
		for(String line: lines) {
			String[] tokens=line.split("\t");
			String date=(tokens[7]);
			if(date.startsWith(dateStart) || dateStart.equals("all")) {
				Collection<String> nsp1Mutants=parseNSP1(tokens[5]);
				for(String mut: nsp1Mutants) {
					int count=0;
					if(counts.containsKey(mut)) {count=counts.get(mut);}
					count++;
					counts.put(mut, count);
				}
				if(!nsp1Mutants.isEmpty()) {allNSP1++;}
				total++;
			}
		}
		
		mwh.set("total", dateStart, total);
		mwh.set("allNSP1", dateStart, allNSP1);
		
		for(String mutation: counts.keySet()) {
			double fraction=(double)counts.get(mutation)/(double)total;
			mwh.set(mutation, dateStart, counts.get(mutation));
		}
		
		return counts;
	}






	private SAMRecord makeSAM(String[] tokens, SAMFileHeader header, Sequence reference, Collection<String> nsp1Mutants) {
		SAMRecord record1=new SAMRecord(header);
		record1.setReferenceName("NSP1");
		record1.setAlignmentStart(1);
		record1.setReadName(tokens[0]);
		String trimmed=getSeq(reference, nsp1Mutants);
		record1.setReadString(trimmed);
		record1.setBaseQualityString(makeArtificialQuality(trimmed.length()));
		record1.setCigarString(getCigar(nsp1Mutants, reference));
		return record1;
	}
	
	private String getSeq(Sequence reference, Collection<String> nsp1Mutants) {
		int count=0;
		for(String mut: nsp1Mutants) {
			if(mut.startsWith(proteinName)){
				String sub=mut.split("_")[1];
				if (sub.contains("del")){
					count++;
				}
			}
		}
		
		int size=reference.getLength()-count;
		
		String rtrn="";
		for(int i=0; i<size; i++) {rtrn+="N";}
		
		return rtrn;
	}

	
	private String getCigar(Collection<String> nsp1Mutants, Sequence reference) {
		//TreeSet<Integer> order=new TreeSet<Integer>();
		char[] chars=new char[reference.getLength()];
		
		for(int i=0; i<chars.length; i++) {chars[i]='M';}
		
		for(String mut: nsp1Mutants) {
			//System.err.println(mut);
			if(mut.startsWith(proteinName)){
				String sub=mut.split("_")[1];
				if (sub.contains("del")){
					int pos=Integer.parseInt(sub.substring(1, sub.length()-3));
					//order.add(pos);
					chars[pos]='D';
				}
				
			}
		}
		
		int start=0;
		String rtrn="";
		for(int i=0; i<chars.length; i++) {
			if(chars[i]=='D') {
				rtrn+=(i-start)+"M";
				int num=getDRun(chars, i);
				rtrn+=num+"D";
				start=i+num;
				i=i+num;
			}
		}
		rtrn+=(chars.length-start)+"M";
				
		//System.err.println(rtrn);
		
		return rtrn;
	}
	
	
	private Collection<SingleInterval> getDeletionDomains(Collection<String> nsp1Mutants, Sequence reference) {
		//TreeSet<Integer> order=new TreeSet<Integer>();
		char[] chars=new char[reference.getLength()];
		
		for(int i=0; i<chars.length; i++) {chars[i]='M';}
		
		for(String mut: nsp1Mutants) {
			//System.err.println(mut);
			if(mut.startsWith(proteinName)){
				String sub=mut.split("_")[1];
				if (sub.contains("del")){
					int pos=Integer.parseInt(sub.substring(1, sub.length()-3));
					//order.add(pos);
					chars[pos]='D';
				}
				
			}
		}
		
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		for(int i=0; i<chars.length; i++) {
			if(chars[i]=='D') {
				int start=i;
				int num=getDRun(chars, i);
				int end=start+num;
				SingleInterval temp=new SingleInterval("NSP1", start, end);
				rtrn.add(temp);
			
				i=i+num;
			}
		}
		
		return rtrn;
	}
	
	
	

	private int getDRun(char[] chars, int start) {
		for(int i=start; i<chars.length; i++) {
			if(chars[i]!='D' || i==chars.length-1) {
				return i-start;
			}
		}
		return 0;
	}

	private String makeArtificialQuality(int readLength2) {
		String rtrn="";
		for(int i=0; i<readLength2; i++){
			rtrn+="I";
		}
		return rtrn;
	}
	
	
	private SAMFileHeader getHeader(Sequence reference) {
		SAMFileHeader header = new SAMFileHeader();
		
	       int size = reference.getLength();
	       SAMSequenceRecord seq = new SAMSequenceRecord("NSP1", size);
	       header.addSequence(seq);
	       
	        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
	                
	        return header;
	}
	
	
	private static void write(String string, Map<String, Integer> counts, int total) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		Map<Integer, Integer> substitutions=getSubs(counts);
		Map<Integer, Integer> deletions=getDeletions(counts);
		
		Collection<Integer> allPos=new TreeSet<Integer>();
		allPos.addAll(substitutions.keySet());
		allPos.addAll(deletions.keySet());
		
		for(Integer pos: allPos) {
			double sub=get(substitutions, pos)/(double)total;
			double del=get(deletions, pos)/(double)total;
			writer.write(pos+"\t"+sub+"\t"+del+"\n");
		}
		
		
		/*for(String mut: counts.keySet()) {
			int pos=getPos(mut);
			//String type=getType(mut);
			double ratio=(double)counts.get(mut)/(double)total;
			writer.write(mut+"\t"+pos+"\t"+counts.get(mut)+"\t"+ratio+"\n");
		}*/
		
		writer.close();
	}


	private static int get(Map<Integer, Integer> substitutions, Integer pos) {
		int count=0;
		if(substitutions.containsKey(pos)) {count=substitutions.get(pos);}
		return count;
	}



	private static Map<Integer, Integer> getDeletions(Map<String, Integer> counts) {
		Map<Integer, Integer> rtrn=new TreeMap<Integer, Integer>();
		
		for(String mut: counts.keySet()) {
			String sub=mut.split("_")[1];
				if (sub.contains("del")){
					int pos=Integer.parseInt(sub.substring(1, sub.length()-3));
					int count=0;
					if(rtrn.containsKey(pos)) {count=rtrn.get(pos);}
					count+=counts.get(mut);
					rtrn.put(pos, count);
				}
		}
		
		
		return rtrn;
	}
	
	
	private static Map<Integer, Integer> getSubs(Map<String, Integer> counts) {
		Map<Integer, Integer> rtrn=new TreeMap<Integer, Integer>();
		
		for(String mut: counts.keySet()) {
			String sub=mut.split("_")[1];
				if (!sub.contains("del") && !sub.contains("ins") && !sub.contains("stop")){
					int pos=Integer.parseInt(sub.substring(1, sub.length()-1));
					int count=0;
					if(rtrn.containsKey(pos)) {count=rtrn.get(pos);}
					count+=counts.get(mut);
					rtrn.put(pos, count);
				}
		}
		
		
		return rtrn;
	}



	private int getPos(String mut) {
		int pos=-1;
		if(mut.startsWith(proteinName)){
			String sub=mut.split("_")[1];
			if (sub.contains("del")){
				pos=Integer.parseInt(sub.substring(1, sub.length()-3));
			}
			else if(sub.contains("ins")) {
				
			}
			else if(sub.contains("stop")) {
				pos=Integer.parseInt(sub.substring(1, sub.length()-4));
			}
			
			else {
				pos=Integer.parseInt(sub.substring(1, sub.length()-1));
			}	
		}
		return pos;
	}



	private Collection<String> parseNSP1(String string) {
		Collection<String> rtrn=new TreeSet<String>();
		
		if(string.equals("()") || string.isEmpty()){return rtrn;}
		else if(!string.contains("(")) {return rtrn;}
		
		String trimmed=string.split("\\(")[1].split("\\)")[0];
			
		
		
		
		String[] tokens=trimmed.split(",");
		
		for(int i=0; i<tokens.length; i++) {
			if(tokens[i].startsWith(proteinName)) {rtrn.add(tokens[i]);}
		}
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException {
		new ParseVariants(args[0], args[1], args[2]);
	}
	
}

package guttmanlab.core.sars.leader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

public class PullFastaFileRegion {

	public PullFastaFileRegion(File fasta, int start, int end, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		//Collection<Sequence> seqs=FastaFileIOImpl.readFromFile(fasta);
		
		FastaFileIOImpl f=new FastaFileIOImpl(fasta);
		
		/*SAMFileHeader header=getHeader();
		header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		SAMFileWriter alignmentWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(save+".bam"));
		*/
		
		int counter=0;
		while(f.hasNext()) {
			Sequence seq=f.next();
			String sub=seq.getSubSequence(start, end);
			counter++;
			//SAMRecord r=makeSAM(sub, seq.getName(), "ORF1a", header);
			//alignmentWriter.addAlignment(r);
			writer.write(">"+seq.getName()+"\n"+sub+"\n");
			if(counter%10000==0) {System.err.println(counter+" "+seq.getName());}
		}
		writer.close();
		//alignmentWriter.close();
	}
	
	
	private SAMFileHeader getHeader() {
		SAMFileHeader header = new SAMFileHeader();
		
	       int size = 13202;
	       SAMSequenceRecord seq = new SAMSequenceRecord("ORF1a", size);
	       header.addSequence(seq);
	       
	        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
	                
	        return header;
	}


	private SAMRecord makeSAM(String seq, String name, String geneName, SAMFileHeader header) {
		SAMRecord record1=new SAMRecord(header);
		record1.setReferenceName(geneName);
		record1.setAlignmentStart(1);
		record1.setReadName(name);
		String trimmed=trim(seq);
		record1.setReadString(trimmed);
		record1.setBaseQualityString(makeArtificialQuality(trimmed.length()));
		record1.setCigarString(getCigar(seq));
		return record1;
	}
	
	private String trim(String seq) {
		String rtrn="";
		char[] chars=seq.toCharArray();
		
		for(int i=0; i<chars.length; i++) {
			if(chars[i]!='-') {rtrn+=chars[i];}
		}
		
		//System.err.println(seq+"\n"+rtrn);
		return rtrn;
	}


	private String getCigar(String seq) {
		char[] chars=seq.toCharArray();
		LinkedList<Integer> posOfIndel=new LinkedList<Integer>();
		for(int i=0; i<chars.length; i++) {
			if(chars[i]=='-') {
				posOfIndel.add(i);
			}
		}
		
		String rtrn="";
		int start=0;
		//Iterator<Integer> iter=posOfIndel.iterator();
		
		if(posOfIndel.isEmpty()) {rtrn=seq.length()+"M";}
		
		else {
			while(!posOfIndel.isEmpty()) {
				int pos=posOfIndel.pop();
				//int pos=iter.next();
				rtrn+=(pos-start)+"M";
				int numD=getD(pos, posOfIndel);
				start=pos+numD;
				rtrn+=numD+"D";
			}
			if(seq.length()-start>0){
				rtrn+=(seq.length()-start)+"M";
			}
		}
		
		//System.err.println(rtrn);
		
		return rtrn;
	}


	private int getD(int prev, LinkedList<Integer> posOfIndel) {
		int counter=1;
		if(!posOfIndel.isEmpty()) {
			Integer next=posOfIndel.peek();
			while(next==(prev+1)) {
				prev=posOfIndel.pop();
				if(!posOfIndel.isEmpty()) {next=posOfIndel.peek();}
				counter++;
			}
		}
		return counter;
	}


	private String makeArtificialQuality(int readLength2) {
		String rtrn="";
		for(int i=0; i<readLength2; i++){
			rtrn+="I";
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>3) {
		File fasta=new File(args[0]);
		int start=Integer.parseInt(args[2]);
		int end=Integer.parseInt(args[3]);
		String save=args[1];
		new PullFastaFileRegion(fasta, start, end, save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=fasta \n args[1]=save \n args[2]=start \n args[3]=end";
	
}

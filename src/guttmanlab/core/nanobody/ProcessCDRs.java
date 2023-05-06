package guttmanlab.core.nanobody;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;

import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class ProcessCDRs {

	public static void processCDRs(File fastq, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		FastqReader f=new FastqReader(fastq);
		Iterator<FastqRecord> iter=f.iterator();
		
		int counter=0;
		while(iter.hasNext()) {
			FastqRecord seq=iter.next();
			String seqBases=seq.getReadString();
			String CDR3=getCDR3(seqBases);
			
			if(CDR3!=null) {
				String CDR2=getCDR2(seqBases, CDR3);
				if(CDR2!=null) {
					writer.write(seq.getReadHeader().split(" ")[0]+"\t"+CDR2+"\t"+CDR3+"\t"+translate(CDR2,CDR3)+"\n");
				}
			}
			
			
			
			/*if(CDR3!=null && CDR2!=null) {
				writer.write(seq.getReadHeader()+"\t"+CDR2+"\t"+CDR3+"\n");
			}
			else {
				//System.err.println("Full seq: "+seqBases);
			}*/
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		f.close();
		writer.close();
	}
	
	private static String translate(String CDR2, String CDR3) {
		String rCDR2=Sequence.reverseComplement(CDR2);
		String rCDR3=Sequence.reverseComplement(CDR3);
		
		String aa2=Sequence.translate(rCDR2);
		String aa3=Sequence.translate(rCDR3);
		return aa2+aa3;
	}

	private static String getCDR2(String seqBases, String cDR3) {
		//String barcode=getVal1_2(seqBases);
		
		
		int commonPos=seqBases.indexOf(cDR3);
		
		if(commonPos<0 || (commonPos+cDR3.length()+123+15)>seqBases.length()) {return null;}
		
		//System.err.println(seqBases+" "+commonPos+" "+(commonPos+cDR3.length())+" "+ (commonPos+cDR3.length()+15)+" "+seqBases.substring(commonPos+cDR3.length(), commonPos+cDR3.length()+15));
		
		return seqBases.substring(commonPos+cDR3.length()+123, commonPos+cDR3.length()+15+123);
	}

	

	private static String getVal1_2(String seqBases) {
		String commonSeq="AGAAGGC";
		String subseq=seqBases.substring(seqBases.length()-35, seqBases.length());
		
		int commonPos=subseq.indexOf(commonSeq);
		
		if(commonPos<0) {return null;}
		return seqBases.substring(commonPos+commonSeq.length(), commonPos+commonSeq.length()+15);
	}

	private static String getCDR3(String seqBases) {
		String val1=getVal1(seqBases);
		String val2=getVal2(seqBases);
		String barcode=consensus(val1, val2);
		return barcode;
	}

	private static String consensus(String val1, String val2) {
		if(val1==null) {return val2;}
		if(val2==null) {return val1;}
		if(val1.equals(val2)) {return val1;}
		//System.err.println(val1+" "+ val2);
		return val2;
	}

	private static String getVal2(String seqBases) {
		String commonSeq="AATAGTC";
		int commonPos=seqBases.indexOf(commonSeq);
		if(commonPos<0 || (commonPos+commonSeq.length()+30)>=seqBases.length()) {return null;}
		return seqBases.substring(commonPos+commonSeq.length(), commonPos+commonSeq.length()+30);
	}

	private static String getVal1(String seqBases) {
		String commonSeq="CGCAGCA";
		int commonPos=seqBases.indexOf(commonSeq);
		if((commonPos-30)<0) {return null;}
		return seqBases.substring(commonPos-30, commonPos);
	}

	public static void main(String[] args) throws IOException {
		File fastq=new File(args[0]);
		String save=args[1];
		processCDRs(fastq, save);
	}
	
}

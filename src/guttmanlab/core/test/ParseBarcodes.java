package guttmanlab.core.test;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class ParseBarcodes {
	
	public ParseBarcodes(FastqReader reader, Map<String, String> primersToName, Map<String, String> barcodeToName){
		Iterator<FastqRecord> iter1=reader.iterator();
		
		/*int matchesR1=0;
		int totalR1=0;
		
		int matchesR2=0;
		int totalR2=0;
		
		int matchesR3=0;
		int totalR3=0;
		
		int allMatch=0;
		int allTotal=0;*/
		
		int count=0;
		int total=0;
		
		while(iter1.hasNext()){
			FastqRecord record=iter1.next();
			String primer=record.getReadString().substring(0,15);
			
			String R1=record.getReadString().substring(22,30);
			/*String R2=record.getReadString().substring(37,45);
			String R3=record.getReadString().substring(52,60);
			String R4=record.getReadString().substring(67,75);*/
			
			String primerName=primersToName.get(primer);
			String R1Name=barcodeToName.get(R1);
			/*String R2Name=barcodeToName.get(R2);
			String R3Name=barcodeToName.get(R3);
			String R4Name=barcodeToName.get(R4);*/
			/*String R5Name=barcodeToName.get(R5);
			String R6Name=barcodeToName.get(R6);
			String R6Name=barcodeToName.get(R7);*/
			
			if(has(primerName, R1Name)){
				System.out.println(record.getReadHeader()+"\t"+record.getReadString()+"\t"+primer+"\t"+R1+"\t"+primerName+"\t"+R1Name);
				if(matches(primerName, R1Name)){
					count++;
				}
				total++;
			}
			
			
			/*if(has(primerName, R2Name)){
				if(matchesR2(primerName, R2Name)){
					matchesR2++;
				}
				else{
					match=false;
				}
				totalR2++;
			}
			
			if(has(primerName, R3Name)){
				if(matchesR3(primerName, R3Name)){
					matchesR3++;
				}
				else{
					match=false;
				}
				totalR3++;
			}
			
			if(has(primerName, R1Name, R2Name, R3Name)){
				if(match){allMatch++;}
				else{
					System.out.println(primerName+"\t"+R1Name+"\t"+R2Name+"\t"+R3Name+"\t"+R4Name);
				}
				allTotal++;
			}*/
			//System.out.println(primer+" "+primerName+"\t"+R1Name+"\t"+R2Name+"\t"+R3Name+"\t"+R4Name);
			
		}
		
		double ratio = (double)count/(double)total;
		System.err.println(count+" "+total+" "+ratio);
		
		
		/*double ratioR1=(double)matchesR1/(double)totalR1;
		double ratioR2=(double)matchesR2/(double)totalR2;
		double ratioR3=(double)matchesR3/(double)totalR3;
		double allRatio=(double)allMatch/(double)allTotal;
		System.err.println(ratioR1+" "+ratioR2+" "+ratioR3+ " "+allRatio);*/
		
		reader.close();
	}

	
	private boolean has(String primerName, String r1Name, String r2Name, String r3Name) {
		return primerName!=null && r1Name!=null && r2Name!=null && r3Name!=null;
	}


	private boolean matchesR3(String primerName, String r2Name) {
		if((primerName.equals("P1") || primerName.equals("P2") || primerName.equals("P3") || primerName.equals("P4")) && (r2Name.equals("A1") || r2Name.equals("A2")|| r2Name.equals("A3") || r2Name.equals("A4"))){return true;}
		if((primerName.equals("P5") || primerName.equals("P6") || primerName.equals("P7") || primerName.equals("P8")) && (r2Name.equals("A5") || r2Name.equals("A6")|| r2Name.equals("A7") || r2Name.equals("A8"))){return true;}
		if((primerName.equals("P9") || primerName.equals("P10") || primerName.equals("P11") || primerName.equals("P12")) && (r2Name.equals("A9") || r2Name.equals("A10")|| r2Name.equals("A11") || r2Name.equals("A12"))){return true;}
		
		
		return false;
	}
	
	private boolean matchesR2(String primerName, String r2Name) {
		if((primerName.equals("P1") || primerName.equals("P2")) && (r2Name.equals("B1") || r2Name.equals("B2"))){return true;}
		if((primerName.equals("P3") || primerName.equals("P4")) && (r2Name.equals("B3") || r2Name.equals("B4"))){return true;}
		if((primerName.equals("P5") || primerName.equals("P6")) && (r2Name.equals("B5") || r2Name.equals("B6"))){return true;}
		if((primerName.equals("P7") || primerName.equals("P8")) && (r2Name.equals("B7") || r2Name.equals("B8"))){return true;}
		if((primerName.equals("P9") || primerName.equals("P10")) && (r2Name.equals("B9") || r2Name.equals("B10"))){return true;}
		if((primerName.equals("P11") || primerName.equals("P12")) && (r2Name.equals("B11") || r2Name.equals("B12"))){return true;}
		
		
		return false;
	}



	private boolean has(String primerName, String r1Name) {
		return primerName!=null && r1Name!=null;
	}



	private boolean matches(String primerName, String r1Name) {
		return primerName.replace("P", "").equals(r1Name.replace("A", ""));
	}

	

	private static Map<String, String> parseBarcodes(String file) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		List<String> lines=BEDFileIO.loadLines(file);
		
		for(String line: lines){
			
			String name=line.split("\t")[0];
			String seq=line.split("\t")[1];
			rtrn.put(seq, name);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		FastqReader reader1=new FastqReader(new File(args[0]));
		
		
		new ParseBarcodes(reader1, parseBarcodes(args[1]), parseBarcodes(args[2]));
		
		
		
		
	}
	
	
}

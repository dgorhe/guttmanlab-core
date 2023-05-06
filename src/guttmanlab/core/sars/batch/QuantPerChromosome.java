package guttmanlab.core.sars.batch;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.datastructures.MatrixWithHeaders;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class QuantPerChromosome {

	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		List<String> columns=new ArrayList<String>();
		MatrixWithHeaders matrix=null;
		for(int i=0; i< files.length; i++){
			if(files[i].getName().endsWith(".bam")){
				System.err.println(files[i].getName());
				columns.add(files[i].getName());
			}
		}
		
		
		for(int i=0; i< files.length; i++){
			if(files[i].getName().endsWith(".bam")){
				String column=files[i].getName();
				Map<String, Integer> countsPerChr=getCountsPerChr(files[i]);
				if(matrix==null){List<String> rows=getRows(countsPerChr); matrix=new MatrixWithHeaders(rows, columns);}
				for(String row: countsPerChr.keySet()){matrix.set(row, column, countsPerChr.get(row));}
			}
		}
		
		matrix.write(args[1]);
	}

	private static List<String> getRows(Map<String, Integer> countsPerChr) {
		List<String> rtrn=new ArrayList<String>();
		rtrn.addAll(countsPerChr.keySet());
		return rtrn;
	}

	private static Map<String, Integer> getCountsPerChr(File bam) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag()){
				String chr=read.getReferenceName();
				int score=0;
				if(rtrn.containsKey(chr)){score=rtrn.get(chr);}
				score++;
				rtrn.put(chr, score);
			}
			
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
	
		return rtrn;
	}
	
}

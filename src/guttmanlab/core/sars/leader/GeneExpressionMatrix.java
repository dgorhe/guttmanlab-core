package guttmanlab.core.sars.leader;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.datastructures.MatrixWithHeaders;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class GeneExpressionMatrix {

	
	public GeneExpressionMatrix(File[] bamFiles, String save) throws IOException {
		
		Map<String, Double>[] scores=new Map[bamFiles.length];
		
		for(int i=0; i<bamFiles.length; i++) {
			System.err.println(bamFiles[i].getName());
			scores[i]=getScores(bamFiles[i]);
		}
		
		MatrixWithHeaders matrix=makeMatrix(scores, bamFiles);
		
		matrix.write(save);
		
	}

	private MatrixWithHeaders makeMatrix(Map<String, Double>[] scores, File[] bamFiles) {
		List<String> columns=getNames(bamFiles);
		List<String> rows=getNames(scores);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(int i=0; i<bamFiles.length; i++) {
			Map<String, Double> map=scores[i];
			add(rtrn, map, bamFiles[i].getName());
		}
		
		return rtrn;
	}

	private List<String> getNames(Map<String, Double>[] scores) {
		TreeSet<String> temp=new TreeSet<String>();
		
		for(int i=0; i<scores.length; i++) {temp.addAll(scores[i].keySet());}
		
		List<String> rtrn= new ArrayList<String>();
		rtrn.addAll(temp);
		return rtrn;
	}

	private List<String> getNames(File[] bamFiles) {
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<bamFiles.length; i++) {
			rtrn.add(bamFiles[i].getName());
		}
		
		return rtrn;
	}

	private void add(MatrixWithHeaders rtrn, Map<String, Double> map, String columnName) {
		for(String row: map.keySet()) {
			rtrn.set(row, columnName, map.get(row));
		}
		
	}

	private Map<String, Double> getScores(File file) {
		SAMFileReader inputReader= new SAMFileReader(file);
		Map<String, Double> geneCount=new TreeMap<String, Double>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				String gene=read.getReferenceName();
				
				double score=0;
				if(geneCount.containsKey(gene)){score=geneCount.get(gene);}
				score++;
				geneCount.put(gene, score);
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		return geneCount;
	}
	
	
	
	public static void main(String[] args)throws IOException{
		File[] bamFiles=new File(args[0]).listFiles();
		String save=args[1];
		new GeneExpressionMatrix(bamFiles, save);
	}
	
	
}

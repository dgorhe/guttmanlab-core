package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.rnasprite.Kmer;

public class ProteinKmers {

	public ProteinKmers(BarcodingDataStreaming data, String save) throws IOException {
		Map<String, Integer> proteinCounts=score(data);
		
		List<String> rows=makeList(proteinCounts.keySet());
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, rows);
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Collection<String> proteins=c.getProteinSet();
			for(String p1:proteins) {
				for(String p2: proteins) {
					mwh.incrementCount(p1, p2);
				}
			}
			
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		data.close();
	
		mwh=norm(mwh, proteinCounts, counter);
		
		mwh.write(save);
		
		//write(save, kmerCounts, proteinCounts, counter);
	}

	private MatrixWithHeaders norm(MatrixWithHeaders mwh, Map<String, Integer> proteinCounts, int total) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		for(String row: mwh.getRowNames()) {
			for(String column: mwh.getColumnNames()) {
				double score=mwh.get(row, column);
				double freq1=(double)proteinCounts.get(row)/(double)total;
				double freq2=(double)proteinCounts.get(column)/(double)total;
				double exepected=freq1*freq2*total;
				rtrn.set(row, column, score/exepected);
			}
		}
			
		return rtrn;
	}

	private Map<String, Integer> score(BarcodingDataStreaming data) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			
			for(String p: c.getProteinSet()) {
				int count=0;
				if(rtrn.containsKey(p)) {count=rtrn.get(p);}
				count++;
				rtrn.put(p, count);
			}
			
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		data.close();
		
		return rtrn;
	}

	private List<String> makeList(Set<String> keySet) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String k: keySet) {rtrn.add(k);}
		
		return rtrn;
	}

	private Map<Kmer, Integer> cumulativeTotals(Map<Kmer, Integer> kmerCounts) {
		Map<Kmer, Integer> rtrn=new TreeMap<Kmer, Integer>();
		
		for(Kmer k: kmerCounts.keySet()) {
			int sum=kmerCounts.get(k);
		
			for(Kmer sub: kmerCounts.keySet()) {
				if(sub(sub, k)) {
					//System.err.println(k+" "+sub);
					int count=kmerCounts.get(sub);
					count=count+sum;
					rtrn.put(sub, count);
				}
			}
		}
		
		return rtrn;
	}

	private boolean sub(Kmer sub, Kmer k) {
		if(sub.getSize()>k.getSize()) {return false;}
		
		boolean missing=false;
		for(String r: sub.getRegions()) {
			if(!k.containsRegion(r)) {missing=true;}
		}
		
		return !missing;
		
	}

	private Map<String, Integer> score(Map<Kmer, Integer> kmerCounts) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		for(Kmer k: kmerCounts.keySet()) {
			int count=kmerCounts.get(k);
			for(String p: k.getRegions()) {
				int sum=0;
				if(rtrn.containsKey(p)) {sum=rtrn.get(p);}
				sum+=count;
				rtrn.put(p, sum);
			}
		}
		
		return rtrn;
	}

	private Kmer makeProteinKmer(Cluster c) {
		Kmer k=new Kmer();
		
		k.addRegions(c.getProteinSet());
		
		return k;
	}

	private void write(String save, Map<Kmer, Integer> counts, Map<String, Integer> proteinCounts, int total) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Kmer k: counts.keySet()) {
			writer.write(k.toString()+"\t"+counts.get(k));
			
			for(String p: k.getRegions()) {
				int count=proteinCounts.get(p);
				writer.write("\t"+count);
			}
			
			writer.write("\t"+total+"\n");		
		}
		
		writer.close();
	}

	private void add(Kmer k, Map<Kmer, Integer> counts) {
		int count=0;
		if(counts.containsKey(k)) {count=counts.get(k);}
		count++;
		counts.put(k, count);
	}
	
	public static void main(String[] args) throws IOException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		new ProteinKmers(data, save);
	}
	
}

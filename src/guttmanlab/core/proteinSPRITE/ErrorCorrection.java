package guttmanlab.core.proteinSPRITE;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.Kmer;

public class ErrorCorrection {

	public ErrorCorrection(Collection<Kmer> referenceBarcodes) {
		
		Map<Kmer, Integer> subToFull=subToFull(referenceBarcodes);
		//Map<String, Collection<String>> subToFull2=subToFull(unmatched);
		
		write(subToFull);
		
	}

	private void write(Map<Kmer, Integer> subToFull) {
		for(Kmer sub: subToFull.keySet()) {
			System.out.println(sub+" "+subToFull.get(sub));
		}
		
	}

	private Map<Kmer, Integer> subToFull(Collection<Kmer> referenceBarcodes) {
		Map<Kmer, Integer> rtrn=new TreeMap<Kmer, Integer>();
		
		int counter=0;
		for(Kmer barcode: referenceBarcodes) {
			Collection<Kmer> subs=barcode.subtractOne();
			for(Kmer sub: subs) {
				int count=0;
				if(rtrn.containsKey(sub)) {count=rtrn.get(sub);}
				count++;
				//if(count>1) {System.err.println(sub+" "+count);}
				rtrn.put(sub, count);
			}
			counter++;
			if(counter%10000==0) {System.err.println(counter);}
		}
		
		
		return rtrn;
	}
	
	private String makeString(String[] b, int skip) {
		String rtrn="";
		
		int start=0;
		for(int i=0; i<b.length; i++) {
			if(i!=skip) {
				if(start>0) {rtrn+=".";}
				rtrn+=b[i];
				start++;
			}
		}
		return rtrn;
	}

	private Collection<Kmer> removeOne(Kmer k) {
		return k.enumerateSubK(k.getSize()-1);
	}
	
	public static void main(String[] args) throws IOException {
		Collection<Kmer> barcodes=parse(args[0]);
		new ErrorCorrection(barcodes);
	}

	private static Collection<Kmer> parse(String string) throws IOException {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		List<String> lines=BEDFileIO.loadLines(string);
		
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			String[] b=tokens[0].split("\\.");
			Kmer k=new Kmer();
			for(String n: b) {k.addRegion(n);}
			rtrn.add(k);
		}
		
		return rtrn;
	}
	
}

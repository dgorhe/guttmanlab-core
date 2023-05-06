package guttmanlab.core.rnasprite.hubs;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.rnasprite.Kmer;

public class MissingKmers {

	
	private static void write(String string, Collection<Kmer> unscored) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(Kmer k: unscored) {writer.write(k+"\tmissing"+"\n");}
		
		
		writer.close();
	}
	
	
	private static Collection<Kmer> parse(String string) throws IOException {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines) {
			rtrn.add(new Kmer(line.split("\t")[0]));
		}
		
		return rtrn;
	}
	
	private static Collection<Kmer> diff(Collection<Kmer> allKmers, Collection<Kmer> scoredKmers) {
		Collection<Kmer> rtrn=new TreeSet<Kmer>();
		
		for(Kmer k: allKmers) {
			if(!scoredKmers.contains(k)) {rtrn.add(k);}
		}
		
		
		return rtrn;
	}

	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
		Collection<Kmer> allKmers=parse(args[0]);
		Collection<Kmer> scoredKmers=parse(args[1]);
		
		Collection<Kmer> unscored=diff(allKmers, scoredKmers);
		write(args[2], unscored);
		}
		else {System.err.println(usage);}
	}

	
	


	static String usage=" args[0]=all kmers \n args[1]=scored kmers \n args[2]=save";
	
}

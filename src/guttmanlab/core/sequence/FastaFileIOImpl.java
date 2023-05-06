package guttmanlab.core.sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import guttmanlab.core.annotation.Gene;

public class FastaFileIOImpl {
	
	private static Logger logger = Logger.getLogger(FastaFileIOImpl.class.getName());
	
	private BufferedReader reader;
	private String currSeqID;
	private Sequence currentSequence;
	private boolean done;
	
	
	public FastaFileIOImpl() {}
	
	public FastaFileIOImpl(File file) throws IOException {
		//currSeqID = null;
		reader = new BufferedReader(new FileReader(file));
		String line = reader.readLine();
		currSeqID=line.substring(1);
		done=false;
	}
	
	
	
	/**
	 * Read sequences from fasta file and return by name
	 * @param fileName Input fasta
	 * @return Map of sequence name to sequence
	 */
	public static Map<String, Sequence> readFromFileByName(String fileName) {
		Map<String, Sequence> rtrn = new TreeMap<String, Sequence>();
		Collection<Sequence> seqs = readFromFile(fileName);
		for(Sequence seq : seqs) {
			rtrn.put(seq.getName().split(" ")[0], seq);
		}
		return rtrn;
	}
	
	
	public static Map<String, Sequence> readFromFilesByName(File[] files) {
		Map<String, Sequence> rtrn = new TreeMap<String, Sequence>();
		
		for(int i=0; i<files.length; i++) {
			Collection<Sequence> seqs = readFromFile(files[i].getAbsolutePath());
			for(Sequence seq : seqs) {
				rtrn.put(seq.getName().split(" ")[0], seq);
			}

		}
		return rtrn;
	}
	
	public static Map<String, Sequence> readFromFilesByName(File[] files, Collection<Gene> genes) {
		Map<String, Sequence> rtrn = new TreeMap<String, Sequence>();
		
		Collection<String> chromosomes=new TreeSet<String>();
		for(Gene g: genes) {chromosomes.add(g.getReferenceName());}
		
		for(int i=0; i<files.length; i++) {
			String chr=files[i].getName().split("\\.")[0];
			//System.err.println(chr+" "+files[i].getName());
			if(chromosomes.contains(chr)) {
				Collection<Sequence> seqs = readFromFile(files[i].getAbsolutePath());
				for(Sequence seq : seqs) {
					rtrn.put(seq.getName().split(" ")[0], seq);
				}
			}
		}
		return rtrn;
	}
	

	public static Collection<Sequence> readFromFile(String fileName) {
		logger.info("Reading sequences from fasta file " + fileName + "...");
		Collection<Sequence> rtrn = new ArrayList<Sequence>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(fileName));
			boolean started = false;
			String currSeqID = null;
			StringBuilder currSeq = null;
			while(reader.ready()) {
				String line = reader.readLine();
				if(line.startsWith(">")) {
					if(started) {
						rtrn.add(new Sequence(currSeqID, currSeq.toString()));
						logger.info("Added " + currSeqID + " " + currSeq.length());
					}
					currSeqID = line.substring(1);
					currSeq = new StringBuilder();
					continue;
				}
				currSeq.append(line);
				started = true;
			}
			Sequence lastSeq = new Sequence(currSeqID, currSeq.toString());
			rtrn.add(lastSeq);
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		logger.info("Got " + rtrn.size() + " sequences.");
		return rtrn;
	}

	public boolean hasNext() throws IOException {
		if(done) {return false;}
		StringBuilder currSeq = new StringBuilder();
		while(reader.ready()) {
			String line = reader.readLine();
			if(line.startsWith(">")) {
				//Done
				currentSequence=new Sequence(currSeqID, currSeq.toString());
				//logger.info("Added " + currSeqID + " " + currSeq.length());
				currSeqID = line.substring(1);
				return true;
			}
			else{
				currSeq.append(line);
			}	
		}
		currentSequence = new Sequence(currSeqID, currSeq.toString());
		reader.close();
		done=true;
		return true;
	}
	
	public Sequence next() throws IOException {
		return this.currentSequence;
	}
	
	
	

}

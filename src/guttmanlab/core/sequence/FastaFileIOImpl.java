package guttmanlab.core.sequence;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

public class FastaFileIOImpl implements FastaFileIO {
	
	private static Logger logger = Logger.getLogger(FastaFileIOImpl.class.getName());
	
	public FastaFileIOImpl() {}
	
	/**
	 * Read sequences from fasta file and return by name
	 * @param fileName Input fasta
	 * @return Map of sequence name to sequence
	 */
	public Map<String, Sequence> readFromFileByName(String fileName) {
		Map<String, Sequence> rtrn = new TreeMap<String, Sequence>();
		Collection<Sequence> seqs = readFromFile(fileName);
		for(Sequence seq : seqs) {
			rtrn.put(seq.getName(), seq);
		}
		return rtrn;
	}
	
	@Override
	public Collection<Sequence> readFromFile(String fileName) {
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
	

	@Override
	public Iterator<Sequence> iterateThroughFile(String fileName) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void writeToFile(Collection<Sequence> seqs, String fileName, int basesPerLine) {
		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter(fileName));
			for(Sequence seq : seqs) {
				write(seq, bw, basesPerLine);
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	/**
	 * Write one sequence to an output stream
	 * @param seq Sequence
	 * @param bw Buffered writer
	 * @param lineLength Bases per line
	 * @throws IOException
	 */
	public static void write(Sequence seq, BufferedWriter bw, int lineLength) throws IOException {
		
		StringBuilder sequenceBuilder = new StringBuilder(seq.getSequenceBases());
		
		if(seq == null || sequenceBuilder.length() == 0) {
			return;
		}
		if(sequenceBuilder.length() == 0) {
			return;
		}

		int currentIndex = 0;
		bw.write(">" + seq.getName());
		bw.newLine();
		while(currentIndex < sequenceBuilder.length()) {
			int toWrite = Math.min(lineLength, sequenceBuilder.length() - currentIndex) - 1;
			bw.write(sequenceBuilder.substring(currentIndex, currentIndex + toWrite + 1));
			bw.newLine();
			currentIndex = currentIndex + toWrite + 1;
		}
	}

}

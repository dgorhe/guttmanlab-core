package guttmanlab.core.pipeline.util;

import java.io.*;
import java.util.*;

import org.apache.log4j.Logger;


public class FastqParser implements Iterator<FastqSequence>{
	static Logger logger = Logger.getLogger(FastqParser.class.getName());
	Collection<FastqSequence> sequences;
	File fastqFile;
	BufferedReader reader;
	private int numberOfSeq;
	String nextLine = null;
	
	/**
	 * @deprecated This constructure is highly discouraged as it opens a reader. Use the empty constructor instead and set the 
	 * file to read, then call start.
	 * @param fastqFile
	 * @throws IOException 
	 */
	public FastqParser(File fastqFile) throws IOException{
		this.fastqFile=fastqFile;
		reader=new BufferedReader(new InputStreamReader(new FileInputStream(fastqFile)));
		nextLine = reader.readLine();
	}
	
	/**
	 * Empty constructor. Call before setting the file.
	 */
	public FastqParser() {
		super();
	}
	
	/**
	 * Set file and start reader
	 * @param fastqParser The fastq file
	 * @throws IOException
	 */
	public void start(File fastqParser) throws IOException {
		this.fastqFile = fastqParser;
		reader=new BufferedReader(new InputStreamReader(new FileInputStream(fastqFile)));
		nextLine = reader.readLine();
	}
	
	/**
	 * Set reader to the passed reader and start
	 * @param br Reader to set
	 * @throws IOException
	 */
	public void start (BufferedReader br) throws IOException {
		reader=br;
		nextLine = reader.readLine();
	}
	
	/**
	 * Convert to fasta format and save to file
	 * @param save Output file
	 * @throws IOException
	 */
	public void convertToFasta(String save)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
        		String sequence=seq.getSequence();
        		//int num=10;
        		//String lastNBps=getLastBps(sequence, num);
    			//String polyN=polyN("T", num);
    			//if(lastNBps.equalsIgnoreCase(polyN)){System.err.println(sequence);}
    			writer.write(seq.toFasta());
        	}
        }
        writer.close();
	}
	
	/**
	 * Convert to fasta file where read names just differ by a number in the name
	 * @param save Output file
	 * @throws IOException
	 */
	public void convertToNumberedFasta(String save)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		String nextLine;
    	int i=1;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
    			writer.write(seq.toFasta(i));
    			i++;
        	}
        }
        writer.close();
	}
	
	/**
	 * Read an entire fastq file into memory
	 * @param file File
	 * @return Collection of fastq records in the file
	 * @throws IOException
	 */
	public Collection<FastqSequence> parse(File file)throws IOException{
		Collection<FastqSequence> rtrn=new ArrayList<FastqSequence>();
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine,thirdLine, fourthLine);
  			
        		rtrn.add(seq);
        	}
        	
        	
        }
        return rtrn;
	}
	
	/**
	 * Get all records in the file
	 * @return Collection of records in the file
	 * @throws IOException
	 */
	public Collection<FastqSequence> getSequences() throws IOException{
		if(this.sequences!=null){
		return this.sequences;
		}
		else{this.sequences=this.parse(fastqFile); return this.sequences;}
	}

	/**
	 * Divide into multiple fastq files
	 * @param save Output prefix
	 * @param chunkSize Number of records per file
	 * @return The output files
	 * @throws IOException
	 */
	public File[] writeChunks(String save, int chunkSize) throws IOException {
		Collection<File> rtrn=new TreeSet<File>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(this.fastqFile)));
		String nextLine;
		FileWriter writer=new FileWriter(save+".0.fq");
    	int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		if(i%chunkSize ==0){
        			writer.close();
        			rtrn.add(new File(save+"."+(i/chunkSize)+".fq")); 
        			writer=new FileWriter(save+"."+(i/chunkSize)+".fq"); 
        			System.err.println(i);
        		}
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
        		//String sequence=seq.getSequence();
        		writer.write(seq.toFastq()+"\n");
        		i++;
        		if(i % 100000 == 0){System.err.println("Iterating.. "+ i);}
        	}
        	
        	
        }
        this.numberOfSeq=i;
        File[] files=new File[rtrn.size()];
        int counter=0;
        for(File file: rtrn){files[counter++]=file;}
		return files;
	}
	
	/**
	 * @return Number of records in file
	 */
	public int getNumberOfSequences(){
		if(this.numberOfSeq>0){return this.numberOfSeq;}
		else{
			try{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(this.fastqFile)));
			String nextLine;
			int i=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	        	if(nextLine.startsWith("@")){i++;}
	        }
	        this.numberOfSeq=i;
	        reader.close();
			}catch(IOException ex){}
			return this.numberOfSeq;
		}
	}

	public boolean hasNext() {
		return nextLine != null;
	}

	public FastqSequence next() {
		FastqSequence seq = null;
		try{

			
	        if (nextLine  != null) {
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
				String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
	        }
	        nextLine = reader.readLine() ;

		}catch(Exception ex){ 
			logger.error("Exception thrown while reading fastq file",ex);
		}

		return seq;
	}

	/**
	 * Close the reader
	 * @throws IOException
	 */
	public void close() throws IOException{
		reader.close();
	}

	public void remove() {}
	
}

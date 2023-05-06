package guttmanlab.core.proteinSPRITE;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PullBAMFile {
	
	Collection<String> proteins;
	
	public PullBAMFile(File bamFile, File clusters, String save) throws IOException{
		this.proteins=new TreeSet<String>();
		Map<String, Collection<String>> readsToProtein=parseClusters(clusters);
		System.err.println(this.proteins.size());
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
	
		
		Map<String, SAMFileWriter> fileWriter=makeFileWriters(proteins, reader.getFileHeader(), save);
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			if(readsToProtein.containsKey(record.getReadName())){
				Collection<String> proteinNames=readsToProtein.get(record.getReadName());
				for(String protein: proteinNames){
					SAMFileWriter writer=fileWriter.get(protein);
					writer.addAlignment(record);
				}
			}
			
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		
		
		reader.close();
		reads.close();
		close(fileWriter);
	}
	
	private Map<String, Collection<String>> parseClusters(File fileName) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			Collection<String> proteins=new TreeSet<String>();
			Collection<String> reads=new TreeSet<String>();
			String[] tokens=nextLine.split("\t");
			for(int i=1; i<tokens.length; i++){
				if(tokens[i].startsWith("BEAD")){
					proteins.add(tokens[i].split(":")[0]);
				}
				else{reads.add(tokens[i]);}
			}
			this.proteins.addAll(proteins);
			
			
			for(String read: reads){
				rtrn.put(read, proteins);
			}
			counter++;
			if(counter%10000==0){System.err.println(counter);}
		}
		reader.close();
		return rtrn;
	}
	
	
	private Map<String, SAMFileWriter> makeFileWriters(Collection<String> proteins, SAMFileHeader fileHeader, String save) {
		Map<String, SAMFileWriter> rtrn=new TreeMap<String, SAMFileWriter>();
		
		for(String protein: proteins){
			SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(save+"."+protein+".bam"));
			rtrn.put(protein, writer1);
		}
		
		return rtrn;
	}

	

	private void close(Map<String, SAMFileWriter> fileWriter) {
		for(String name: fileWriter.keySet()){fileWriter.get(name).close();}
		
	}

	public PullBAMFile(File bamFile, File clusters, String protein, double percent, String save) throws IOException{
		Collection<String> readNames= parseClusters(clusters, percent, protein);
		
		System.err.println(readNames.size());
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(reader.getFileHeader(), false, new File(save+".bam"));
		
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(overlaps(record, readNames)){
				writer1.addAlignment(record);
				//add(record, map);
			}
			counter++;
			if(counter%100000==0){System.err.println(counter + " "+record.getReadName()+" "+readNames.size());}
		}
		
		
		reader.close();
		reads.close();
		writer1.close();
		
		//write(save, map);
	}

	private Collection<String> parseClusters(File fileName, double percent, String protein) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			List<String> proteins=new ArrayList<String>();
			Collection<String> reads=new TreeSet<String>();
			String[] tokens=nextLine.split("\t");
			for(int i=1; i<tokens.length; i++){
				if(tokens[i].startsWith("BEAD")){
					proteins.add(tokens[i].split(":")[0].split("_")[1]);
				}
				else{reads.add(tokens[i]);}
			}
			
			if(percent(proteins, protein)>percent){rtrn.addAll(reads);}
		}
		reader.close();
		return rtrn;
	}

	private double percent(List<String> proteins, String protein) {
		int count=0;
		int total=proteins.size();
		
		for(String p: proteins){
			if(p.equalsIgnoreCase(protein)){count++;}
		}
		
		return (double)count/(double)total;
	}

	

	private boolean overlaps(SAMRecord record, Collection<String> readNames) {
		if(readNames.contains(record.getReadName())){return true;}
		return false;
	}

	
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File bamFile=new File(args[0]);
			File clusters=new File(args[1]);
			String save=args[2];
			new PullBAMFile(bamFile, clusters, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam \n args[1]=clusters \n args[2]=save (base name)";
}

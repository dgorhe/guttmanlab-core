package guttmanlab.core.scrna;

import java.io.File;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SplitBAMByCell {

	public SplitBAMByCell(File bam, String save) {
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		
		Map<String, SAMFileWriter> writers=new TreeMap<String, SAMFileWriter>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			//System.err.println(record.getReadName());
			String barcode=record.getReadName().split("::")[1];
			barcode=toFileName(barcode);
			//System.err.println(record.getReadName()+" "+barcode);
			//String umi=record.getReadName().split(":")[0].split("_")[1];
			add(writers, barcode, record, save, reader.getFileHeader());
			
			counter++;
			if(counter%100000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		close(writers);
	}

	private String toFileName(String barcode) {
		String rtrn="";
		barcode=barcode.replaceAll("\\[", "");
		String[] tokens=barcode.split("\\]");
		for(int i=0; i<tokens.length; i++) {
			if(i!=0) {rtrn+="_";}
			rtrn+=tokens[i];
			
		}
		return rtrn;
	}

	private void close(Map<String, SAMFileWriter> writers) {
		for(String b: writers.keySet()) {
			writers.get(b).close();
		}
		
	}

	private void add(Map<String, SAMFileWriter> writers, String barcode, SAMRecord record, String fileName, SAMFileHeader fileHeader) {
		
		if(!writers.containsKey(barcode)) {
			SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(fileName+"."+barcode+".bam"));
			writers.put(barcode, writer);
		}
		
		SAMFileWriter writer=writers.get(barcode);
		writer.addAlignment(record);
	}
	
	
	public static void main(String[] args) {
		File bam=new File(args[0]);
		String save=args[1];
		new SplitBAMByCell(bam, save);
	}
	
	
}

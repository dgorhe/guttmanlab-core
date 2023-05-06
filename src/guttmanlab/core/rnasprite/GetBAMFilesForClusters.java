package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class GetBAMFilesForClusters {

	
	public GetBAMFilesForClusters(BarcodingDataStreaming data, String gene, File bamFile, String save){
		
		Collection<String> barcodes=getBarcodes(data, gene);
		
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMFileHeader header=reader.getFileHeader();
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(header, false, new File(save));
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			String barcode=getBarcode(record);
			if(barcodes.contains(barcode)){
				//System.err.println(barcode);
				writer.addAlignment(record);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter+" "+barcode);}
		}
		writer.close();
		reads.close();
		reader.close();
		
	}

	private String getBarcode(SAMRecord record) {
		String name=record.getReadName();
		String barcode=name.split("::")[1];
		String newBarcode=splitOnBrackets(barcode);
		
		return newBarcode+name.split(":")[2];
	}

	private String splitOnBrackets(String barcode) {
		String rtrn="";
		String[] tokens=barcode.split("\\[");
		for(int i=1; i<tokens.length-1; i++){
			String name=tokens[i].split("\\]")[0];
			rtrn=rtrn+name+".";
		}
		rtrn.replaceAll("RPM.", "");
		rtrn.replaceAll("DPM.", "");
		return rtrn;
	}

	private Collection<String> getBarcodes(BarcodingDataStreaming data, String gene) {
		Collection<String> rtrn=new TreeSet<String>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.containsRNA(gene)){
				rtrn.add(c.getBarcode());
			}
		}
		
		data.close();
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String gene=args[1];
		File bam=new File(args[2]);
		String save=args[3];
		new GetBAMFilesForClusters(data, gene, bam, save);
	}
	
}

package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class FastqToClusters {

	/*
	    @FS10000829:63:BPC29511-2024:1:1101:1010:1000::[DAP_K27me3][NYBot4_Stg][Even2Bo88][D10][Even2Bo18][F3]
		AAGGTAGCTAGTTCGTACCTTGTCTAATGGCACCAGAGTATGTGTTGTCA
		+
		FFFFFFFF:FFFFFFFFF:F,F:F:FFFFFFFFFFFFFFFFF,FF:F:F,
		@FS10000829:63:BPC29511-2024:1:1101:1040:1000::[DAP_K36me3][NYBot6_Stg][Even2Bo88][C5][Even2Bo26][C9]
		ATTCGGTGGCCGCTTCCCATGTCATGCGATCCAGATTGAGAAGTTGTCAC
		+
		F:FFFF,FFFFFFFFFFFFF:FFF,FFFFFFFFFFFFFFFFFFFFFFFF:
		@FS10000829:63:BPC29511-2024:1:1101:1050:1000::[DAP_K27me3][NYBot18_Stg][Even2Bo64][D8][Even2Bo22][F6]
		AAGGTAGCTACCGTCCTACTTGTCTACAATTATATTGCGGTGTGTTGTCA

	I typically split the header at the "::" delimiter. I only create clusters based off of barcodes 2-6. Since barcode1 is the oligo ID.
	For example, this oligo @FS10000829:63:BPC29511-2024:1:1101:1010:1000::[DAP_K27me3][NYBot4_Stg][Even2Bo88][D10][Even2Bo18][F3] belongs to this cluster [NYBot4_Stg][Even2Bo88][D10][Even2Bo18][F3].
	In place of chromosome information I add the oligo name, which is Barcode#1 to the cluster file and the UMI, which is in position 10-16 of the read.
	Here is an example cluster:
	NYBot9_Stg.Even2Bo64.D6.Even2Bo30.F2	DAP_CTCF:GACTCAA

	I also do not populate clusters if they have a NOTFOUND in them since those are incomplete barcode calls.
	Here is file location:
	/groups/guttman/primarydata/sequencingruns/20200927/A1/A1_S1_L001_R1_001.fastq.gz.bcid.fastq.gz
	 */
	
	public FastqToClusters(File fastq, String save) throws IOException{
		Map<String, Collection<String>> clusters=new TreeMap<String, Collection<String>>();
		FastqReader p1=new FastqReader(fastq);
		int counter=0;
		while(p1.hasNext()){
			FastqRecord record=p1.next();
			String barcode=getBarcode(record);
			String umi=getUMI(record);
			String protein=getProteinID(record);
			Collection<String> temp=new TreeSet<String>();
			if(clusters.containsKey(barcode)){temp=clusters.get(barcode);}
			temp.add(protein+":"+umi);
			clusters.put(barcode, temp);
			counter++;
			if(counter%100000==0){System.err.println(counter+" "+record.getReadString()+"\t"+umi+"\t"+record.getReadHeader()+"\t"+barcode+" "+protein);}
		}
		p1.close();
		write(save, clusters);
	}
	
	private String getBarcode(FastqRecord record) {
		//	I typically split the header at the "::" delimiter. I only create clusters based off of barcodes 2-6. Since barcode1 is the oligo ID.
		//For example, this oligo @FS10000829:63:BPC29511-2024:1:1101:1010:1000::[DAP_K27me3][NYBot4_Stg][Even2Bo88][D10][Even2Bo18][F3] belongs to this cluster [NYBot4_Stg][Even2Bo88][D10][Even2Bo18][F3].
		//In place of chromosome information I add the oligo name, which is Barcode#1 to the cluster file and the UMI, which is in position 10-16 of the read.
		//Here is an example cluster:
		String rtrn="";
		String barcodes=record.getReadHeader().split("::")[1];
		String[] tokens=barcodes.split("]");
		for(int i=1; i<tokens.length; i++){
			rtrn+=tokens[i].replace("[", "");
			if(i!=tokens.length-1){rtrn+=".";}
		}
		return rtrn;
	}
	
	private String getProteinID(FastqRecord record) {
		//	I typically split the header at the "::" delimiter. I only create clusters based off of barcodes 2-6. Since barcode1 is the oligo ID.
		//For example, this oligo @FS10000829:63:BPC29511-2024:1:1101:1010:1000::[DAP_K27me3][NYBot4_Stg][Even2Bo88][D10][Even2Bo18][F3] belongs to this cluster [NYBot4_Stg][Even2Bo88][D10][Even2Bo18][F3].
		//In place of chromosome information I add the oligo name, which is Barcode#1 to the cluster file and the UMI, which is in position 10-16 of the read.
		//Here is an example cluster:
		
		String barcodes=record.getReadHeader().split("::")[1];
		String[] tokens=barcodes.split("]");
		return tokens[0].replace("[", "");
	}

	private String getUMI(FastqRecord record) {
		//UMI is in position 10-16 of the read
		
		return record.getReadString().substring(9, 16);
	}

	private void write(String save, Map<String, Collection<String>> clusters) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: clusters.keySet()){
			Collection<String> names=clusters.get(barcode);
			writer.write(barcode+"\t"+names.size());
			for(String name: names){writer.write("\t"+name);}
			writer.write("\n");
		}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		File fastq=new File(args[0]);
		String save=args[1];
		new FastqToClusters(fastq, save);
	}
	
	
}

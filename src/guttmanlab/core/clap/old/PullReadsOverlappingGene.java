package guttmanlab.core.clap.old;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.test.BedToBedgraph;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PullReadsOverlappingGene {

	public PullReadsOverlappingGene(File bam, Gene gene, String save) throws IOException {
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator iter=inputReader.queryOverlapping(gene.getReferenceName(), gene.getReferenceStartPosition(), gene.getReferenceEndPosition());
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(inputReader.getFileHeader(), false, new File(save));
		
		
		Collection<SingleInterval> regions=new TreeSet<SingleInterval>();
		
		int counter=0;
		while(iter.hasNext()) {
			SAMRecord record=iter.next();
			if(!record.getNotPrimaryAlignmentFlag() && record.getMappingQuality()>1){
				SAMFragment f=new SAMFragment(record);
				if(record.getFirstOfPairFlag()) {
					SingleInterval interval=getRegion(record, f.getOrientation());
					if(interval.getGenomicLength()<1000) {
						regions.add(interval);
						counter++;
						//System.out.println(interval.toBED());
					}
				}
				
				
				if(f.getOrientation().equals(gene.getOrientation())) {
					writer.addAlignment(record);
				}
			}
			counter++;
		}
		
		
		double expected=(double)counter/(double)gene.size();
		System.err.println(counter+" "+gene.size()+" "+expected);		
		iter.close();
		inputReader.close();
		writer.close();
		
		
		
		new BedToBedgraph(regions, save+".bedgraph", expected);
		
	}
	
	private SingleInterval getRegion(SAMRecord record, Strand o) {
		return new SingleInterval(record.getReferenceName(), Math.min(record.getAlignmentStart(), record.getMateAlignmentStart()), Math.max(record.getAlignmentStart(), record.getMateAlignmentStart())+35, o, record.getReadName());
		
		
	}

	public static void main(String[] args) throws IOException {
		File bam=new File(args[0]);
		Gene gene=BEDFileIO.loadRegionsFromFile(args[1]).iterator().next();
		String save=args[2];
		new PullReadsOverlappingGene(bam, gene, save);
	}
	
}

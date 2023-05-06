package guttmanlab.core.snps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class QuantifyGeneExpression {

	public QuantifyGeneExpression(Collection<File> files, File genes, String save) throws IOException{
		//Map<String, Pair<File>> snpFiles=getSNPFiles(files);
		
		MatrixWithHeaders mwh=quantify(files, genes);
		
		mwh.write(save);
		
	}

	

	private Map<String, Pair<File>> getSNPFiles(File[] files) {
		Map<String, Pair<File>> rtrn=new TreeMap<String, Pair<File>>();
		
		for(int i=0; i<files.length; i++){
			File file=files[i];
			if(file.getName().endsWith(".bam")){
				String name=file.getName().split("\\.")[0];
				String snp=file.getName().split("\\.")[2];
				Pair<File> pair=new Pair<File>();
				if(rtrn.containsKey(name)){
					pair=rtrn.get(name);
				}
				pair.setName(name);
				System.err.println(file.getName()+" "+name+" "+snp);
				if(snp.equalsIgnoreCase("SNP1")){pair.setValue1(file);}
				else if(snp.equalsIgnoreCase("SNP2")){pair.setValue2(file);}
				rtrn.put(name, pair);
			}
			
		}
		
		return rtrn;
	}



	private MatrixWithHeaders quantify(Collection<File> bamFiles, File geneFile) throws IOException {
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(geneFile.getAbsolutePath());
		List<String> rowNames=getNames(genes);
		Map<String, String> position=getPositions(genes);
		List<String> columnNames=getColumnNames(bamFiles);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rowNames, columnNames);
		mwh.setPIDToName(position);
		
		Map<String, IntervalTree<Gene>> tree=BEDFileIO.loadTree(geneFile.getAbsolutePath());
		
		for(File bam: bamFiles){
			System.err.println(bam.getName());
			quantify(bam, tree, mwh);
		}
		return mwh;
	}

	

	private void quantify(File bam, Map<String, IntervalTree<Gene>> genes, MatrixWithHeaders mwh) {
		String name1=bam.getName();
		double total=quant(name1, genes, bam, mwh);
		
		for(String row: mwh.getRowNames()){
			double val=1000000.0*(mwh.get(row, name1))/total;
			mwh.set(row, name1, val);
		}
		
	}

	private double quant(String name, Map<String, IntervalTree<Gene>> genes, File bamFile, MatrixWithHeaders mwh) {
		System.err.println(bamFile.getAbsolutePath());
		SAMFileReader reader=new SAMFileReader(bamFile);
		double total=0;
		SAMRecordIterator reads=reader.iterator();
		
		
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			Collection<String> list=getGenes(record, genes);
			for(String g: list){
				mwh.incrementCount(g, name);
			}
			total++;
		}
		
		reads.close();
		reader.close();
		return total;
	}

	private Collection<String> getGenes(SAMRecord record, Map<String, IntervalTree<Gene>> genes) {
		Collection<String> rtrn=new TreeSet<String>();
		
		if(genes.containsKey(record.getReferenceName())){
			IntervalTree<Gene> tree= genes.get(record.getReferenceName());
			Iterator<Gene> iter=tree.overlappingValueIterator(record.getAlignmentStart()-1, record.getAlignmentEnd()+1);
			while(iter.hasNext()){
				SAMFragment frag=new SAMFragment(record);
				Gene g=iter.next();
				if(frag.getOrientation().equals(g.getOrientation())){
					rtrn.add(g.getName());
				}
			
			}
		}
		
		return rtrn;
	}

	private List<String> getColumnNames(Collection<File> snpFiles) {
		List<String> rtrn=new ArrayList<String>();
		
		
		for(File file: snpFiles){
			String name=file.getName();
			rtrn.add(name);
		}
		
		return rtrn;
	}

	private Map<String, String> getPositions(Collection<Gene> genes) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(Gene gene: genes){
			rtrn.put(gene.getName(), gene.toUCSC());
		}
		return rtrn;
	}

	private List<String> getNames(Collection<Gene> genes) {
		List<String> rtrn=new ArrayList<String>();
		
		for(Gene gene: genes){
			if(!rtrn.contains(gene.getName())){
				rtrn.add(gene.getName());
			}
		}
		return rtrn;
	}
	
	
	private static File getBam(File dir) {
		File[] files=dir.listFiles();
		
		for(int i=0; i<files.length; i++){
			if(files[i].getName().endsWith(".bam")){return files[i];}
		}
		return null;
	}
	
	
	private static Collection<File> getFiles(String dir) {
		List<File> files=new ArrayList<File>();
		//String dir="/groups/guttman/jamie/data/sharp/heard_data_nature_2020/correctFq/alignments/mm10_default/star/";
		File[] dirs=new File(dir).listFiles();
		for(int i=0; i<dirs.length; i++){
			String name=dirs[i].getName();
			File bam=getBam(dirs[i]);
			if(bam!=null){files.add(bam);}
		}
		return files;
	}
	
	public static void main(String[] args) throws IOException{
		Collection<File> files=getFiles(args[0]);
		File bedFile=new File(args[1]);
		String save=args[2];
		new QuantifyGeneExpression(files, bedFile, save);
	}



	
	
}

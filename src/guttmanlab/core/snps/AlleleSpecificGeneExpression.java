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

public class AlleleSpecificGeneExpression {

	public AlleleSpecificGeneExpression(File[] files, File genes, String save) throws IOException{
		Map<String, Pair<File>> snpFiles=getSNPFiles(files);
		
		MatrixWithHeaders mwh=quantify(snpFiles, genes);
		
		mwh.write(save);
		
	}

	

	private Map<String, Pair<File>> getSNPFiles(File[] files) {
		Map<String, Pair<File>> rtrn=new TreeMap<String, Pair<File>>();
		
		for(int i=0; i<files.length; i++){
			File file=files[i];
			if(file.getName().endsWith(".bam")){
				//System.err.println(file.getName());
				String name=file.getName().split("\\.")[0];
				String snp=file.getName().split("\\.")[1];
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



	private MatrixWithHeaders quantify(Map<String, Pair<File>> snpFiles, File geneFile) throws IOException {
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(geneFile.getAbsolutePath());
		List<String> rowNames=getNames(genes);
		Map<String, String> position=getPositions(genes);
		List<String> columnNames=getColumnNames(snpFiles);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rowNames, columnNames);
		mwh.setPIDToName(position);
		
		Map<String, IntervalTree<Gene>> tree=BEDFileIO.loadTree(geneFile.getAbsolutePath());
		
		for(String name: snpFiles.keySet()){
			Pair<File> file=snpFiles.get(name);
			quantify(file, tree, mwh);
		}
		return mwh;
	}

	

	private void quantify(Pair<File> file, Map<String, IntervalTree<Gene>> genes, MatrixWithHeaders mwh) {
		String name1=file.getName()+"_snp1";
		quant(name1, genes, file.getValue1(), mwh);
		
		String name2=file.getName()+"_snp2";
		quant(name2, genes, file.getValue2(), mwh);
		
	}

	private void quant(String name, Map<String, IntervalTree<Gene>> genes, File bamFile, MatrixWithHeaders mwh) {
		System.err.println(bamFile.getAbsolutePath());
		SAMFileReader reader=new SAMFileReader(bamFile);
		
		SAMRecordIterator reads=reader.iterator();
		
		
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			Collection<String> list=getGenes(record, genes);
			for(String g: list){
				mwh.incrementCount(g, name);
			}
		}
		
		reads.close();
		reader.close();
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

	private List<String> getColumnNames(Map<String, Pair<File>> snpFiles) {
		List<String> rtrn=new ArrayList<String>();
		
		
		for(String s: snpFiles.keySet()){
			Pair<File> file=snpFiles.get(s);
			String name=file.getName();
			rtrn.add(name+"_snp1");
			rtrn.add(name+"_snp2");
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
	
	
	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		File bedFile=new File(args[1]);
		String save=args[2];
		new AlleleSpecificGeneExpression(files, bedFile, save);
	}
	
}

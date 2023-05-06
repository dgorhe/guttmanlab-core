package guttmanlab.core.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class GeneExpressionUtils {

	
	public static MatrixWithHeaders quantify(File[] bams, String genes) throws IOException{
		MatrixWithHeaders rtrn=null;
		Map<String, IntervalTree<Gene>> trees= BEDFileIO.loadTree(genes);
		List<String> geneNames=BEDFileIO.getGeneNames(genes);
		
		for(int i=0; i<bams.length; i++){
			MatrixWithHeaders column=quantify(bams[i], trees, geneNames);
			rtrn=merge(rtrn, column);
		}
		
		return rtrn;
	}

	private static MatrixWithHeaders quantify(File file, Map<String, IntervalTree<Gene>> trees, List<String> geneNames) {
		List<String> column=new ArrayList<String>();
		column.add(file.getName());
		MatrixWithHeaders rtrn=new MatrixWithHeaders(geneNames, column);
		
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			if(trees.containsKey(record.getReferenceName())){
				IntervalTree<Gene> tree=trees.get(record.getReferenceName());
				Iterator<Gene> genes=tree.overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
				Collection<String> names=new TreeSet<String>();
				while(genes.hasNext()){
					Gene gene=genes.next();
					names.add(gene.getName());
				}
				for(String name: names){
					rtrn.incrementCount(name, file.getName());
				}
			}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		return rtrn;
	}

	

	private static MatrixWithHeaders merge(MatrixWithHeaders rtrn, MatrixWithHeaders column) {
		if(rtrn==null){return column;}
		rtrn.appendColumns(column);
		return rtrn;
	}
	
}

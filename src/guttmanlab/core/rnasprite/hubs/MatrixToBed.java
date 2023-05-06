package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class MatrixToBed {

	public static void main(String[] args) throws IOException {
		MatrixWithHeaders matrix=new MatrixWithHeaders(new File(args[0]));
		String saveDir=args[1];
		Map<String, IntervalTree<String>> geneNames=BEDFileIO.loadGeneNamesFromRefFlat(args[2]);
		
		MatrixWithHeaders geneMatrix=convertToGeneMatrix(matrix, geneNames);
		geneMatrix.write(saveDir+"/lncRNAxgene.matrix");
		
		writeGeneToLncRNA(saveDir+"/geneTolncRNA.list", geneMatrix);
		writeLncRNAToGene(saveDir+"/lncRNAToGene.list", geneMatrix);
		
		//FileWriter writerFull=new FileWriter(saveDir+"/merged.genes");
		
		for(String rna: matrix.getRowNames()) {
			String save=saveDir+"/"+rna+".bed";
			String saveBedGraph=saveDir+"/"+rna+".bedgraph";
			FileWriter writer=new FileWriter(save);
			
			FileWriter writerBedgraph=new FileWriter(saveBedGraph);
			
			List<SingleInterval> regions=new ArrayList<SingleInterval>();
			for(String pos: matrix.getColumnNames()) {
				double score=matrix.get(rna, pos);
				SingleInterval r=new SingleInterval(pos);
				writerBedgraph.write(r.toBedgraph(score)+"\n");
				if(score>5) {
					writer.write(r.toShortBED()+"\n");
					//System.err.println(rna+" "+pos+" "+score);
					regions.add(r);
				}
				
			}
			writer.close();
			writerBedgraph.close();
			
			/*Collection<String> genes=getGenes(geneNames, regions);
			
			writerFull.write(rna+"\t"+genes.size());
			for(String gene: genes) {writerFull.write("\t"+gene);}
			writerFull.write("\n");*/
			
			
		}
		
		//writerFull.close();
	}

	private static void writeLncRNAToGene(String save, MatrixWithHeaders geneMatrix) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String lncRNA: geneMatrix.getRowNames()) {
			Collection<Score> list=new TreeSet<Score>();
			for(String gene: geneMatrix.getColumnNames()) {
				double score=geneMatrix.get(lncRNA, gene);
				if(score>5) {
					Score newScore=new Score(gene, score);
					list.add(newScore);
				}
			}
			
			writer.write(lncRNA);
			for(Score score: list) {
				writer.write("\t"+score.name);
			}
			
			/*for(Score score: list) {
				writer.write("\t"+score.score);
			}*/
			
			writer.write("\n");
		}
		
		
		writer.close();
		
	}

	private static void writeGeneToLncRNA(String save, MatrixWithHeaders geneMatrix) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String gene: geneMatrix.getColumnNames()) {
			Collection<Score> list=new TreeSet<Score>();
			for(String lncRNA: geneMatrix.getRowNames()) {
				double score=geneMatrix.get(lncRNA, gene);
				if(score>5) {
					Score newScore=new Score(lncRNA, score);
					list.add(newScore);
					//writer.write("\t"+lncRNA+"_"+score);
				}
			}
			
			writer.write(gene);
			for(Score score: list) {
				writer.write("\t"+score.name);
			}
			
			for(Score score: list) {
				writer.write("\t"+score.score);
			}
			
			writer.write("\n");
		}
		
		
		writer.close();
	}

	private static MatrixWithHeaders convertToGeneMatrix(MatrixWithHeaders matrix, Map<String, IntervalTree<String>> geneNames) {
		Map<String, Collection<String>> geneToBin=getGeneToBin(matrix.getColumnNames(), geneNames);
		
		List<String> genes=getList(geneNames);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(matrix.getRowNames(), genes);
		
		for(String row: rtrn.getRowNames()) {
			for(String column: rtrn.getColumnNames()) {
				Collection<String> posBins=new TreeSet<String>();
				if(geneToBin.containsKey(column)) {
					posBins=geneToBin.get(column);
				}
				double maxScore=0;
				for(String pos: posBins) {
					if(matrix.containsColumn(pos)) {
						maxScore=Math.max(maxScore, matrix.get(row, pos));
					}
				}
				rtrn.set(row, column, maxScore);
			}
		}
		
		int count1=0;
		int count2=0;
		List<String> sub=new ArrayList<String>();
		
		for(String gene: rtrn.getColumnNames()) {
			double max=Statistics.max(rtrn.getColumn(gene));
			if(max<5) {count1++;}
			else {count2++; sub.add(gene);}
		}
		
		rtrn=rtrn.submatrixByColumnNames(sub);
		System.err.println(count1+" "+count2);
		return rtrn;
	}

	private static Map<String, Collection<String>> getGeneToBin(List<String> positions, Map<String, IntervalTree<String>> geneNames) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		int resolution=new SingleInterval(positions.iterator().next()).getGenomicLength();
		
		for(String chr: geneNames.keySet()) {
			IntervalTree<String> tree=geneNames.get(chr);
			Iterator<Node<String>> iter=tree.iterator();
			while(iter.hasNext()) {
				Node<String> val=iter.next();
				SingleInterval r=new SingleInterval(chr, val.getStart(), val.getEnd());
				Collection<SingleInterval> regions=r.allBins(resolution);
				Collection<String> list=convert(regions);
				rtrn.put(val.getValue(), list);
			}
		}
		
		return rtrn;
	}

	private static Collection<String> convert(Collection<SingleInterval> regions) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(SingleInterval r: regions) {
			rtrn.add(r.toUCSC());
		}
		
		return rtrn;
	}

	private static List<String> getList(Map<String, IntervalTree<String>> geneNames) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String chr: geneNames.keySet()) {
			IntervalTree<String> tree=geneNames.get(chr);
			Iterator<String> iter=tree.valueIterator();
			while(iter.hasNext()) {
				String name=iter.next();
				if(!rtrn.contains(name)) {
					rtrn.add(name);
				}
			}
		}
		
		return rtrn;
	}

	private static Collection<String> getGenes(Map<String, IntervalTree<String>> geneNames, List<SingleInterval> regions) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(SingleInterval region: regions) {
			if(geneNames.containsKey(region.getReferenceName())) {
				IntervalTree<String> tree=geneNames.get(region.getReferenceName());
				Iterator<String> names=tree.overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
				while(names.hasNext()) {
					rtrn.add(names.next());
				}
			}
		}
		
		return rtrn;
	}
	
	
	
	
}

package guttmanlab.core.util;

import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;

public class RefFlatTest {

	public static void main(String[] args) throws IOException {
		Map<String, IntervalTree<String>> trees=BEDFileIO.loadRefFlatTree(args[0]);
		
		for(String chr: trees.keySet()) {
			IntervalTree<String> tree=trees.get(chr);
			Iterator<Node<String>> iter=tree.iterator();
			while(iter.hasNext()) {
				Node<String> node=iter.next();
				String name=node.getValue();
				System.out.println(chr+"\t"+node.getStart()+"\t"+node.getEnd()+"\t"+name);
			}
		}
		
	}
	
}

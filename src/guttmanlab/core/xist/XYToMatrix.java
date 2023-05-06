package guttmanlab.core.xist;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class XYToMatrix {

	public XYToMatrix(File file, String save) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file.getAbsolutePath(),1);
		
		List<String> rows=getRows(lines);
		List<String> columns=getColumns(lines);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			if(tokens.length>2) {
				mwh.set(tokens[0], tokens[1], Double.parseDouble(tokens[2]));
			}
		}
		
		mwh.write(save);
	}
	
	private List<String> getColumns(List<String> lines) {
		Set<String> rtrn=new TreeSet<String>();
			
		for(String line: lines) {
			if(line.split("\t").length>2) {
				rtrn.add(line.split("\t")[1]);
			}
		}
			
		ArrayList<String> list=new ArrayList<String>();
		list.addAll(rtrn);
			
		return list;
	}

	private List<String> getRows(List<String> lines) {
		Set<String> rtrn=new TreeSet<String>();
		
		for(String line: lines) {
			if(line.split("\t").length>2) {
				rtrn.add(line.split("\t")[0]);
			}
		}
		
		ArrayList<String> list=new ArrayList<String>();
		list.addAll(rtrn);
		
		return list;
	}

	public static void main(String[] args) throws IOException {
		new XYToMatrix(new File(args[0]), args[1]);
		
		/*File[] files=new File(args[0]).listFiles();
		
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i].getAbsolutePath());
			new XYToMatrix(files[i], files[i].getAbsolutePath()+".matrix");
		}*/
	}
	
}

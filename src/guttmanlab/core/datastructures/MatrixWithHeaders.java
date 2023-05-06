package guttmanlab.core.datastructures;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import Jama.Matrix;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class MatrixWithHeaders {
	private Matrix data;
	private LinkedHashMap<String, Integer> columnIndexMap;
	private LinkedHashMap<String, Integer> rowIndexMap;
	private LinkedHashMap<String, List<Integer>> rowDescrIndexMap;
	private List<String> columns;
	private List<String> rows;
	private Map<String, String> pidToName;
	private TreeMap<String, String> probeClass;
	private TreeMap<String, Collection<Integer>> replicateMap;
	
	protected MatrixWithHeaders() {
		super();
	}
	
	
	public MatrixWithHeaders(List<String> rows, List<String> columns) {
		this.data = new Matrix(rows.size(), columns.size());
		initNameIndexMaps(rows, columns);
	}
	
	
	public MatrixWithHeaders(File dataFile, Map<String, String> geneToChr) throws IOException {
		initFromTxt(dataFile, geneToChr);
	}
	
	public MatrixWithHeaders(File dataFile) throws IOException {
		initFromTxt(dataFile);
	}
	
	public MatrixWithHeaders(BufferedReader reader, String sep, int linesToSkip) throws IOException {
		initFromTxt(reader, sep, linesToSkip);
	}
	
	public MatrixWithHeaders(BufferedReader reader, String sep, int linesToSkip, int timeShift) throws IOException {
		initFromTxt(reader, sep, linesToSkip, timeShift);
	}

	public MatrixWithHeaders(File dataFile, String sep) throws IOException {
		initFromTxt(dataFile, sep);
	}

	public Map<String, String> getPIDToName(){return this.pidToName;}

	public void write(BufferedWriter bw) throws IOException {
		//write("row");
		bw.write("");
		//if(pidToName!=null && !pidToName.isEmpty()){bw.write("\tCoordinate");}
		for(String columnName : columnIndexMap.keySet()) {
			//bw.write("Name\t");
			
			bw.write("\t"+columnName);
		}
		bw.newLine();
		bw.flush();
		for(String rowName : rowIndexMap.keySet()) {
			bw.write(rowName);
			//if(pidToName!=null && !pidToName.isEmpty()){bw.write("\t"+pidToName.get(rowName));}
			for(String colName : columnIndexMap.keySet()) {
				bw.write("\t");
				bw.write(String.valueOf(get(rowName, colName)));
			}
			bw.newLine();
		}
	}
	
	public void write(String fileName) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		write(bw);
		bw.close();
		
	}
	
	public void writeGCT(BufferedWriter bw) throws IOException {
		bw.write("#1.2");
		bw.newLine();
		bw.write(data.getRowDimension() +"\t"+data.getColumnDimension());
		bw.newLine();
		
		bw.write("name\tdescription");
		for(String columnName : columnIndexMap.keySet()) {
			bw.write("\t");
			bw.write(columnName);
		}
		bw.newLine();
		for(String rowName : rowIndexMap.keySet()) {
			bw.write(rowName);
			bw.write("\t");
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){bw.write(rowName);}
			else{bw.write(pidToName.get(rowName));}
			for(String colName : columnIndexMap.keySet()) {
				bw.write("\t");
				bw.write(String.valueOf(get(rowName, colName)));
			}
			bw.newLine();
		}
	}
	
	public void writeGCT(String save) throws IOException {
		writeGCT(save, this.getColumnNames(), this.getRowNames());
	}
	
	public void writeGCT(String save, Collection<String> columns) throws IOException {
		writeGCT(save, columns, this.getRowNames());
	}
	
	
	public void writeGCT(String save, Map<String, Collection<String>> experimentInfo, boolean ordered) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		
		int columnSize=0;
		for(String group: experimentInfo.keySet()){
			Collection<String> experiment=experimentInfo.get(group);
			for(String e: experiment){
				columnSize++;
			}
		}
		
		writer.write(rows.size() +"\t"+columnSize+"\n");
				
		writer.write("name\tdescription");
		
		
		for(String group: experimentInfo.keySet()){
			Collection<String> experiment=experimentInfo.get(group);
			for(String e: experiment){
				writer.write("\t"+e);
			}
		}
		writer.write("\n");
		
		
		for(String rowName : this.getRowNames()) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			
			for(String group: experimentInfo.keySet()){
				Collection<String> experiment=experimentInfo.get(group);
				for(String e: experiment){
					writer.write("\t"+get(rowName, e));
				}
			}
			
			writer.write("\n");
		}
		writer.close();
	}
	
	
	public void writeGCT(String save, Collection<String> columns, Collection<String> rows) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		for(String columnName : columns) {
			writer.write("\t"+columnName);
		}
		writer.write("\n");
		for(String rowName : rows) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			
			for(String colName : columns) {
				writer.write("\t"+get(rowName, colName));
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	public void writeGCTWithHeaders(String save, Map<String, String> columnMapping) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		for(String columnName : columns) {
			String info=columnMapping.get(columnName);
			writer.write("\t"+info);
		}
		writer.write("\n");
		for(String rowName : rows) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			
			for(String colName : columns) {
				writer.write("\t"+get(rowName, colName));
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	//SPECIFY THE PRECISION FOR PRINTING DOUBLE
	public void writeGCT(String save, Collection<String> columns, Collection<String> rows, int precision) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		for(String columnName : columns) {
			writer.write("\t"+columnName);
		}
		writer.write("\n");
		for(String rowName : rows) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			
			for(String colName : columns) {
				writer.write("\t");
				String s= new Double(get(rowName, colName)).toString();
				if (s.length()<2+precision) writer.write (s); 
				else writer.write (s,0,2+precision);
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	
	public void writeGCT(String save, Map<String, String> columns, Map<String, String> rows) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		for(String columnName : columns.keySet()) {
			writer.write("\t"+columns.get(columnName));
		}
		writer.write("\n");
		for(String rowName : rows.keySet()) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			for(String colName : columns.keySet()) {
				writer.write("\t"+get(rowName, colName));
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	public void writeGCT(String save, Map<String, String> columns, Collection<String> rows) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		int counter=0;
		for(String columnName : columns.keySet()) {
			writer.write("\t"+columns.get(columnName)+"_"+counter);
			counter++;
		}
		writer.write("\n");
		for(String rowName : rows) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			for(String colName : columns.keySet()) {
				writer.write("\t"+get(rowName, colName));
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	public void writeGCT(String save, Map<String, String> columns) throws IOException {
		writeGCT(save, columns, this.getRowNames());
	}
	
	
	
	/**
	 * Takes the inverse (or pseudo inverse if the underlying data is not a squared matrix)
	 * @return
	 */
	public Matrix dataInverse() {
		return data.inverse();		
	}
	
	public int columnDimension() { return data.getColumnDimension();}
	public int rowDimension() { return data.getRowDimension();}
	public double get(int i, int j) {return data.get(i,j);}
	
	public double[] getColumn(int j){
		double[] rtrn=new double[rowDimension()];
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=get(i,j);
		}
		return rtrn;
	}
	
	public double[] getColumn(String columnName){
		double[] rtrn=new double[this.rowDimension()];
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=get(i,columnName);
		}
		return rtrn;
	}
	
	public double[] getRow(String rowName){
		double[] rtrn=new double[this.columnDimension()];
		//System.err.println(rowName);
		
		for(int j=0; j<rtrn.length; j++){
			rtrn[j]=get(rowName, j);
		}
		
		return rtrn;
	}
	
	public void set(int i, int j, double val) {data.set(i, j, val);}
	public void set(String row, String column, double value) {
		if(!rowIndexMap.containsKey(row)) {
			throw new IllegalArgumentException("Trying to add a value using a row ("+row+") not in the matrix");
		}
		
		if(!columnIndexMap.containsKey(column)) {
			throw new IllegalArgumentException("Trying to add a value using a column ("+column+") not in the matrix");
		}
		data.set(rowIndexMap.get(row), columnIndexMap.get(column), value);
	}
	
	public void setRow(String row, double[] vals){
		for(int i=0; i<vals.length; i++){
			set(row, i, vals[i]);
		}
	}
	
	public void setColumn(double[] vals, String column){
		for(int i=0; i<vals.length; i++){
			set(i, column, vals[i]);
		}
	}
	
	public void setColumn(double[] vals, int column){
		for(int i=0; i<vals.length; i++){
			set(i, column, vals[i]);
		}
	}
	
	public void set(String row, int colIdx, double value) {
		data.set(rowIndexMap.get(row), colIdx, value);
	}
	
	public void set(int rowIdx, String column, double value) {
		data.set(rowIdx, columnIndexMap.get(column), value);
		
	}
	
	public double get(String row, String column) {
		if(!rowIndexMap.containsKey(row) ) {
			//System.err.println("Row " + row + " not found");
		}
		
		if(!columnIndexMap.containsKey(column) ) {
			//System.err.println("Column " + column + " not found");
		}
		
		return data.get(rowIndexMap.get(row), columnIndexMap.get(column));
	}
	public boolean containsColumn (String colName) {
		return columnIndexMap.containsKey(colName);
	}
	
	public boolean containsRow (String rowName) {
		return rowIndexMap.containsKey(rowName);
	}
	public double get(String rowName, int colIdx) {return data.get(rowIndexMap.get(rowName), colIdx);}
	public double get(int rowIdx, String colName) {return data.get(rowIdx, columnIndexMap.get(colName));}
	
	public List<String> getColumnNames() { return new ArrayList<String>(columnIndexMap.keySet());}
	public List<String> getRowNames() { return new ArrayList<String>(rowIndexMap.keySet());}
	public String getRowName(int i){return rows.get(i);}
	public String getColoumnName(int i){return columns.get(i);}
	public List<String> getRowDescriptions() { return new ArrayList<String>(rowDescrIndexMap.keySet());}
	public boolean hasColumn(String column) { return columnIndexMap.containsKey(column);}
	public boolean hasRow(String row) { return rowIndexMap.containsKey(row);}

	public List<Integer> getIndecesForRowDescription(String rowDescription) {
		return rowDescrIndexMap.get(rowDescription);
	}

	public void setRowDescription(int rowIdx, String rowDescription) {
		List<Integer> descrIdxs =  rowDescrIndexMap.get(rowDescription);
		if(descrIdxs == null) {
			descrIdxs = new ArrayList<Integer>();
			rowDescrIndexMap.put(rowDescription, descrIdxs);
		}
		descrIdxs.add(rowIndexMap.get(rowIdx));
	}
	
	public void setRowDescription(String row, String rowDescription) {
		setRowDescription(rowIndexMap.get(row), rowDescription);
	}
	
	
	
	
	
	
	
	
	
	
	public MatrixWithHeaders submatrixByRowNames(Collection<String> rowNames) {
		ArrayList<String> names=new ArrayList<String>();
		
		for(String rowName: rowNames){
			if(this.rowIndexMap.containsKey(rowName)){names.add(rowName);}
		}
		// return null if no names are found
		if (names.isEmpty())
			return null;
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(names, getColumnNames());
		if(this.pidToName!=null){rtrn.setPIDToName(this.pidToName);}
		
		for(String rowName: rowNames){
			if(this.rowIndexMap.containsKey(rowName)){rtrn.setRow(rowName, getRow(rowName));}
		}
		
		return rtrn;
	}
	
	public MatrixWithHeaders submatrixByRowNames(String rowName) {
		ArrayList<String> s= new ArrayList<String>();
		s.add(rowName);
		return (submatrixByRowNames(s));
	}
	public MatrixWithHeaders submatrixByRowNames(String [] rowNames) {
		List<String> rowList = new ArrayList<String>(rowNames.length);
		for(String row: rowNames) {
			rowList.add(row);
		}

		return submatrixByRowNames(rowList);
	}
	
	public MatrixWithHeaders submatrixByColumnNames(Collection<String> columnNames) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(getRowNames(), new ArrayList<String>(columnNames));
		if(this.pidToName!=null){rtrn.setPIDToName(this.pidToName);}
		
		for(String columnName: columnNames){
			//System.err.println(columnName);
			double[] array=getColumn(columnName);
			rtrn.setColumn(array, columnName);
		}
		
		return rtrn;
	}
	
	public MatrixWithHeaders submatrixByColumnIndex(Collection<Integer> columnIndex) {
		ArrayList<String> columnNames=new ArrayList();
		for(Integer index: columnIndex){columnNames.add(this.getColumnNames().get(index));}
		
		return this.submatrixByColumnNames(columnNames);
	}
	
	public MatrixWithHeaders submatrixByColumnNames(String [] columnNames) {
		List<String> columnList = new ArrayList<String>(columnNames.length);
		for(String col: columnNames) {
			columnList.add(col);
		}
		return submatrixByColumnNames(columnList);
	}
	
	
	
	
	
	
		
	
	
	
	public void append(MatrixWithHeaders other) {
		if(other == null || other.rowDimension() != rowDimension()) {
			throw new IllegalArgumentException ("To append a matrix both matrices must have same number of rows");
		}
		
		Matrix newData = new Matrix(rowDimension(), columnDimension() + other.columnDimension());
		List<String> columnNames = getColumnNames();
		List<String> rowNames    = getRowNames();
		columnNames.addAll(other.getColumnNames());
		for(int i = 0; i < rowDimension(); i++) {
			for(int j = 0; j < columnDimension(); j++) {
				newData.set(i, j, data.get(i,j));
			}
			for(int j = 0; j < other.columnDimension(); j++) {
				newData.set(i, j + columnDimension(), other.get(i,j));
			}
		}
		this.data = newData;
		initNameIndexMaps(rowNames, columnNames);
	}
	
	public void appendColumns(MatrixWithHeaders other){
		append(other);
	}
	
	public void appendRows(MatrixWithHeaders other) {
		if(other == null || other.columnDimension() != columnDimension()) {
			throw new IllegalArgumentException ("To append a matrix both matrices must have same number of columns");
		}
		
		Matrix newData = new Matrix(rowDimension() + other.rowDimension(), columnDimension());
		List<String> columnNames = getColumnNames();
		List<String> rowNames    = getRowNames();
		rowNames.addAll(other.getRowNames());
		for(int j = 0; j < columnDimension(); j++) {
			for(int i = 0; i < rowDimension(); i++) {
				newData.set(i, j, data.get(i,j));
			}
			for(int i = 0; i < other.rowDimension(); i++) {
				newData.set(i + rowDimension(), j, other.get(i,j));
			}
		}
		this.data = newData;
		initNameIndexMaps(rowNames, columnNames);
	}
	
	public void add(MatrixWithHeaders m2) {
		data.plusEquals(m2.data);
	}
	
	
	public void log10() {
		log10(0);
	}
	
	public void log10(double fudge) {
		for(int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				data.set(i, j, Math.log10(data.get(i,j) + fudge));
			}
		}
	}
	
	public void pow() {
		for(int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				data.set(i, j, Math.pow(data.get(i,j),2));
			}
		}
	}
	
	public void round() {
		for(int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				data.set(i, j, Math.round(data.get(i,j)));
			}
		}
	}
	

	
	
	public void minus(MatrixWithHeaders other) {
		if(data.getColumnDimension() != other.columnDimension() || data.getRowDimension() != other.rowDimension()) {
			throw new IllegalArgumentException ("Trying to substract non compatible matrices  this has dimension " + rowDimension() + "x" + columnDimension() + "  other is " + other.rowDimension() + "x" +other.columnDimension());
		}
		
		data.minusEquals(other.data);
	}
	
	
	public Matrix getData() { return data;}
	
	protected void initFromRegularMatrix(BufferedReader br, String header)
			throws IOException {
		String [] columnNames = header.split("\t");
		List<String> columnNameList = new ArrayList<String> (columnNames.length);
		List<String> rowNameList = new ArrayList<String>();
		
		for(int i = 1; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
		}
		
		String line = null;
		List<List<Double>> rawData = new ArrayList<List<Double>>();
		int lineNum = 1;
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			if(info.length != columnNames.length) {
				throw new IllegalArgumentException("Line " + lineNum + " has " + info.length + " columns but header had " + columnNames.length + " columns");
			}
			rowNameList.add(info[0]);
			List<Double> lineData = new ArrayList<Double>(info.length - 1);
			rawData.add(lineData);
			for(int i = 1 ; i < info.length; i++) {
				lineData.add(Double.parseDouble(info[i]));
			}
			lineNum++;
		}
		
		init(rawData, rowNameList, columnNameList);
	}
	
	protected void initGCTFromReader(BufferedReader br) throws IOException {
		String line = br.readLine();
		if(line.startsWith("#1.2")) {
			line = br.readLine();
		}
		String [] dimensionsStr = line.split("\t");
		int expectedRowDimension = Integer.parseInt(dimensionsStr[0]);
		int expectedColumnDimension = Integer.parseInt(dimensionsStr[1]);
		
		String header = br.readLine();
		String [] columnNames = header.split("\t");
		List<String> columnNameList = new ArrayList<String> (columnNames.length - 2);
		if(columnNames.length - 2 != expectedColumnDimension) {
			System.err.println("WARNING: expected "+expectedColumnDimension+ " columns but read " +  (columnNames.length - 2) );
		}
		List<String> rowNameList = new ArrayList<String>(expectedRowDimension);
		List<String> rowDescrList = new ArrayList<String>(expectedRowDimension);
		
		for(int i = 2; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
		}
		
		
		this.pidToName=new TreeMap<String, String>();
		List<List<Double>> rawData = new ArrayList<List<Double>>(expectedRowDimension);
		int lineNum = 0;
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			if(info.length != columnNames.length) {
				throw new IllegalArgumentException("Line " + lineNum + " has " + info.length + " columns but header had " + columnNames.length + " columns");
			}
			rowNameList.add(info[0].toUpperCase().intern());
			rowDescrList.add(info[1].toUpperCase().intern());
			this.pidToName.put(info[0].toUpperCase(), info[1].toUpperCase());
			List<Double> lineData = new ArrayList<Double>(info.length - 2);
			rawData.add(lineData);
			for(int i = 2 ; i < info.length; i++) {
				lineData.add(Double.parseDouble(info[i]));
			}
			if(lineNum % 1000 == 0) {System.out.println("Line  " +  lineNum + " Used Mem " + ((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024)  );}
			lineNum++;
		}
		if(lineNum != expectedRowDimension) {
			System.err.println("WARNING: Expected " + expectedRowDimension + " but  read " + lineNum);
		}
		init(rawData, rowNameList, columnNameList, rowDescrList);
	}
	
	
	protected void initFromTxt(File file) throws IOException {
		initFromTxt(file, "\t");
	}
	
	protected void initFromTxt(BufferedReader br) throws IOException {
		initFromTxt(br, "\t", 0);
	}


	protected void initFromTxt(File file, String sep) throws IOException {
		BufferedReader br=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		initFromTxt(br, sep, 0);
	}
	
	protected void initFromTxt(BufferedReader br, String sep, int skip) throws IOException {
		for(int i=0; i<skip; i++) {br.readLine();}
		
		String line = br.readLine();
		
		this.pidToName=new TreeMap<String, String>();
		String header = line;
		String [] columnNames = header.split(sep); //Cells
		List<String> columnNameList = new ArrayList<String> (columnNames.length - 1);
		
		List<String> rowNameList = new ArrayList<String>();
		List<String> rowDescrList = new ArrayList<String>();
		
		for(int i = 1; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
		}
		
		
		List<List<Double>> rawData = new ArrayList<List<Double>>();
		int lineNum = 0;
		while( (line = br.readLine()) != null) {
			String [] info = line.split(sep);
				
				String name=info[0];
				rowNameList.add(name);
				rowDescrList.add(name);
				this.pidToName.put(name, name);
				
				List<Double> lineData = new ArrayList<Double>(info.length - 1);
				rawData.add(lineData);
				for(int i = 1 ; i < info.length; i++) {
					lineData.add(Double.parseDouble(info[i]));
				}
				//if(lineNum % 1000 == 0) {System.out.println("Line  " +  lineNum + " Used Mem " + ((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024)  );}
			
			lineNum++;
			if(lineNum%100000==0) {System.err.println(lineNum);}
		}
		br.close();
		init(rawData, rowNameList, columnNameList, rowDescrList);
	}
	
	protected void initFromTxt(BufferedReader br, String sep, int skip, int timeShift) throws IOException {
		for(int i=0; i<skip; i++) {br.readLine();}
		
		String line = br.readLine();
		
		this.pidToName=new TreeMap<String, String>();
		String header = line;
		String [] columnNames = header.split(sep); //Cells
		List<String> columnNameList = new ArrayList<String> (columnNames.length - 1);
		
		List<String> rowNameList = new ArrayList<String>();
		List<String> rowDescrList = new ArrayList<String>();
		
		for(int i = 1; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
		}
		
		
		List<List<Double>> rawData = new ArrayList<List<Double>>();
		int lineNum = 0;
		while( (line = br.readLine()) != null) {
			String [] info = line.split(sep);
				
				int val=Integer.parseInt(info[0])-timeShift;
				String name=new Integer(val).toString();
				rowNameList.add(name);
				rowDescrList.add(name);
				this.pidToName.put(name, name);
				
				List<Double> lineData = new ArrayList<Double>(info.length - 1);
				rawData.add(lineData);
				for(int i = 1 ; i < info.length; i++) {
					lineData.add(Double.parseDouble(info[i]));
				}
				//if(lineNum % 1000 == 0) {System.out.println("Line  " +  lineNum + " Used Mem " + ((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024)  );}
			
			lineNum++;
			if(lineNum%100000==0) {System.err.println(lineNum);}
		}
		br.close();
		init(rawData, rowNameList, columnNameList, rowDescrList);
	}
	
	
	
	protected void initFromTxt(File file, Map<String, String> nameToChr) throws IOException {
		BufferedReader br=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String line = br.readLine();
		
		this.pidToName=new TreeMap<String, String>();
		String header = line;
		String [] columnNames = header.split("\t"); //Cells
		List<String> columnNameList = new ArrayList<String> (columnNames.length - 1);
		
		List<String> rowNameList = new ArrayList<String>();
		List<String> rowDescrList = new ArrayList<String>();
		
		for(int i = 1; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
		}
		
		
		List<List<Double>> rawData = new ArrayList<List<Double>>();
		int lineNum = 0;
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			
			String name=info[0];
			String chr=nameToChr.get(name);
			if(chr==null){chr=name;}
			rowNameList.add(name);
			rowDescrList.add(chr);
			this.pidToName.put(name, chr);
			
			List<Double> lineData = new ArrayList<Double>(info.length - 1);
			rawData.add(lineData);
			for(int i = 1 ; i < info.length; i++) {
				lineData.add(Double.parseDouble(info[i]));
			}
			if(lineNum % 1000 == 0) {System.out.println("Line  " +  lineNum + " Used Mem " + ((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024)  );}
			lineNum++;
		}
		
		init(rawData, rowNameList, columnNameList, rowDescrList);
	}
	
	public Map<String, Collection<Integer>> getNanostringReplicateMap(){return this.replicateMap;}
	
	protected void initFromNanostring(BufferedReader br, String header) throws IOException {
		this.replicateMap=new TreeMap<String, Collection<Integer>>();
		
		String line = header;
		System.err.println(line);
		String [] columnNames = line.split("\t");
		List<String> columnNameList = new ArrayList<String> ();
		
		List<String> rowNameList = new ArrayList<String>();
		List<String> rowDescrList = new ArrayList<String>();
		
		for(int i = 3; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
			Collection<Integer> indexes=new ArrayList<Integer>();
			if(replicateMap.containsKey(columnNames[i])){indexes=replicateMap.get(columnNames[i]);}
			indexes.add(i-3);
			this.replicateMap.put(columnNames[i], indexes);
		}
		
		
		this.pidToName=new TreeMap<String, String>();
		this.probeClass=new TreeMap<String, String>();
		List<List<Double>> rawData = new ArrayList<List<Double>>();
		int lineNum = 0;
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			rowNameList.add(info[1]);
			rowDescrList.add(info[2]);
			this.pidToName.put(info[1], info[2]);
			this.probeClass.put(info[1], info[0]);
			
			List<Double> lineData = new ArrayList<Double>();
			rawData.add(lineData);
			for(int i = 3 ; i < info.length; i++) {
				lineData.add(Double.parseDouble(info[i]));
			}
			lineNum++;
		}
		
		init(rawData, rowNameList, columnNameList, rowDescrList);
	}
	
	public Map<String, String> getNanostringProbeClasses(){return this.probeClass;}

	private void init(List<List<Double>> rawData, List<String> rowNameList, List<String> columnNameList, List<String> rowDescriptionList ) {
		int m = rowNameList.size();
		int n = columnNameList.size();
		
		data = new Matrix(m,n);
		for(int i = 0; i < m; i++){
			for(int j = 0; j < n; j++) {
				data.set(i, j, rawData.get(i).get(j));
			}
		}
		initNameIndexMaps( rowNameList, columnNameList, rowDescriptionList);
	}
	
	private void init(List<List<Double>> rawData, List<String> rowNameList, List<String> columnNameList ) {
		init(rawData, rowNameList, columnNameList, new ArrayList<String>());
	}

	protected void initNameIndexMaps(List<String> rowNameList,List<String> columnNameList) {
		initNameIndexMaps(rowNameList, columnNameList, new ArrayList<String>());
	}
	
	protected void setData(Matrix data) {
		this.data = data;
	}
	
	protected void initRowIndexMaps(List<String> rowNameList) {
		int m = rowNameList.size();
		rowIndexMap    = new LinkedHashMap<String, Integer>(m);
		HashMap<String, Integer> occurrence = new HashMap<String, Integer>();
		occurrence = new HashMap<String, Integer>();
		for(int i = 0; i < rowNameList.size(); i ++) {
			String rowName = rowNameList.get(i);
			String uniqueRowName = getUniqueName(occurrence, rowName);
			rowIndexMap.put(uniqueRowName, i);
		}
	}
	
	protected void initColIndexMaps(List<String> columnNameList) {
		int n = columnNameList.size();
		columnIndexMap = new LinkedHashMap<String, Integer>(n);
		HashMap<String, Integer> occurrence = new HashMap<String, Integer>();
		for(int i = 0; i < columnNameList.size(); i ++) {
			String colName = columnNameList.get(i);
			String uniqueColName = getUniqueName(occurrence, colName);
			columnIndexMap.put(uniqueColName, i);
		}
	}
	
	protected void initRowDescrIndexMaps(List<String> rowDescriptionList) {
		rowDescrIndexMap = new LinkedHashMap<String, List<Integer>>();
	
		for(int i = 0; i < rowDescriptionList.size(); i ++) {
			String rowDescr = rowDescriptionList.get(i).toUpperCase();
			List<Integer> rowDescrIndeces = rowDescrIndexMap.get(rowDescr);
			if(rowDescrIndeces == null) {
				rowDescrIndeces = new ArrayList<Integer>();
				rowDescrIndexMap.put(rowDescr, rowDescrIndeces);
			}
			rowDescrIndeces.add(i);
		}
	}
	
	private void initNameIndexMaps(List<String> rowNameList,List<String> columnNameList, List<String> rowDescriptionList) {
		initRowIndexMaps(rowNameList);
		initColIndexMaps(columnNameList);
		initRowDescrIndexMaps(rowDescriptionList);
		rows = rowNameList;
		columns = columnNameList;
		
	}
	private String getUniqueName(HashMap<String, Integer> occurrence, String key) {
		Integer colNameOccurrence = occurrence.get(key) ;
		String uniqueKey = key;
		if(colNameOccurrence == null){
			occurrence.put(key,1);
		} else {
			occurrence.put(key,colNameOccurrence++);
			uniqueKey = getUniqueName(occurrence, key + "_" + colNameOccurrence);
		}
		
		return uniqueKey;
	}

	public int getNumberColumns(){return this.columnDimension();}
	public int getNumberRows(){return this.rowDimension();}
	
	public Map<String, double[]> toMap(){
		Map<String, double[]> rtrn=new TreeMap();
		
		for(String geneName: this.getRowNames()){
			double[] array=this.getRow(geneName);
			//System.err.println(geneName+" "+Statistics.average(array)+" "+array.length);
			rtrn.put(geneName, array);
		}
		
		return rtrn;
	}



	public void writeGMT(String save, double fold) throws IOException {
		FileWriter writer=new FileWriter(save);
		System.err.println("here");
		Collection<String> columnNames=getColumnNames();
		
		for(String name: columnNames){
			List<String> rowNames=getRowNames();
			System.err.println(name);
			double[] vals=getColumn(name);
			writer.write(name+"_UP"+"\t"+name+"_UP");
			for(int i=0; i<vals.length; i++){
				if(vals[i]>fold){writer.write("\t"+rowNames.get(i));}
			}
			writer.write("\n");
			
			writer.write(name+"_DOWN"+"\t"+name+"_DOWN");
			for(int i=0; i<vals.length; i++){
				if(vals[i]<(-fold)){writer.write("\t"+rowNames.get(i));}
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	//Moran's version
	/*public MatrixWithHeaders medianNorm(){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(this.rows, this.columns);
		
		for(int i=0; i<data.getColumnDimension(); i++){
			double median=Statistics.median(data.getColumn(i));
			double[] vals=medianNormColumn(data.getColumn(i), median);
			rtrn.setColumn(vals, i);
			System.err.println(this.columns.get(i)+" "+median+" "+this.getColumn(i).length);
		}
		
		return rtrn;
	}*/
	
	private double[] shiftArrayByConstant(double[] column, double median) {
		double[] rtrn=new double[column.length];
		for(int i=0; i<column.length; i++){
			rtrn[i]=column[i]-median;
		}
		return rtrn;
	}

	
	
	public MatrixWithHeaders copy(){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(this.rows, this.columns);
		
		for(int i=0; i<this.columnDimension(); i++){
			rtrn.setColumn(this.getColumn(i), i);
		}
		
		return rtrn;
	}


	public void setPIDToName(Map<String, String> rowDescriptions){
		if(rowDescriptions==null){
			System.err.println("Row description is null");
			return;
		}
		Map<String, String> allUpper=new TreeMap<String, String>();
		for(String pid: rowDescriptions.keySet()){
			String name=rowDescriptions.get(pid);
			allUpper.put(pid, name);

			List<Integer> descriptionRowIdxs = rowDescrIndexMap.get(name);
			if(descriptionRowIdxs == null) {
				descriptionRowIdxs = new ArrayList<Integer>();
				rowDescrIndexMap.put(name, descriptionRowIdxs);
			}
			descriptionRowIdxs.add(rowIndexMap.get(pid));
		}
		
		this.pidToName=allUpper;
	}


	public void writeCLS(String string, Collection<String> subset, Collection<String> negatives) throws IOException {
		FileWriter writer=new FileWriter(string);
		writer.write((subset.size()+negatives.size())+"\t2\t1\n");
		writer.write("# sample\tcontrol\n");
		for(String sample: subset){writer.write("0\t");}
		writer.write("\n");
		for(String sample: negatives){writer.write("1\t");}
		writer.write("\n");
		writer.write("control\n");
		writer.close();
	}
	
	
	
   

	//Only works is the matrix was defines to add rows
	public void addRow(String row, String description) {
			
		int index=this.getNumberRows();
		rowIndexMap.put(row, index);
			
		List<Integer> rowDescrIndeces = rowDescrIndexMap.get(description);
		if(rowDescrIndeces == null) {
			rowDescrIndeces = new ArrayList<Integer>();
			rowDescrIndexMap.put(description, rowDescrIndeces);
		}
		rowDescrIndeces.add(index);
			
		for(int i=0; i<this.getNumberColumns(); i++){
			set(row, i, 0.0);
		}
			
	}
	

	//Only works is the matrix was defines to add rows
	public void addRow(String row, String description, double[] vals) {
		
		int index=this.getNumberRows();
		rowIndexMap.put(row, index);
		
		List<Integer> rowDescrIndeces = rowDescrIndexMap.get(description);
		if(rowDescrIndeces == null) {
			rowDescrIndeces = new ArrayList<Integer>();
			rowDescrIndexMap.put(description, rowDescrIndeces);
		}
		rowDescrIndeces.add(index);
		
		for(int i=0; i<vals.length; i++){
			set(row, i, vals[i]);
		}
		
	}

	public void addColumn(String column) {
		
		Matrix newData = new Matrix(rowDimension(), columnDimension()+1);
		for(int i=0; i< rowDimension(); i++){
			for(int j = 0; j < columnDimension(); j++) {
				newData.set(i, j, get(i,j));
			}
		}
		
		setData(newData);	
		int index=columns.size();
		columnIndexMap.put(column, index);
		columns.add(column);
	}

	
	
	public void addRow(String column) {
		
		Matrix newData = new Matrix(rowDimension()+1, columnDimension());
		for(int i=0; i< rowDimension(); i++){
			for(int j = 0; j < columnDimension(); j++) {
				newData.set(i, j, get(i,j));
			}
		}
		
		setData(newData);	
		int index=columns.size();
		columnIndexMap.put(column, index);
		columns.add(column);
	}

	


	public void writeBox(String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(String rowName : rowIndexMap.keySet()) {
			for(String colName : columnIndexMap.keySet()) {
				writer.write(String.valueOf(get(rowName, colName))+"\t");
			}
			writer.write("\n");
		}
		
		writer.close();
	}


	public MatrixWithHeaders excludeByRowNames(Collection<String> flaggedGenes) {
		Collection<String> include=new TreeSet();
		
		
		for(String row: this.getRowNames()){
			if(!flaggedGenes.contains(row)){include.add(row);}
		}
		
		return this.submatrixByRowNames(include);
	}

	
	public MatrixWithHeaders excludeByColumnNames(Collection<String> flaggedColumns) {
		Collection<String> include=new TreeSet<String>();
		
		
		for(String row: this.getColumnNames()){
			if(!flaggedColumns.contains(row)){include.add(row);}
		}
		
		return this.submatrixByColumnNames(include);
	}

	
	
	/**
	 * Resets all values below floor to floor.
	 * @param floor
	 */
	public void floor(double floor) {
		for(int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				if(get(i,j) < floor) {
					set(i,j, floor);
				}
			}
		}
	}


	public void writeRowNamesToFile(String outFile) throws IOException {
		
		FileWriter writer=new FileWriter(outFile);
		for (int i=0; i< this.getNumberRows(); i++ )
			{writer.write(this.getRowName(i)+"\n");}
		writer.close();
		
	}
    
 public void writeColoumnNamesToFile(String outFile) throws IOException {
		
		FileWriter writer=new FileWriter(outFile);
		for (int i=0; i< this.getNumberColumns(); i++ )
			{writer.write(this.getColoumnName(i)+"\n");}
		writer.close();
		
	}


public double[] getRow(String gene, Collection<String> group1) {
	double[] rtrn=new double[group1.size()];
	
	int i=0;
	for(String column: group1){
		rtrn[i]=get(gene, column);
		i++;
	}
	
	return rtrn;
}


public double[] getValues(String gene, Collection<String> controls) {
	double[] rtrn=new double[controls.size()];
	
	int i=0;
	for(String control: controls){
		rtrn[i++]=get(gene, control);
	}
	
	return rtrn;
}


public boolean hasNanostringProbeClasses() {
	if(this.getNanostringProbeClasses()==null || this.getNanostringProbeClasses().isEmpty()){return false;}
	return true;
}

/**
 * @author skadri
 * Multiplies each column with a separate constant. That is, given a vector (dimension same as number of columns) multiplies all rows of column with same constant 
 * @return
 */
public MatrixWithHeaders multiplyColumnsWithConstants(double[] constants){
	
	MatrixWithHeaders resultMat = new MatrixWithHeaders(this.getRowNames(),this.getColumnNames());
	if(resultMat.columnDimension() != constants.length) {
		throw new IllegalArgumentException ("Trying to multiply non compatible matrix and vector. Columns on matrix " + resultMat.columnDimension() + " Vector Dimensions on right " + constants.length);
	}
	else{

		for(int j=0;j<resultMat.columnDimension();j++){
			for(int i=0;i<resultMat.rowDimension();i++){
				resultMat.set(i, j, (this.get(i,j)*constants[j]));
			}
		}
	}
	return resultMat;
}


/**
 * @author skadri
 * Divides each column with a separate constant. That is, given a vector (dimension same as number of columns) multiplies all rows of column with same constant 
 * @return
 */
public MatrixWithHeaders divideColumnsWithConstants(double[] constants){
	
	MatrixWithHeaders resultMat = new MatrixWithHeaders(this.getRowNames(),this.getColumnNames());
	if(resultMat.columnDimension() != constants.length) {
		throw new IllegalArgumentException ("Trying to multiply non compatible matrix and vector. Columns on matrix " + resultMat.columnDimension() + " Vector Dimensions on right " + constants.length);
	}
	else{

		for(int j=0;j<resultMat.columnDimension();j++){
			for(int i=0;i<resultMat.rowDimension();i++){
				resultMat.set(i, j, (this.get(i,j)/constants[j]));
			}
		}
	}
	return resultMat;
}


/**
 * @author skadri
 * Multiplies each column with a separate constant. That is, given a vector (dimension same as number of columns) multiplies all rows of column with same constant 
 * @return
 */
public MatrixWithHeaders multiplyColumnsWithConstants(Map<String,Double> constants){
	
	MatrixWithHeaders resultMat = new MatrixWithHeaders(this.getRowNames(),this.getColumnNames());
	if(resultMat.columnDimension() != constants.keySet().size()) {
		throw new IllegalArgumentException ("Trying to multiply non compatible matrix and vector. Columns on matrix " + resultMat.columnDimension() + " Vector Dimensions on right " + constants.size());
	}
	else{

		for(int j=0;j<resultMat.columnDimension();j++){
			double constant = constants.get(resultMat.getColoumnName(j));
			for(int i=0;i<resultMat.rowDimension();i++){
				double value = constant*this.get(i, j);
				//System.out.println(constants.get(resultMat.getColoumnName(j))+" * "+ this.get(i, j));
				resultMat.set(i, j, value);
			}
		}
	}
	return resultMat;
}


public MatrixWithHeaders transpose() {
	MatrixWithHeaders rtrn=new MatrixWithHeaders(this.getColumnNames(), this.getRowNames());
	
	for(String row: this.getRowNames()){
		rtrn.setColumn(this.getRow(row), row);
	}
	
	return rtrn;
}


public void incrementCount(String row, String column) {
	double score=get(row, column);
	score++;
	set(row, column,score);	
}

public void incrementCount(String row, String column, double increment) {
	double score=get(row, column);
	score+=increment;
	set(row, column,score);	
}












	
	


	
	
	
}
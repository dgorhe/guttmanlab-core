package guttmanlab.core.annotation.io;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BEDFileRecord;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.FeatureCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

public class BEDFileIO implements AnnotationFileIO<BEDFileRecord> {

	private CoordinateSpace referenceSpace;
	
	/**
	 * @param referenceSizes Reference coordinate space
	 */
	public BEDFileIO(CoordinateSpace referenceSizes) {
		this.referenceSpace = referenceSizes;
	}
	
	/**
	 * @param fileName The BED file
	 */
	public BEDFileIO(String fileName) {
		this(new CoordinateSpace(fileName));
	}
	
	/**
	 * Write features to a file
	 * @param regions Features to write
	 * @param outputBedFile Output file
	 * @throws IOException
	 */
	public static void writeToFile(AnnotationCollection<? extends Annotation> regions, File outputBedFile) throws IOException {
		try (FileWriter w = new FileWriter(outputBedFile)) {
			Iterator<? extends Annotation> iter = regions.sortedIterator();
			while(iter.hasNext()) {
				w.write((new BEDFileRecord(iter.next()).toFormattedString()));
			}	
		}
	}
	
	/**
	 * Write features to a file
	 * @param regions Features to write
	 * @param outputBedFile Output file name
	 * @throws IOException
	 */
	public static void writeToFile(AnnotationCollection<? extends Annotation> regions, String outputBedFile) throws IOException {
		writeToFile(regions, new File(outputBedFile));
	}
	
	/**
	 * Loads the contents of a BED file into memory. The contents are organized by reference name.
	 * @param file is the BED file to open 
	 * @param refSpace is the reference coordinate space which contains reference names and sizes
	 * @return A map from reference name to features on that reference described in the BED file
	 * @throws IOException
	 */
	public static Map<String, FeatureCollection<BEDFileRecord>> loadFromFileByReferenceName(File file, CoordinateSpace refSpace) throws IOException {
		Map<String, FeatureCollection<BEDFileRecord>> rtrn = new TreeMap<String, FeatureCollection<BEDFileRecord>>();
		try (BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)))) {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				BEDFileRecord annot = BEDFileRecord.fromFormattedString(line);
				String refName = annot.getReferenceName();
				rtrn.putIfAbsent(refName, new FeatureCollection<BEDFileRecord>(refSpace));
				rtrn.get(refName).addAnnotation(annot);
			}
		}	
		return rtrn;
	}

	/**
	 * Loads the contents of a BED file into memory.
	 * @param inputBedFileName is the name of the BED file to read
	 * @param coordinateSpace is the reference coordinate information containing reference names and sizes
	 * @return a collection of BED file records
	 * @throws IOException
	 */
	public static AnnotationCollection<BEDFileRecord> loadFromFile(String inputBedFileName, CoordinateSpace coordinateSpace) throws IOException {
		BEDFileIO bfio = new BEDFileIO(coordinateSpace);
		return bfio.loadFromFile(new File(inputBedFileName));
	}
	
	/**
	 * Loads the contents of a BED file into memory.
	 * @param inputBedFile is the BED file to read
	 * @param coordinateSpace is the reference coordinate information containing reference names and sizes
	 * @return a collection of BED file records
	 * @throws IOException
	 */
	public static AnnotationCollection<BEDFileRecord> loadFromFile(File inputBedFile, CoordinateSpace coordinateSpace) throws IOException {
		BEDFileIO bfio = new BEDFileIO(coordinateSpace);
		return bfio.loadFromFile(inputBedFile);
	}
	
	/**
	 * Loads the contents of a BED file into memory.
	 * @param inputBedFileName is the name of the BED file to read
	 * @return a collection of BED file records
	 * @throws IOException
	 */
	public AnnotationCollection<BEDFileRecord> loadFromFile(String inputBedFileName) throws IOException {
		return loadFromFile(new File(inputBedFileName));
	}
	
	/**
	 * Loads the contents of a BED file into memory.
	 * @param inputBedFile is the BED file to read
	 * @return a collection of BED file records
	 * @throws IOException
	 */
	@Override
	public AnnotationCollection<BEDFileRecord> loadFromFile(File inputBedFile) throws IOException {
		FeatureCollection<BEDFileRecord> collection = new FeatureCollection<BEDFileRecord>(referenceSpace);
		try (BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(inputBedFile)))) {
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				BEDFileRecord annot = BEDFileRecord.fromFormattedString(line);
				collection.addAnnotation(annot);
			}
		}
		return collection;
	}
}
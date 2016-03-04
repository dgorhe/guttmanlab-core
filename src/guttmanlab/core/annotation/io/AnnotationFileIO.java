package guttmanlab.core.annotation.io;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotationcollection.AnnotationCollection;

import java.io.File;
import java.io.IOException;

/**
 * Methods to read from and write to standard annotation file formats
 * @author prussell
 *
 */
public interface AnnotationFileIO<T extends Annotation> {
	
	/**
	 * Read the file and get the annotations as annotation objects
	 * @param fileName File name
	 * @return The collection of annotations described in the file
	 * @throws IOException 
	 */
	public AnnotationCollection<T> loadFromFile(File fileName) throws IOException;
	
}

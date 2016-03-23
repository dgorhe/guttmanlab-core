package guttmanlab.core.annotationcollection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.util.Collection;

import net.sf.samtools.util.CloseableIterator;

import org.apache.commons.collections15.Predicate;

/**
 * A container of annotations
 * @author prussell
 *
 */
public interface AnnotationCollection<T extends Annotation> {
	
	/**
	 * Merge the annotations in this collection
	 * @return A new annotation collection consisting of the merged annotations
	 */
	public AnnotationCollection<T> merge();
	
	/**
	 * @param other An annotation
	 * @return True iff the annotation overlaps some annotation in this collection
	 */
	public boolean overlaps(Annotation other);
	
	/**
	 * Get the number of annotations in this collection that overlap the region
	 * @param region Region to check for overlappers
	 * @param fullyContained Whether the overlappers must be fully contained in the region in order to count
	 * @return The number of overlappers
	 */
	public int numOverlappers(Annotation region, boolean fullyContained);
	
	/**
	 * @return A sorted iterator over all the annotations in this collection
	 */
	public CloseableIterator<T> sortedIterator();
	
	/**
	 * A sorted iterator that collapses annotations with the same name
	 * Returns exactly one representative for each set of annotations sharing all attributes except name
	 * Any representative may be returned
	 * Never returns multiple annotations that differ by name only
	 * @return A sorted iterator over all the annotations in this collection, never returning
	 * multiple annotations that differ by name only
	 */
	public default CloseableIterator<T> sortedIteratorIgnoreName() {
		throw new UnsupportedOperationException();
	}
	
	
	/**
	 * A sorted iterator over the set of annotations in this collection that overlap the region (must be sorted)
	 * @param region Region to check for overlappers
	 * @param fullyContained Whether the overlappers must be fully contained in the region
	 * @return Iterator over the set of overlappers
	 */	
	public CloseableIterator<T> sortedIterator(Annotation region, boolean fullyContained);
	
	/**
	 * Add a filter that will be applied when getting subsets or iterators
	 * @param filter A predicate that is true if the annotation passes the filter (should be kept)
	 */
	public void addFilter(Predicate<T> filter);
	
	/**
	 * Get all the filters passed to this class
	 * @return A collection of all predicate filters
	 */
	public Collection<Predicate<T>> getFilters();
	
	
	/**
	 * Write the collection of annotations (using filters) to a file
	 * @param fileName The file to write to
	 */
	public void writeToBAM(String fileName);
	
	/**
	 * Write the reads overlapping this region to a file
	 * @param fileName The file to write to
	 * @param region The region to overlap
	 * @param fullyContained whether to use only reads fully contained in this region (true) or merely overlapping (false)
	 */
	public void writeToBAM(String fileName, Annotation region, boolean fullyContained);
	
	/**
	 * @return The reference coordinate space
	 */
	public CoordinateSpace getReferenceCoordinateSpace();

	/**
	 * 
	 * @return The feature coordinate space
	 */
	public CoordinateSpace getFeatureCoordinateSpace();
	
	/** TODO Consider returning an AnnotationCollection
	 * Convert the annotation from the reference space
	 * If the reference equals the reference of the FeatureSpace it will return feature space
	 * If the reference equals the feature in the FeatureSpace it will return reference space
	 * If the reference doesn't equal either, it will throw an Exception
	 * @param annotations An iterator of the annotations to convert
	 * @param referenceSpaceForAnnotations The reference coordinate space of the annotation to map (must match either reference or feature space)
	 * @param fullyContained whether to only convert a fragment if it is fully within the feature (true) or return the trimmed partial (false)
	 * @return An iterator of annotations in the new feature space (A window object but really can be anything with hasParent())
	 */
	public <X extends Annotation> AnnotationCollection<DerivedAnnotation<X>> convertCoordinates(AnnotationCollection<X> readCollection, CoordinateSpace referenceSpaceForAnnotations, boolean fullyContained);
		
	/**
	 * Get an iterator over all populated windows of size widowLength overlapping a region
	 * @param region The region to scan
	 * @param windowLength The length of the window
	 * @return Iterator of windows and all of their overlapping reads of type T
	 */
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(Annotation region, int windowLength);
	
	/**
	 * Get an iterator over all populated windows of size widowLength overlapping a region
	 * @param region The region to scan
	 * @param windowLength The length of the window
	 * @param stepSize Step size
	 * @return Iterator of windows and all of their overlapping reads of type T
	 */
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(Annotation region, int windowLength, int stepSize);
	
	/**
	 * Get an iterator over all populated windows of size widowLength
	 * @param windowLength The length of the window
	 * @param stepSize Step size
	 * @return Iterator of windows and all of their overlapping reads of type T
	 */
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(int windowLength, int stepSize);
	
	/**
	 * Get an iterator over all populated windows of size widowLength overlapping a region
	 * @param region The region to scan
	 * @param windowLength The length of the window
	 * @param stepSize Step size
	 * @param includeEmpties Also return windows with zero overlappers
	 * @return Iterator of windows and all of their overlapping reads of type T
	 */
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(Annotation region, int windowLength, int stepSize, boolean includeEmpties);
	
			
	/**
	 * Iterates through the annotations in sortedIterator()
	 * @return the number of annotations in the collection
	 */
	public int getNumAnnotations();
}

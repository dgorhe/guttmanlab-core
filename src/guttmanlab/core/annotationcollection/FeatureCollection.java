package guttmanlab.core.annotationcollection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.lang.builder.HashCodeBuilder;
import org.apache.log4j.Logger;

import sun.nio.cs.ext.TIS_620;
import net.sf.samtools.util.CloseableIterator;

public class FeatureCollection<T extends Annotation> extends AbstractAnnotationCollection<T> implements Collection<T> {
	
	/**
	 * A loose index that keeps track of one object of class T per window.
	 * This is to be used to get reduced views of the TreeSet of annotations
	 * with TreeSet methods like headSet, lower, etc.
	 * 
	 * The window size is specified as a parameter to the constructor.
	 * As the FeatureCollection is built up one feature at a time, the
	 * first feature to be completely contained in each window becomes
	 * the representative of that window.
	 * 
	 * Windows with no completely contained features are absent from the index,
	 * but that fact is invisible to client code.
	 * 
	 * 
	 * @author prussell
	 *
	 */
	private class Index {
		
		/**
		 * A class that represents windows of the index
		 * @author prussell
		 *
		 */
		private class Interval implements Comparable<Interval> {
			
			private String chr;
			private int start;
			private int end;
			
			public Interval(String chr, int start, int end) {
				this.chr = chr;
				this.start = start;
				this.end = end;
			}
			
			public boolean equals(Object o) {
				if(!o.getClass().equals(Interval.class)) return false;
				Interval oi = (Interval)o;
				return chr.equals(oi.chr) && start == oi.start && end == oi.end;
			}
			
			public int hashCode() {
				HashCodeBuilder builder = new HashCodeBuilder();
				builder.append(chr);
				builder.append(start);
				builder.append(end);
				return builder.toHashCode();
			}

			/**
			 * Overlapping intervals return 0
			 */
			@Override
			public int compareTo(FeatureCollection<T>.Index.Interval o) {
				int compareChr = chr.compareTo(o.chr);
				if(compareChr != 0) return compareChr;
				if(o.end < start) return 1;
				if(o.start > end) return -1;
				return 0;
			}
			
			public boolean completelyContains(Annotation annot) {
				return annot.getReferenceName().equals(chr) && annot.getReferenceStartPosition() > start && annot.getReferenceEndPosition() < end;
			}
			
		}
		
		/**
		 * This is the index. One representative of class T per interval.
		 */
		private TreeMap<Interval, T> oneCompletelyContained;
		
		/**
		 * The set of intervals that have a representative.
		 */
		private TreeSet<Interval> indexedIntervals;
		
		/**
		 * The set of intervals that do not have a representative.
		 */
		private Collection<Interval> noRepresentativeYet;
		
		/**
		 * Initialize a new index
		 * @param intervalLength Interval length
		 */
		public Index(int intervalLength) {
			oneCompletelyContained = new TreeMap<Interval, T>();
			noRepresentativeYet = new TreeSet<Interval>();
			indexedIntervals = new TreeSet<Interval>();
			for(String chr : referenceCoordinateSpace.getRefSeqLengths().keySet()) {
				boolean addedChr = false;
				int chrLen = referenceCoordinateSpace.getRefSeqLengths().get(chr);
				// Construct the intervals for the chromosome
				for(int pos = 0; pos < chrLen - intervalLength; pos += intervalLength) {
					noRepresentativeYet.add(new Interval(chr, pos, pos + intervalLength));
					addedChr = true;
				}
			}
		}
		
		/**
		 * Update the index with a new annotation
		 * @param annot Annotation
		 * @return True if the annotation was incorporated as a representative in the index
		 */
		public boolean update(T annot) {
			for(Interval interval : noRepresentativeYet) {
				if(interval.completelyContains(annot)) {
					noRepresentativeYet.remove(interval);
					oneCompletelyContained.put(interval, annot);
					indexedIntervals.add(interval);
					return true;
				}
			}
			return false;
		}
		
		/**
		 * Construct an Interval representing the span of an annotation
		 * @param annot Annotation
		 * @return Interval representing the span only
		 */
		private Interval fromAnnotation(Annotation annot) {
			return new Interval(annot.getReferenceName(), annot.getReferenceStartPosition(), annot.getReferenceEndPosition());
		}
		
		/**
		 * Get some annotation in the FeatureCollection that is strictly less than 
		 * the given annotation with respect to the compareTo method of class T
		 * @param annot Annotation
		 * @return Some annotation that is guaranteed to be less than annot, but not
		 * guaranteed to be the greatest element of the FeatureCollection that is less
		 * than annot
		 */
		public T getSomeLowerBound(Annotation annot) {
			Interval annotAsInterval = fromAnnotation(annot);
			Interval lowerBoundInterval = indexedIntervals.lower(annotAsInterval);
			if(lowerBoundInterval == null) return null;
			return oneCompletelyContained.get(lowerBoundInterval);
		}
		
		/**
		 * Get some annotation in the FeatureCollection that is strictly greater than 
		 * the given annotation with respect to the compareTo method of class T
		 * @param annot Annotation
		 * @return Some annotation that is guaranteed to be greater than annot, but not
		 * guaranteed to be the least element of the FeatureCollection that is greater
		 * than annot
		 */
		public T getSomeUpperBound(Annotation annot) {
			Interval annotAsInterval = fromAnnotation(annot);
			Interval upperBoundInterval = indexedIntervals.higher(annotAsInterval);
			if(upperBoundInterval == null) return null;
			return oneCompletelyContained.get(upperBoundInterval);
		}
		
		/**
		 * Get a subset of the FeatureCollection that is guaranteed to contain all
		 * overlappers of an annotation
		 * @param annot Annotation
		 * @return A subset of the FeatureCollection that contains all overlappers of
		 * annot, but may also contain non-overlappers
		 */
		public SortedSet<T> getSupersetOfOverlappers(Annotation annot) {
			T lowerBound = getSomeLowerBound(annot);
			T upperBound = getSomeUpperBound(annot);
			if(lowerBound == null && upperBound == null) return annotations;
			if(lowerBound == null) return annotations.headSet(upperBound, true);
			if(upperBound == null) return annotations.tailSet(lowerBound, true);
			return annotations.subSet(lowerBound, true, upperBound, true);
		}
		
		
		
	}
	
	/**
	 * The reference coordinate system that features are mapped to
	 */
	private CoordinateSpace referenceCoordinateSpace;
	private TreeSet<T> annotations;
	private int featureCount;
	private Index index;
	private static final int INDEX_INTERVAL_LENGTH = 500000;
	private static Logger logger = Logger.getLogger(FeatureCollection.class.getName());
	
	/**
	 * @param referenceSpace Reference coordinate space
	 */
	public FeatureCollection(CoordinateSpace referenceSpace){
		super();
		this.referenceCoordinateSpace=referenceSpace;
		this.index = new Index(INDEX_INTERVAL_LENGTH);
		annotations = new TreeSet<T>();
	}
	
	/**
	 * Add annotation to the collection
	 * @param annotation to add
	 * @return true iff the collection changed
	 */
	public boolean addAnnotation(T annotation){
		boolean rtrn = annotations.add(annotation);
		if(rtrn) {
			index.update(annotation);
			featureCount++;
		}
		return rtrn;
	}

	/**
	 * Get the number of features in this collection
	 * @return The number of features
	 */
	public int getCount() {
		return this.featureCount;
	}

	/**
	 * Write bed file
	 * @param fileName Output file
	 */
	public void writeToFile(String fileName) {
		CloseableIterator<T> iter=sortedIterator();
		try{writeToFile(fileName, iter);}catch(IOException ex){ex.printStackTrace();}
	}
	
	private void writeToFile(String fileName, CloseableIterator<T> iter) throws IOException{
		FileWriter writer=new FileWriter(fileName);
		while(iter.hasNext()){
			T next=iter.next();
			writer.write(next.toString()+"\n");
		}
		writer.close();
		iter.close();
	}

	@Override
	public CloseableIterator<T> sortedIterator() {
		return new FilteredIterator<T>(annotations.iterator(), getFilters());
	}

	@Override
	public CloseableIterator<T> sortedIterator(Annotation region, boolean fullyContained) {
		throw new UnsupportedOperationException("Not implemented");
		//TODO write this iterator
		//Iterator<T> iter = null;
		//return new FilteredIterator<T>(iter, getFilters());
	}
	
	@Override
	public FeatureCollection<T> merge() {
		CloseableIterator<T> old = sortedIterator();
		FeatureCollection<T> merged = new FeatureCollection<T>(referenceCoordinateSpace);
		T current = null;
		if(old.hasNext())
			current = old.next();
		while(old.hasNext())
		{
			T next = old.next();
			if(current.overlaps(next))
			{
				current = (T) current.merge(next);
				//System.out.println(current.toBED());
			}
			else
			{
				merged.add(current);
				current = next;
			}
		}
		if(current != null)
			merged.add(current);
		return merged;
	}
	/**
	 * @return The coordinateSpace of the reference for this annotation collection
	 */
	public CoordinateSpace getReferenceCoordinateSpace(){return this.referenceCoordinateSpace;}
	
	
	public class WrappedIterator implements CloseableIterator<T>{

		Map<String, IntervalTree<T>> annotationTree;
		Iterator<String> referenceIterator;
		String currentReference;
		Iterator<T> currentTreeIterator;
		
		public WrappedIterator(Map<String, IntervalTree<T>> annotationTree){
			this.annotationTree=annotationTree;
			this.referenceIterator=annotationTree.keySet().iterator();
		}
		
		@Override
		public boolean hasNext() {
			if(currentTreeIterator==null || !currentTreeIterator.hasNext()){
				if(referenceIterator.hasNext()){
					String chr=referenceIterator.next();
					currentTreeIterator=annotationTree.get(chr).valueIterator();
					return hasNext();
				}
				else{return false;}
			}
			return true;
		}

		@Override
		public T next() {
			return currentTreeIterator.next();
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub
			throw new UnsupportedOperationException();
		}

		@Override
		public void close() {
			// TODO Auto-generated method stub
			//throw new UnsupportedOperationException();
		}
	}

	@Override
	public int size() {
		return getCount();
	}

	@Override
	public boolean isEmpty() {
		return size() == 0;
	}

	@Override
	public boolean contains(Object o) {
		T annot = (T)o;
		return annotations.contains(annot);
	}

	@Override
	public Iterator<T> iterator() {
		return sortedIterator();
	}

	@Override
	public Object[] toArray() {
		throw new UnsupportedOperationException();
	}

	@Override
	public <T1> T1[] toArray(T1[] annotations) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean add(T annotation) {
		return addAnnotation(annotation);
	}
	
	
	@Override
	public boolean overlaps(Annotation other) {
		Collection<T> possibleOverlappers = index.getSupersetOfOverlappers(other);
		for(T annot : possibleOverlappers) {
			if(annot.overlaps(other)) return true;
		}
		return false;
	}
	
	/**
	 * Get all features in this collection that overlap another feature
	 * @param other Other feature
	 * @return Collection of overlappers from this collection
	 */
	public FeatureCollection<T> overlappers(Annotation other) {
		FeatureCollection<T> rtrn = new FeatureCollection<T>(referenceCoordinateSpace);
		index.getSupersetOfOverlappers(other).forEach(t -> {if(t.overlaps(other)) rtrn.add(t);});
		return rtrn;
	}
	
	@Override
	public boolean remove(Object o) {
		T annot = (T)o;
		boolean rtrn = annotations.remove(annot);
		if(rtrn) featureCount--;
		return rtrn;
	}

	@Override
	public boolean containsAll(Collection<?> annotations) {
		for(Object o : annotations) {
			if(!contains(o)) return false;
		}
		return true;
	}

	@Override
	public boolean addAll(Collection<? extends T> annotations) {
		boolean rtrn = false;
		for(T annotation : annotations) {
			boolean changed = add(annotation);
			if(changed) rtrn = true;
		}
		return rtrn;
	}

	@Override
	public boolean removeAll(Collection<?> annotations) {
		boolean rtrn = false;
		for(Object o : annotations) {
			boolean r = remove(o);
			if(r) rtrn = true;
		}
		return rtrn;
	}

	@Override
	public boolean retainAll(Collection<?> annotations) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void clear() {
		annotations.clear();
		featureCount = 0;
	}
	
}

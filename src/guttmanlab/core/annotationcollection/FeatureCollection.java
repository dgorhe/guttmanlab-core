package guttmanlab.core.annotationcollection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import sun.nio.cs.ext.TIS_620;
import net.sf.samtools.util.CloseableIterator;

public class FeatureCollection<T extends BlockedAnnotation> extends AbstractAnnotationCollection<T> implements Collection<T> {

	/**
	 * The reference coordinate system that features are mapped to
	 */
	private CoordinateSpace referenceCoordinateSpace;
	private Map<String, IntervalTree<T>> annotationTree;
	private int featureCount;
	
	public FeatureCollection(CoordinateSpace referenceSpace){
		super();
		this.referenceCoordinateSpace=referenceSpace;
		this.annotationTree=new TreeMap<String, IntervalTree<T>>();
	}
	
	/**
	 * Add annotation to the collection
	 * @param annotation to add
	 * @return true iff the collection changed
	 */
	public boolean addAnnotation(T annotation){
		boolean alreadyContains = contains(annotation);
		IntervalTree<T> tree=new IntervalTree<T>();
		if(annotationTree.containsKey(annotation.getReferenceName())){
			tree=annotationTree.get(annotation.getReferenceName());
		}
		tree.put(annotation.getReferenceStartPosition(), annotation.getReferenceEndPosition(), annotation);
		annotationTree.put(annotation.getReferenceName(), tree);
		featureCount++;
		return !alreadyContains;
	}

	/**
	 * Get the number of features in this collection
	 * @return The number of features
	 */
	public int getCount() {
		return this.featureCount;
	}

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
		return new FilteredIterator<T>(new WrappedIterator(this.annotationTree), getFilters());
	}

	@Override
	public CloseableIterator<T> sortedIterator(Annotation region, boolean fullyContained) {
		IntervalTree<T> tree=this.annotationTree.get(region.getReferenceName());
		if(tree==null)  //tree may be null if region's reference was not found in annotations CB
			return null;
		Iterator<T> iter=tree.overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
		return new FilteredIterator<T>(iter, getFilters());
	}
	
	/*@Override
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
	}*/
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

	public void writeToFile(String fileName, Annotation region) {
		try{writeToFile(fileName, sortedIterator(region, false));
		}catch(IOException ex){ex.printStackTrace();}
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
		//TODO: This cast causes exceptions; consider checking type first or changing signature.
		T annot = (T)o;
		String chr = annot.getReferenceName();
		if(!annotationTree.containsKey(chr)) {
			return false;
		}
		int start = annot.getReferenceStartPosition();
		int end = annot.getReferenceEndPosition();
		Node<T> node = annotationTree.get(chr).find(start, end);
		if(node == null) {
			return false;
		}
		return node.getContainedValues().contains(annot);
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
		if(!annotationTree.containsKey(other.getReferenceName())) {
			return false;
		}
		Iterator<T> overlappers = annotationTree.get(other.getReferenceName()).overlappingValueIterator(other.getReferenceStartPosition(), other.getReferenceEndPosition());
		while(overlappers.hasNext()) {
			T annot = overlappers.next();
			if(other.overlaps(annot)) {
				return true;
			}
		}
		return false;
	}
	
	@Override
	public boolean remove(Object o) {
		T annot = (T)o;
		String chr = annot.getReferenceName();
		int start = annot.getReferenceStartPosition();
		int end = annot.getReferenceEndPosition();
		Node<T> node = annotationTree.get(chr).find(start, end);
		if(node == null) {
			return false;
		}
		boolean rtrn = node.getContainedValues().remove(annot);
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
		annotationTree.clear();
		featureCount = 0;
		for(IntervalTree<T> tree : annotationTree.values()) {
			featureCount += tree.size();
		}
	}
	
}

package guttmanlab.core.annotationcollection;

import guttmanlab.core.annotation.Annotation;
import net.sf.samtools.util.CloseableIterator;

public class NameIgnoringIterator<T extends Annotation> implements CloseableIterator<T> {
	
	private CloseableIterator<T> sortedIter;
	
	private T next;
	
	private NameIgnoringIterator(CloseableIterator<T> iter) {
		if(!iter.hasNext()) {
			throw new IllegalArgumentException("Iterator is empty");
		}
		sortedIter = iter;
		next = sortedIter.next();
	}
	
	/**
	 * Factory method
	 * @param annotations Annotation collection
	 * @return Name ignoring iterator over the collection
	 */
	public static <T extends Annotation> NameIgnoringIterator<T> forAnnotations(AnnotationCollection<T> annotations) {
		return new NameIgnoringIterator<T>(annotations.sortedIterator());
	}
	
	@Override
	public boolean hasNext() {
		return next != null;
	}
	
	private void assertSorted(T first, T second) {
		if(first.compareToIgnoreName(second) > 0) {
			throw new IllegalStateException("\nAnnotations are not in sorted order:\n" + first.toBED() + "\n" + second.toBED());
		}
	}
	
	@Override
	public T next() {
		T rtrn = next;
		next = null;
		T prev = rtrn;
		while(sortedIter.hasNext()) {
			T curr = sortedIter.next();
			assertSorted(prev, curr);
			prev = curr;
			if(!curr.equalsIgnoreName(rtrn)) {
				next = curr;
				break;
			}
		}
		return rtrn;
	}

	@Override
	public void close() {
		sortedIter.close();
	}
		
}

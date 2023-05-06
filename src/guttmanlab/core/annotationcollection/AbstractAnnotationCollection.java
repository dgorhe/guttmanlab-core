package guttmanlab.core.annotationcollection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.ContiguousWindow;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.BlockedWindow;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.util.CloseableIterator;

public abstract class AbstractAnnotationCollection<T extends Annotation> implements AnnotationCollection<T>{

	private Collection<Predicate<T>> filters;
	private int numAnnotations;

	public AbstractAnnotationCollection(){
		filters=new ArrayList<Predicate<T>>();
	}

	@Override
	public void addFilter(Predicate<T> filter) {
		filters.add(filter);
		numAnnotations = 0;
	}
	
	@Override
	public void addFilter(Collection<Predicate<T>> filters) {
		Iterator<Predicate<T>> iter = filters.iterator();
		while(iter.hasNext())
			addFilter(iter.next());
	}

	@Override
	public Collection<Predicate<T>> getFilters(){
		return filters;
	}
	
	@Override
	public void clearFilters(){
		filters=new ArrayList<Predicate<T>>();
	}

	@Override
	public AnnotationCollection<T> merge() {
		// FIXME Auto-generated method stub
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public boolean overlaps(Annotation other) {
		// FIXME Auto-generated method stub
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public <X extends Annotation> AnnotationCollection<DerivedAnnotation<X>> convertCoordinates(AnnotationCollection<X> readCollection, CoordinateSpace referenceSpaceForAnnotations, boolean fullyContained){
		//TODO Should check the coordinate space and convert appropriately
		return convertFromReferenceSpace(readCollection, referenceSpaceForAnnotations, fullyContained);
	}

	private <X extends Annotation> AnnotationCollection<DerivedAnnotation<X>> convertFromReferenceSpace(AnnotationCollection<X> readCollection, CoordinateSpace referenceSpace, boolean fullyContained){
		return new ConvertedSpace<X>(readCollection, this, referenceSpace, fullyContained);
	}

	
	
	@Override
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(Annotation region, int windowLength){
		return getPopulatedWindows(region, windowLength, 1);
	}

	@Override
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(Annotation region, int winSize, int stepSize) {
		return getPopulatedWindows(region, winSize, stepSize, false);
	}	
	
	@Override
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(int winSize, int stepSize) {
		CloseableIterator<T> iter=sortedIterator();
		return new WindowIterator<T>(iter,winSize,stepSize);
	}	
	
	@Override
	public CloseableIterator<? extends PopulatedWindow<T>> getPopulatedWindows(Annotation region, int winSize, int stepSize, boolean includeEmpties) {
		CloseableIterator<T> iter=sortedIterator(region, false);
		return new WindowRegionIterator<T>(region, iter,winSize,stepSize);
	}
	
	@Override
	public int numOverlappers(Annotation region, boolean fullyContained) {
		int counter=0;
		CloseableIterator<T> iter=sortedIterator(region, fullyContained);
		try {
			while(iter.hasNext()){
				iter.next();
				counter++;
			}
			iter.close();
			return counter;
		} catch(Exception e) {
			iter.close();
			throw (e);
		}
	}
	
	
	
	
	

	public class WindowRegionIterator<T1 extends Annotation> implements CloseableIterator<PopulatedWindow<T1>>{

		CloseableIterator<PopulatedWindow<T1>> windowIter;
		
		
		public WindowRegionIterator(Annotation region, CloseableIterator<T1> reads, int windowLength, int stepSize){
			FeatureCollection<PopulatedWindow<T1>> windows=convert(region.getWindows(windowLength, stepSize));
			formWindows(windows, reads);
		}
		
		private FeatureCollection<PopulatedWindow<T1>> convert(AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows) {
			FeatureCollection<PopulatedWindow<T1>> rtrn=new FeatureCollection<PopulatedWindow<T1>>();
			
			CloseableIterator<DerivedAnnotation<? extends Annotation>> iter=windows.sortedIterator();
			while(iter.hasNext()){
				DerivedAnnotation<? extends Annotation> window=iter.next();
				PopulatedWindow<T1> pw=new BlockedWindow<T1>(window);
				rtrn.add(pw);
			}
			return rtrn;
		}

		private void formWindows(FeatureCollection<PopulatedWindow<T1>> windows, CloseableIterator<T1> reads) {
			while(reads.hasNext()){
				T1 read=reads.next();
				CloseableIterator<PopulatedWindow<T1>> iter=windows.sortedIterator(read, false);
				while(iter.hasNext()){
					PopulatedWindow<T1> window=iter.next();
					if(window.getOrientation().equals(read.getOrientation()) || window.getOrientation().equals(Strand.BOTH)){
						window.addAnnotation(read);
					}
				}
			}
			this.windowIter=windows.sortedIterator();
		}

		@Override
		public boolean hasNext() {
			return windowIter.hasNext();
		}

		@Override
		public PopulatedWindow<T1> next() {
			return windowIter.next();
		}

		@Override
		public void close() {
			windowIter.close();
		}
		
		
		
	}
	

	/**
	 * This class requires that you have a sorted iterator of reads
	 * @author mguttman
	 *
	 * @param <T1>
	 */
	/**
	 * This class requires that you have a sorted iterator of reads
	 * @author mguttman
	 *
	 * @param <T1>
	 */
	
	
	public class WindowIterator<T1 extends Annotation> implements CloseableIterator<PopulatedWindow<T1>>{

		IntervalTree<PopulatedWindow<T1>> windows;
		CloseableIterator<T1> iter;
		Iterator<PopulatedWindow<T1>> fullyFormedWindows;
		int windowLength;
		boolean hasNext;
		boolean assumeForward;
		Annotation region;
		int stepSize;
		boolean includeEmpties;
		
		PopulatedWindow<T1> nextWin;
		int nextPos;
		
		public WindowIterator(CloseableIterator<T1> iter, int windowLength, boolean assumeForward){
			this.iter=iter;
			this.windowLength=windowLength;
			this.windows=new IntervalTree<PopulatedWindow<T1>>();
			this.hasNext=false;
			this.assumeForward=assumeForward;
			this.stepSize=1;
			this.includeEmpties=false;
			
		}
		
		public WindowIterator(CloseableIterator<T1> iter, int windowLength)
		{
			this(iter,windowLength,true);
			this.region = null;
		}
		
		public WindowIterator(CloseableIterator<T1> iter, int windowLength, int stepSize)
		{
			this(iter,windowLength,true);
			this.stepSize=stepSize;
			this.region = null;
		}

		public WindowIterator(CloseableIterator<T1> iter, int windowLength,Annotation region) {
			this(iter,windowLength,true);
			this.region = region;
			this.nextPos = region.getReferenceStartPosition();
		}
		
		public WindowIterator(CloseableIterator<T1> iter, int windowLength,Annotation region, int stepSize) {
			this(iter,windowLength,region);
			this.stepSize = stepSize;
		}

		public WindowIterator(CloseableIterator<T1> iter, int windowLength,Annotation region, int stepSize, boolean includeEmpties)
		{
			this(iter,windowLength,region);
			this.stepSize = stepSize;
			this.includeEmpties = includeEmpties;
		}
		
		@Override
		public boolean hasNext() {
			if(fullyFormedWindows!=null && fullyFormedWindows.hasNext()){return true;}
			boolean hasEmpties = (includeEmpties && region!=null && nextPos < region.getReferenceEndPosition());
			//else if(iter.hasNext()){updateWindows(); return hasNext();}
			//else if(!iter.hasNext() && !windows.isEmpty())
			//{
			//	updateRemainingWindows();
			//	return hasNext();
			//}
			
			updateWindows();
			return (hasEmpties || fullyFormedWindows!=null && fullyFormedWindows.hasNext());
			
			//return false;
		}

		private void updateRemainingWindows(){
			fullyFormedWindows = windows.valueIterator();
			windows = new IntervalTree<PopulatedWindow<T1>>();
		}
		
		@Override
		public PopulatedWindow<T1> next() {
			PopulatedWindow<T1> rtrn;
			
			if(nextWin == null)
			{
				if(fullyFormedWindows!=null && fullyFormedWindows.hasNext())
					nextWin = fullyFormedWindows.next();
				else
				{
					nextPos++;
					return new ContiguousWindow<T1>(region.getName(), nextPos, nextPos+windowLength, Strand.BOTH);

				}
			}
			
			if(!includeEmpties || nextWin.getReferenceStartPosition()==nextPos)
			{
				rtrn = nextWin;
				nextWin = null;
			}
			else
				rtrn = new ContiguousWindow<T1>(region.getName(), nextPos, nextPos+windowLength, Strand.BOTH);
			
			nextPos++;
			return rtrn;
		}

		private void updateWindows(){
			while((fullyFormedWindows==null || !fullyFormedWindows.hasNext()) && iter.hasNext())
			{	T1 read=iter.next();
				//all windows with an end position before the start of this window
				fullyFormedWindows=removeFullyFormedWindows(read).iterator();
				addReadToWindows(read);
			}
			if((fullyFormedWindows==null || !fullyFormedWindows.hasNext()) && !windows.isEmpty())
			{
				updateRemainingWindows();
			}
		}

		private void addReadToWindows(T1 read){
			//Create all windows which overlap read blocks
			Iterator<SingleInterval> interval=read.getBlocks();
			while(interval.hasNext()){
				SingleInterval block=interval.next();
				int start = 0;
				if(region!=null)
					start=Math.max(roundUp(region.getReferenceStartPosition()-windowLength+1), roundUp(block.getReferenceStartPosition()-windowLength));
				else
					start = Math.max(0,roundUp(block.getReferenceStartPosition()-windowLength));
				int end = block.getReferenceEndPosition();
				if(region!=null)
					end = Math.min(region.getReferenceEndPosition(), block.getReferenceEndPosition());
				for(int i=start; i<end; i+=stepSize){
					PopulatedWindow<T1> window=windows.remove(i, i+windowLength);
					if(window==null){
						//make a window
						window=new ContiguousWindow<T1>(read.getReferenceName(), i, i+windowLength, Strand.BOTH);
					}
					window.addAnnotation(read);
					windows.put(window.getReferenceStartPosition(), window.getReferenceEndPosition(), window);
				}
			}
		}
		
		private int roundUp(int winStart)
		{
			return (winStart + stepSize -1) / stepSize *stepSize;
		}
		
		private Collection<PopulatedWindow<T1>> removeFullyFormedWindows(T1 read) {
			Iterator<PopulatedWindow<T1>> iter;
			if(assumeForward)
				{iter=windows.getNodesBeforeInterval(read.getReferenceStartPosition(), read.getReferenceStartPosition());}
			else
				{iter=windows.getNodesAfterInterval(read.getReferenceEndPosition(), read.getReferenceEndPosition());}
			ArrayList<PopulatedWindow<T1>> rtrn=new ArrayList<PopulatedWindow<T1>>();
			while(iter.hasNext()){
				PopulatedWindow<T1> w=iter.next();
				if((!w.getReferenceName().equalsIgnoreCase(read.getReferenceName())) || (w.getReferenceEndPosition()<read.getReferenceStartPosition())){
					rtrn.add(w);
					windows.remove(w.getReferenceStartPosition(), w.getReferenceEndPosition());
				}
			}
			this.hasNext=!rtrn.isEmpty();
			Comparator comparator = new WindowComparator();
			Collections.sort(rtrn,comparator);
			return rtrn;
		}

		
		@Override
		public void remove() {
			this.iter.remove();
		}

		@Override
		public void close() {
			this.iter.close();
		}


	}

	public class WindowComparator implements Comparator<PopulatedWindow<T>>{
		@Override
		public int compare(PopulatedWindow<T> win1, PopulatedWindow<T> win2)
		{
			return win1.getReferenceStartPosition()-win2.getReferenceEndPosition();
		}
	}
	
	@Override
	public CoordinateSpace getFeatureCoordinateSpace(){
		//Iterate through all records
		CloseableIterator<T> iter=sortedIterator();
		Map<String, Integer> sizes=new TreeMap<String, Integer>();
		while(iter.hasNext()){
			T annotation=iter.next();
			int size=annotation.size();
			String name=annotation.getName();
			sizes.put(name, size);
		}
		iter.close();
		return new CoordinateSpace(sizes);
	}
	
	@Override
	public int getNumAnnotations(){
		if(this.numAnnotations != 0)
			return this.numAnnotations;

		CloseableIterator<T> iter = this.sortedIterator();
		int count = 0;
		while(iter.hasNext())
		{
			iter.next();
			count++;
		}
		iter.close();
		this.numAnnotations = count;
		return this.numAnnotations;
	}
	
	@Override
	public Collection<Annotation> toCollection(){
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		CloseableIterator<T> iter=this.sortedIterator();
		while(iter.hasNext()) {
			rtrn.add(iter.next());
		}
		
		return rtrn;
	}
}

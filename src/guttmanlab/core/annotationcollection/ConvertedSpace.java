package guttmanlab.core.annotationcollection;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.coordinatespace.CoordinateSpace;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import net.sf.samtools.util.CloseableIterator;

public class ConvertedSpace<T extends Annotation> extends AbstractAnnotationCollection<DerivedAnnotation<T>>{

	private AnnotationCollection<? extends Annotation> featureMapping;
	private AnnotationCollection<T> readMapping;
	private boolean fullyContained;
	
	public ConvertedSpace(AnnotationCollection<T> readMapping, AnnotationCollection<? extends Annotation> featureMapping, CoordinateSpace referenceSpaceForAnnotations, boolean fullyContained){
		this.readMapping=readMapping;
		this.featureMapping=featureMapping;
		this.fullyContained=fullyContained;
	}
	
	@Override
	public CloseableIterator<DerivedAnnotation<T>> sortedIterator() {
		return new CoordinateConverterIterator<T>(readMapping.sortedIterator(), featureMapping, fullyContained);
	}

	@Override
	public CloseableIterator<DerivedAnnotation<T>> sortedIterator(Annotation region, boolean fullyContained) {
		//Step 1: Convert region to old reference space
		Collection<DerivedAnnotation<Annotation>> newCoordinateMapping=convertCoordinates(region, fullyContained);
		CloseableIterator<T> iter = null;
		if(newCoordinateMapping.iterator().hasNext())
		{
			iter=readMapping.sortedIterator(newCoordinateMapping.iterator().next(), fullyContained);
			return new CoordinateConverterIterator<T>(iter, featureMapping, fullyContained);
		}
		else
			return new CoordinateConverterIterator<T>();
	}

	@Override
	public void writeToBAM(String fileName) {
		throw new IllegalArgumentException("not implemented");
	}

	@Override
	public void writeToBAM(String fileName, Annotation region, boolean fullyContained) {
		throw new IllegalArgumentException("not implemented");
	}

	public <X extends Annotation> Collection<DerivedAnnotation<X>> convertCoordinates(X annotation, boolean fullyContained){
		//Check if annotation is in Reference or Feature space
		
		if(featureMapping.getReferenceCoordinateSpace().contains(annotation)){
			//If in Reference space
			return convertFromReference(annotation, fullyContained);
		}
		
		else if(featureMapping.getFeatureCoordinateSpace().contains(annotation)){
			//If in Feature space
			return convertFromFeature(annotation, fullyContained);
		}
		
		else{
			throw new IllegalArgumentException(annotation.getReferenceName()+":"+annotation.getReferenceStartPosition()+"-"+annotation.getReferenceEndPosition()+" annotation is not mapped to either Reference or Feature space");
		}
		
	}
	
	public <X extends Annotation> Collection<DerivedAnnotation<X>> convertFromFeature(X featureAnnotation, boolean fullyContained){
		Collection<DerivedAnnotation<X>> rtrn=new ArrayList<DerivedAnnotation<X>>();
		CloseableIterator<? extends Annotation> iter=featureMapping.sortedIterator();
		while(iter.hasNext()){
			Annotation referenceAnnotation=iter.next();
			if(referenceAnnotation.getName().equals(featureAnnotation.getReferenceName())){
				//if the names are the same then it is the same feature
				//trim to relative start and end	
				Annotation a=referenceAnnotation.convertToReferenceSpace(featureAnnotation);
				DerivedAnnotation<X> dA=new DerivedAnnotation<X>(a, featureAnnotation);
				rtrn.add(dA);	
			}
		}
		iter.close();
		return rtrn;
	}
	
	public <X extends Annotation> Collection<DerivedAnnotation<X>> convertFromReference(X annotation, boolean fullyContained){
		Collection<DerivedAnnotation<X>> rtrn=new ArrayList<DerivedAnnotation<X>>();

		//Find features overlapping the annotation
		CloseableIterator<? extends Annotation> iter=featureMapping.sortedIterator(annotation, fullyContained);
	
			//Adjust the coordinates of the feature as needed in featureSpace (ie as distance from start and end)
		while(iter.hasNext()){
			Annotation feature=iter.next();
			Annotation intersect=feature.intersect(annotation); //TODO Consider whether to remove this, it may not be needed and is expensive per read
			//Is annotation fully contained in feature
			boolean isFullyContained=feature.fullyContained(annotation);
			if(intersect.size()>0){
				if(!fullyContained || isFullyContained){
					//Annotation interval=feature.convertToFeatureSpace(intersect);
					Annotation converted=annotation.convert(feature);
					if(converted!=null){
						DerivedAnnotation<X> dInterval=new DerivedAnnotation<X>(converted, annotation);
						rtrn.add(dInterval);
					}
				}
			}
		}
		iter.close();
		
		return rtrn;
	}
	
	
	@Override
	public CoordinateSpace getReferenceCoordinateSpace() {
		return featureMapping.getFeatureCoordinateSpace(); 
	}
	
	@Override
	public CoordinateSpace getFeatureCoordinateSpace() {
		return featureMapping.getReferenceCoordinateSpace(); 
	}

	@Override
	public CloseableIterator<PopulatedWindow<DerivedAnnotation<T>>> getPopulatedWindows(Annotation region, int windowLength)
	{
		//get read iterator overlapping the feature		
		CloseableIterator<DerivedAnnotation<T>> iter = sortedIterator(region,true); 
		//pass the read iterator to WindowIterator with boolean set by orientation of the feature
		CloseableIterator<PopulatedWindow<DerivedAnnotation<T>>> windows;
		if(region.getOrientation().equals(Strand.NEGATIVE))
			{windows = new WindowIterator<DerivedAnnotation<T>>(iter,windowLength,true);}
		else
			{windows = new WindowIterator<DerivedAnnotation<T>>(iter,windowLength,false);}
		//return the windowIterator
		return windows;
	}
	
	public class CoordinateConverterIterator<X extends Annotation> implements CloseableIterator<DerivedAnnotation<X>>{

		CloseableIterator<X> iter;
		Iterator<DerivedAnnotation<X>> next;
		boolean started;
		AnnotationCollection<? extends Annotation> mapping;
		boolean fullyContained;

		public CoordinateConverterIterator(CloseableIterator<X> iterator, AnnotationCollection<? extends Annotation> mapping, boolean fullyContained){
			this.iter=iterator;
			this.started=false;
			this.mapping=mapping;
			this.fullyContained=fullyContained;
		}

		public CoordinateConverterIterator() {
			//creates an empty iterator
			this.iter = null;
			this.started = false;
			this.mapping = null;
			this.fullyContained = false;
		}

		@Override
		public boolean hasNext() {
			if(iter==null)
			{
				//System.err.println("converted read iterator was empty.");
				return false;
			}
			if((!started || !next.hasNext()) && iter.hasNext())
			{
				started=true;
				findNext();
				return hasNext();
			}
			if(next!=null && next.hasNext()){
				return true;
			}

			return false;
		}

		private void findNext() {
			//We should iterate until we either run out of reads or find next
			boolean done=false;
			while(iter.hasNext() && !done){
				X annotation=iter.next();
				this.next=convertCoordinates(annotation, fullyContained).iterator();
				if(next.hasNext()){done=true;}
			}
		}

		@Override
		public DerivedAnnotation<X> next() {
			return next.next();
		}

		@Override
		public void remove() {
			this.iter.remove();
		}

		@Override
		public void close() {
			this.iter.close();

		}}
	
}

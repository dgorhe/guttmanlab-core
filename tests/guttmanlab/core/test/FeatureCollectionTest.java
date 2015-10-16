package guttmanlab.core.test;

/**
 * Created by cburghard on 10/15/15.
 */
import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.FeatureCollection;
import net.sf.samtools.util.CloseableIterator;
import org.junit.Before;
import org.junit.Test;

import java.io.IOException;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class FeatureCollectionTest {
    private FeatureCollection<Gene> smallCollection;
    private FeatureCollection<Gene> emptyCollection;
    private Map<String, FeatureCollection<Gene>> genesByReference;

    @Before
    public void setUp() throws Exception {
        BEDFileIO io = new BEDFileIO("tests/guttmanlab/core/data/refspace.txt");

        String refSeq = ( "tests/guttmanlab/core/data/refSeq.bed" );
        String small = ( "tests/guttmanlab/core/data/refSeq20.bed" );

        String referenceSizes = "tests/guttmanlab/core/data/sizes";

        genesByReference = BEDFileIO.loadFromFileByReferenceName(refSeq, referenceSizes);
        smallCollection = BEDFileIO.loadFromFileByReferenceName(refSeq, referenceSizes).get("chr1");
        emptyCollection = new FeatureCollection<Gene>(genesByReference.get("chr1").getReferenceCoordinateSpace());
    }

    @Test
    /*
    Tests that getCounts() returns the correct number of features for empty and nonempty FeatureCollections
     */
    public void testGetCounts() {
        FeatureCollection<Gene> chr1Features = genesByReference.get("chr1");
        assertEquals("chr1 has 1624 annotations.",1624,chr1Features.getCount());

        assertEquals("empty collection has no annotations",0,emptyCollection.getCount());
    }

    @Test
    /*
    Tests that sortedIterator() returns all objects in the collection, and that they are in sorted order
     */
    public void testSortedIterator() {
        FeatureCollection<Gene> chr1Features = genesByReference.get("chr1");
        CloseableIterator<Gene> si = chr1Features.sortedIterator();

        int count = 0;
        boolean nullFound = false;
        boolean sorted = true;
        int prevStart = 0;

        while(si.hasNext())
        {
            Gene g = si.next();
            if(g==null)
                nullFound = true;
            if(g.getReferenceStartPosition() < prevStart)
                sorted = false;
            prevStart = g.getReferenceStartPosition();
            count++;
        }

        assertEquals("chr1 has 1624 annotations.",1624,count);
        assertEquals("no annotations were null",false,nullFound);
        assertEquals("annotations were sorted by start position",true,sorted);
    }

    @Test
    /*
    Tests that sortedIterator(Annotation a) returns only features that overlap a for both fullyContained and
    !fullyContained cases.  Checks that returned annotations are in sorted order.
     */
    public void testSortedIteratorOverRegion() {
        FeatureCollection<Gene> chr1Features = genesByReference.get("chr1");
        Annotation a = new SingleInterval("chr1",63492580,63799400, Annotation.Strand.BOTH);
        CloseableIterator<Gene> siFullyContained = chr1Features.sortedIterator(a,true);
        CloseableIterator<Gene> siNotFullyContained = chr1Features.sortedIterator(a,false);

        int count = 0;
        while(siFullyContained.hasNext())
        {
            Gene g = siFullyContained.next();
            System.out.println(g.toBED());
            count++;
        }

        int count2 = 0;
        while(siNotFullyContained.hasNext())
        {
            Gene g = siNotFullyContained.next();
            count2++;
        }
        assertEquals("region has 2 annotations fully contained.",2,count);
        assertEquals("region has 5 annotations overlapping.",5,count2);
    }

    @Test
    /*
    Tests contains() when annotation is/is not in the collection
     */
    public void testContains() {
        FeatureCollection<Gene> chr1Features = genesByReference.get("chr1");

        Annotation inCollection = null;
        CloseableIterator<Gene> si = chr1Features.sortedIterator();
        if(si.hasNext())
            inCollection = si.next();
        SingleInterval block = new SingleInterval("chr2",100,200, Annotation.Strand.BOTH);
        BlockedAnnotation notInCollection = new BlockedAnnotation(block);

        assertEquals("chr1 gene should be found in collection",true,chr1Features.contains(inCollection));
        assertEquals("chr2 gene should not be found in collection",false,chr1Features.contains(notInCollection));
    }

    @Test
    /*
    Tests containsAll for a collection of annotations
     */
    public void testContainsAll() {
        FeatureCollection<Gene> chr1Features = genesByReference.get("chr1");

        assertEquals("all refSeq genes found in collection",true,chr1Features.containsAll(smallCollection));
    }

    @Test
    /*
    Tests removal of an object from the collection. Verifies contains() then returns false and getCounts() reduced by 1
     */
    public void testRemove() {
        FeatureCollection<Gene> chr1Features = genesByReference.get("chr1");
        Annotation removed = null;
        CloseableIterator<Gene> si = chr1Features.sortedIterator();
        if(si.hasNext())
            removed = si.next();
        int oldCount = chr1Features.getCount();
        chr1Features.remove(removed);

        assertEquals("contained() returns false for removed item",false,chr1Features.contains(removed));
        assertEquals("removed decreases count by 1", oldCount - 1, chr1Features.getCount());
    }
}
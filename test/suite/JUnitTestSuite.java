package suite;

import guttmanlab.core.alignment.TestLocalAlignmentUtils;
import guttmanlab.core.alignment.TestSmithWatermanAlignment;
import guttmanlab.core.annotation.TestAbstractAnnotation;
import guttmanlab.core.annotation.TestAnnotation;
import guttmanlab.core.annotation.TestAnnotationFileRecord;
import guttmanlab.core.annotation.TestAnnotationHelper;
import guttmanlab.core.annotation.TestBEDFileRecord;
import guttmanlab.core.annotation.TestBlockedAnnotation;
import guttmanlab.core.annotation.TestContiguousWindow;
import guttmanlab.core.annotation.TestDerivedAnnotation;
import guttmanlab.core.annotation.TestGene;
import guttmanlab.core.annotation.TestMappedFragment;
import guttmanlab.core.annotation.TestPairedMappedFragment;
import guttmanlab.core.annotation.TestPopulatedWindow;
import guttmanlab.core.annotation.TestSAMFragment;
import guttmanlab.core.annotation.TestScore;
import guttmanlab.core.annotation.TestSingleInterval;
import guttmanlab.core.annotation.io.TestAnnotationFileIO;
import guttmanlab.core.annotation.io.TestBEDFileIO;
import guttmanlab.core.annotation.predicate.TestContainedByFilter;
import guttmanlab.core.annotation.predicate.TestIndelFilter;
import guttmanlab.core.annotation.predicate.TestInsertSizeFilter;
import guttmanlab.core.annotation.predicate.TestMappedReadFilter;
import guttmanlab.core.annotation.predicate.TestMaximumLengthFilter;
import guttmanlab.core.annotation.predicate.TestMinimumLengthFilter;
import guttmanlab.core.annotation.predicate.TestOverlapsFilter;
import guttmanlab.core.annotation.predicate.TestPairedFilterWrapper;
import guttmanlab.core.annotation.predicate.TestReadClippedFilter;
import guttmanlab.core.annotation.predicate.TestReadFlag;
import guttmanlab.core.annotation.predicate.TestSAMFragmentNumHitsFilter;
import guttmanlab.core.annotation.predicate.TestSecondReadFilter;
import guttmanlab.core.annotation.predicate.TestStrandFilter;
import guttmanlab.core.annotationcollection.TestAbstractAnnotationCollection;
import guttmanlab.core.annotationcollection.TestBAMFragmentCollectionFactory;
import guttmanlab.core.annotationcollection.TestBAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.TestBAMSingleReadCollection;
import guttmanlab.core.annotationcollection.TestConvertedSpace;
import guttmanlab.core.annotationcollection.TestFeatureCollection;
import guttmanlab.core.annotationcollection.TestFilteredIterator;
import guttmanlab.core.annotationcollection.TestNameIgnoringIterator;
import guttmanlab.core.coordinatespace.TestCoordinateSpace;
import guttmanlab.core.coordinatespace.TestGenomeSize;
import guttmanlab.core.datastructures.TestInterval;
import guttmanlab.core.datastructures.TestIntervalTree;
import guttmanlab.core.datastructures.TestPair;
import guttmanlab.core.math.TestMaximumContiguousSubsequence;
import guttmanlab.core.math.TestScanStat;
import guttmanlab.core.math.TestStatistics;
import guttmanlab.core.pipeline.TestConfigFile;
import guttmanlab.core.pipeline.TestConfigFileOption;
import guttmanlab.core.pipeline.TestConfigFileOptionValue;
import guttmanlab.core.pipeline.TestConfigFileSection;
import guttmanlab.core.pipeline.TestJob;
import guttmanlab.core.pipeline.TestJobUtils;
import guttmanlab.core.pipeline.TestLSFJob;
import guttmanlab.core.pipeline.TestOGSJob;
import guttmanlab.core.pipeline.TestOGSUtils;
import guttmanlab.core.pipeline.TestScheduler;
import guttmanlab.core.pipeline.util.TestBamUtils;
import guttmanlab.core.pipeline.util.TestFastaUtils;
import guttmanlab.core.pipeline.util.TestFastqParser;
import guttmanlab.core.pipeline.util.TestFastqSequence;
import guttmanlab.core.pipeline.util.TestFastqUtils;
import guttmanlab.core.sequence.TestFastaFileIO;
import guttmanlab.core.sequence.TestSequence;
import guttmanlab.core.sequence.TestFastaFileIOImpl;
import guttmanlab.core.serialize.TestAbstractAvroIndex;
import guttmanlab.core.serialize.TestAvroIndex;
import guttmanlab.core.serialize.TestAvroStringIndex;
import guttmanlab.core.serialize.TestBuildAvroIndex;
import guttmanlab.core.serialize.sam.TestAvroSamRecord;
import guttmanlab.core.serialize.sam.TestAvroSamStringIndex;
import guttmanlab.core.serialize.sam.TestSerializeBam;
import guttmanlab.core.util.TestCommandLineParser;
import guttmanlab.core.util.TestCountLogger;
import guttmanlab.core.util.TestMismatchGenerator;
import guttmanlab.core.util.TestSAMFlagDecoder;
import guttmanlab.core.util.TestStringParser;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

@RunWith(Suite.class)
@Suite.SuiteClasses({
	// alignment
	TestLocalAlignmentUtils.class,
	TestSmithWatermanAlignment.class,
	// annotation
	TestAbstractAnnotation.class,
	TestAnnotation.class,
	TestAnnotationFileRecord.class,
	TestAnnotationHelper.class,
	TestBEDFileRecord.class,
	TestBlockedAnnotation.class,
	TestContiguousWindow.class,
	TestDerivedAnnotation.class,
	TestGene.class,
	TestMappedFragment.class,
	TestPairedMappedFragment.class,
	TestPopulatedWindow.class,
	TestSAMFragment.class,
	TestScore.class,
	TestSingleInterval.class,
	// annotation.io
	TestAnnotationFileIO.class,
	TestBEDFileIO.class,
	// annotation.predicate
	TestContainedByFilter.class,
	TestIndelFilter.class,
	TestInsertSizeFilter.class,
	TestMappedReadFilter.class,
	TestMaximumLengthFilter.class,
	TestMinimumLengthFilter.class,
	TestOverlapsFilter.class,
	TestPairedFilterWrapper.class,
	TestReadClippedFilter.class,
	TestReadFlag.class,
	TestSAMFragmentNumHitsFilter.class,
	TestSecondReadFilter.class,
	TestStrandFilter.class,
	// annotationcollection
	TestAbstractAnnotationCollection.class,
	TestBAMFragmentCollectionFactory.class,
	TestBAMPairedFragmentCollection.class,
	TestBAMSingleReadCollection.class,
	TestConvertedSpace.class,
	TestFeatureCollection.class,
	TestFilteredIterator.class,
	TestNameIgnoringIterator.class,
	// coordinatespace
	TestCoordinateSpace.class,
	TestGenomeSize.class,
	// datastructures
	TestInterval.class,
	TestIntervalTree.class,
	TestPair.class,
	// math
	TestMaximumContiguousSubsequence.class,
	TestScanStat.class,
	TestStatistics.class,
	// pipeline
	TestConfigFile.class,
	TestConfigFileOption.class,
	TestConfigFileOptionValue.class,
	TestConfigFileSection.class,
	TestJob.class,
	TestJobUtils.class,
	TestLSFJob.class,
	TestOGSJob.class,
	TestOGSUtils.class,
	TestScheduler.class,
	// pipeline.util
	TestBamUtils.class,
	TestFastaUtils.class,
	TestFastqParser.class,
	TestFastqSequence.class,
	TestFastqUtils.class,
	// sequence
	TestFastaFileIO.class,
	TestFastaFileIOImpl.class,
	TestSequence.class,
	// serialize
	TestAbstractAvroIndex.class,
	TestAvroIndex.class,
	TestAvroStringIndex.class,
	TestBuildAvroIndex.class,
	// serialize.sam
	TestAvroSamRecord.class,
	TestAvroSamStringIndex.class,
	TestSerializeBam.class,
	// util
	TestCommandLineParser.class,
	TestCountLogger.class,
	TestMismatchGenerator.class,
	TestSAMFlagDecoder.class,
	TestStringParser.class
})

public class JUnitTestSuite {}

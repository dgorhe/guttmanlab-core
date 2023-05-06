package guttmanlab.core.math;

public class MaximumContiguousSubsequence {

	
	//static private int seqStart = 0;
    //static private int seqEnd = -1;

	
	/**
     * Linear-time maximum contiguous subsequence sum algorithm.
     * seqStart and seqEnd represent the actual best sequence.
     */
    public static int[] maxSubSum3(int[] a)
    {
        int maxSum = 0;
        int thisSum = 0;
        
        int seqStart=0;
        int seqEnd=-1;

        for( int i = 0, j = 0; j < a.length; j++ )
        {
            thisSum += a[ j ];

            if( thisSum > maxSum )
            {
                maxSum = thisSum;
                seqStart = i;
                seqEnd   = j;
            }
            else if( thisSum < 0 )
            {
                i = j + 1;
                thisSum = 0;
            }
        }

        int[] rtrn={maxSum, seqStart, seqEnd};
        
        return rtrn;
    }

    
    public static int maxFromEnd(double[] vals) {
		//Fix end position and move start
		double[] sum=new double[vals.length];
		
		double sumSoFar=0;
		for(int i=vals.length-1; i>=0; i--) {
			sum[i]=sumSoFar+vals[i];
			//System.err.println(i+" "+sum[i]);
			sumSoFar=sum[i];
		}
		
		double maxSoFar=-1;
		int pos=-1;
		for(int i=0; i<sum.length; i++) {
			if(sum[i]>maxSoFar) {
				maxSoFar=sum[i];
				pos=i;
			}
		}
		
		return pos;
	}
    
    
    public static int maxFromStart(double[] vals) {
		//Fix start position and move end
		double[] sum=new double[vals.length];
		
		double sumSoFar=0;
		for(int i=0; i<vals.length; i++) {
			sum[i]=sumSoFar+vals[i];
			//System.err.println(i+" "+sum[i]);
			sumSoFar=sum[i];
		}
		
		double maxSoFar=-1;
		int pos=-1;
		for(int i=0; i<sum.length; i++) {
			if(sum[i]>maxSoFar) {
				maxSoFar=sum[i];
				pos=i;
			}
		}
		
		return pos;
	}
    
    
    public static int[] maxSubSum3(double[] a)
    {
        double maxSum = 0;
        double thisSum = 0;
        
        int seqStart=0;
        int seqEnd=-1;

        for( int i = 0, j = 0; j < a.length; j++ )
        {
            thisSum += a[ j ];

            if( thisSum > maxSum )
            {
                maxSum = thisSum;
                seqStart = i;
                seqEnd   = j;
            }
            else if( thisSum < 0 )
            {
                i = j + 1;
                thisSum = 0;
            }
        }

        int[] rtrn={seqStart, seqEnd};
        
        return rtrn;
    }
    	 
    
    public static int contiguousStartSubSequenceOverMin(double [] a, double min) {
    	int i = 0;
    	while (i < a.length && a[i++] < min) { 
    		;
    	}
    	
    	return i;
    }
    
    public static int contiguousEndSubSequenceOverMin(double [] a, double min) {
    	int i = a.length;
    	while (--i >=0 && a[i] < min) { 
    		;
    	}
    	
    	return i;
    }
    
    
    public static void main( String [ ] args )
    {
        double a[ ] = {-1,+1,-1,-1,1,1, 1,-1,1,1,-1,-1};
       int pos = maxFromStart( a );
        System.out.println(pos);
    }


	
}

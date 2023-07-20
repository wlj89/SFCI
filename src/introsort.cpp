#include "all.h"

const int threshold = 16; 

int lg(int n) 
{
    int k;
    for(k=0;n>1;n>>=1) ++k; 
    return k;
}

void copy_bckwd(dynVec<DET_TYPE,FLOAT_TYPE>& array, int begin, int i_pos)
{
    
    DET_TYPE det_prev = array.basis[begin]; 
    DET_TYPE det_now; 

    FLOAT_TYPE coeff_prev = array.amp[begin];        
    FLOAT_TYPE coeff_now; 

    for(int i = begin+1; i< i_pos+1; i++) 
    {
        det_now =  array.basis[i];
        coeff_now = array.amp[i];

        array.basis[i] = det_prev;  
        array.amp[i] = coeff_prev; 

        det_prev = det_now;
        coeff_prev = coeff_now; 

    }   
    //copy_backward(array.basis.begin()+begin, array.basis.begin()+i, array.basis.begin()+i+1);

}

/*
    introsort for basis 
*/

int median3(dynVec<DET_TYPE,FLOAT_TYPE>& array,int first,int median,int end)
{
    if(array.basis[first]<array.basis[median])
    {
        if(array.basis[median]<array.basis[end])
            return median;
        else if(array.basis[first]<array.basis[end])
            return end; 
        else
            return first;
    }
    else if(array.basis[first]<array.basis[end])
        return first;
    else if(array.basis[median]<array.basis[end])
        return end;
    else 
        return median;

}

void push_heap_1(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                        int begin, 
                                        int holeIndex, 
                                        int topIndex,
                                   DET_TYPE value_basis, 
                                 FLOAT_TYPE value_amp)
{
    int parent = (holeIndex-1)/2;
    
    while((holeIndex > topIndex) && array.basis[begin+parent] < value_basis)
    {
        array.pasteFromIdx(begin+holeIndex, begin+parent);
        holeIndex = parent; 
        parent = (holeIndex-1)/2;

    }   

    array.assignValue(begin+holeIndex, value_basis, value_amp);
}

void adjust_heap(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                          int begin, 
                                          int holeIndex, 
                                          int len, 
                                     DET_TYPE value_basis,
                                   FLOAT_TYPE value_amp)
{
    int topIndex = holeIndex;
    int secondChild = 2 * holeIndex + 2; 

    while(secondChild < len)
    {
        if (array.basis[begin+secondChild] < array.basis[begin+secondChild-1] )
            secondChild --;

        array.pasteFromIdx(begin+holeIndex, begin+secondChild);
        holeIndex = secondChild;
        secondChild = 2 * (secondChild + 1); 
    }

    if (secondChild == len )    
    {
        array.pasteFromIdx(begin+holeIndex, begin+secondChild-1);
        holeIndex = secondChild - 1; 
    }

    push_heap_1(array,begin, holeIndex, topIndex, value_basis, value_amp); 
}   

void make_heap(dynVec<DET_TYPE,FLOAT_TYPE>& array, int begin, int end)
{
    

    if (end - begin < 2) 
        return;
    int len = end - begin; 
    int parent = (len-2)/2; 

    while (1) 
    {
        adjust_heap(array, 
                    begin, 
                    parent, 
                    len, 
                    array.basis[begin+parent], 
                    array.amp[begin+parent]);

        if (parent == 0)
            return;
        //cout << "parent=" << parent << endl; 
        parent--;
    }
}

void pop_heap(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                       int begin, 
                                       int end, 
                                       int result, 
                                  DET_TYPE value_basis, 
                                FLOAT_TYPE value_amp)
{   
    array.pasteFromIdx(result, begin);
    adjust_heap(array,begin, 0, end-begin, value_basis, value_amp); 
}

void sort_heap(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                       int begin, 
                                       int end)
{
    while(end-begin > 1)
    {
        pop_heap(array, begin, end-1, end-1, array.basis[end-1],array.amp[end-1]); 
        end-=1; 
    }
}

int unguarded_partition(dynVec<DET_TYPE,FLOAT_TYPE>& array, int begin, int end, int pivotPos)
{
    auto pivot = array.basis[pivotPos]; 

    while(1)
    {
        while (array.basis[begin] < pivot) 
            begin++;
        
        end--;

        while (pivot < array.basis[end]) 
            end--; 

        if(!(begin < end)) 
            return begin; 

        array.swapIdx(begin, end);

        begin++; 
    }
}

void partial_sort(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                           int begin, 
                                           int middle, 
                                           int end)
{   
    /*
        the heap sort
        first form a max heap  
    */

    make_heap(array, begin, middle); 

    /* this block does not appear useful*/
    for(int i = middle; i <end; i++)
    {
        if(array.basis[i] < array.basis[begin])
            pop_heap(array, begin, middle, i, array.basis[i],array.amp[i]);
    }

    // before this should be a max heap 
    sort_heap(array, begin, middle); 
}

void introsort_loop(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                             int begin,
                                             int end,
                                             int depthLimit)
{
    while(end - begin > threshold) 
    {   
        // threshold = 16
        if(depthLimit==0)    
        {
            // heap sort starts
            partial_sort(array,begin,end,end);  
            
            return ;
        }   
        
        --depthLimit; 
        
        //int cut = unguarded_partition(array,begin,end,median3(array,begin,begin+(end-begin)/2,end)); 
        
        /*
            a quicksort style partition 
        */
        int cut = unguarded_partition(  array,
                                        begin,
                                        end,
                                        median3(array,begin,begin+(end-begin)/2,end-1)); 
        
        // take care of the segment on the right 
        introsort_loop(array,cut,end,depthLimit);
        
        end=cut;    
    }
}   


void unguarded_linear_insertion(dynVec<DET_TYPE,FLOAT_TYPE>& array,
                                                         int end,
                                                         DET_TYPE value_basis,
                                                       FLOAT_TYPE value_amp)
{
    int next = end; 
    next--; 

    while(value_basis < array.basis[next])
    {
        /*
            "next" should be bounded, but it goes to minus. 
            This means the global minimum is not showing up in the first 16 elelemnts. 

        */
        array.pasteFromIdx(end,next);
        end = next; 
        next--;
        //cout << "next = " << next<<endl;
    }
    //cout << "****************************\n";

    array.assignValue(end,value_basis, value_amp);
}

void insertion_sort(dynVec<DET_TYPE,FLOAT_TYPE>& array,
                                             int begin,
                                             int end)
{   
    if(begin == end)
        return; 

    for (int i = begin + 1; i != end; i++)
    {
        auto value_basis = array.basis[i];
        auto value_amp = array.amp[i];

        if (value_basis < array.basis[begin])
        {
            //copy_backward(array.basis.begin()+begin, array.basis.begin()+i, array.basis.begin()+i+1);
            //copy_backward(array.amp.begin()+begin, array.amp.begin()+i, array.amp.begin()+i+1);
            copy_bckwd(array,begin,i); 
            array.assignValue(begin, value_basis, value_amp);
        }
        else 
        {   
            unguarded_linear_insertion(array, i, value_basis,value_amp); 
        }
    }
}

void final_insertion_sort(dynVec<DET_TYPE,FLOAT_TYPE>& array,int begin, int end)
{
    if(end - begin > threshold)
    {
        insertion_sort(array, begin, begin + threshold);
        
        // unguarded insertion sort 
        for (int i = begin + threshold ; i< end ; i++)
            unguarded_linear_insertion(array, i, array.basis[i], array.amp[i]);     
        //unguarded_insertion_sort(array, begin + threshold,end); 
    }
    else
    {
        insertion_sort(array, begin, end);
    }

}

void introsort(dynVec<DET_TYPE,FLOAT_TYPE>& array,int len)
{
    if(len>1)
    {
        //swap(array.amp[0],array.amp[len-1]);
		//swap(array.basis[0],array.basis[len-1]);
		
        //introsort_loop(array,0,len,lg(len)*2);  
        introsort_loop(array,0,len,lg(len)*2);     

        //cout << "beep\n";
		final_insertion_sort(array,0,len); 
    }
} 

/*
    introsort for amp 
*/

int median3_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array,int first,int median,int end)
{
    if(abs(array.amp[first])<abs(array.amp[median]))
    {
        if(abs(array.amp[median])<abs(array.amp[end]))
            return median;
        else if(abs(array.amp[first])<abs(array.amp[end]))
            return end; 
        else
            return first;
    }
    else if(abs(array.amp[first])<abs(array.amp[end]))
        return first;
    else if(abs(array.amp[median])<abs(array.amp[end]))
        return end;
    else 
        return median;
}

void push_heap_1_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                        int begin, 
                                        int holeIndex, 
                                        int topIndex,
                                   DET_TYPE value_basis, 
                                 FLOAT_TYPE value_amp)

{ 
    int parent = (holeIndex-1)/2;
    
    while((holeIndex > topIndex) && abs(array.amp[begin+parent]) < abs(value_amp))
    {
        array.pasteFromIdx(begin+holeIndex, begin+parent);
        holeIndex = parent; 
        parent = (holeIndex-1)/2;
    }   

    array.assignValue(begin+holeIndex, value_basis, value_amp);

}

void adjust_heap_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                          int begin, 
                                          int holeIndex, 
                                          int len, 
                                     DET_TYPE value_basis,
                                   FLOAT_TYPE value_amp)
{
    int topIndex = holeIndex;
    int secondChild = 2 * holeIndex + 2; 
    while(secondChild < len)
    {
        if (abs(array.amp[begin+secondChild]) < abs(array.amp[begin+secondChild-1]) )
            secondChild --;

        array.pasteFromIdx(begin+holeIndex, begin+secondChild);
        holeIndex = secondChild;
        secondChild = 2 * (secondChild + 1); 
    }

    if (secondChild == len )    
    {
        array.pasteFromIdx(begin+holeIndex, begin+secondChild-1);
        holeIndex = secondChild - 1; 
    }

    push_heap_1_amp(array,begin, holeIndex, topIndex, value_basis, value_amp); 
}   

void make_heap_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array, int begin, int end)
{
    if (begin - end < 2) 
        return;
    int len = end - begin; 
    int parent = (len-2)/2; 

    while (1) 
    {
        adjust_heap_amp(array, begin, parent, len, array.basis[begin+parent], array.amp[begin+parent]);

        if (parent == 0)
            return;
        parent--;
    }
}

void pop_heap_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                       int begin, 
                                       int end, 
                                       int result, 
                                  DET_TYPE value_basis, 
                                FLOAT_TYPE value_amp)
{   
    array.pasteFromIdx(result, begin);
    adjust_heap_amp(array,begin, 0, end-begin, value_basis, value_amp); 
}

void sort_heap_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                       int begin, 
                                       int end)
{
    while(end-begin > 1)
    {
        pop_heap_amp(array, begin, end-1, end-1, array.basis[end-1],array.amp[end-1]); 
        end--; 
    }
}

int unguarded_partition_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array, int begin, int end, int pivotPos)
{
    /*

        this is a trick part. 
    */

    auto pivot = abs(array.amp[pivotPos]); 

    while(1)
    {
        while (abs(array.amp[begin]) < pivot) 
            begin++;
        
        end--;

        while (pivot < abs(array.amp[end])) 
            end--; 

        if(!(begin < end)) 
            return begin; 

        array.swapIdx(begin, end);

        begin++; 
    }
}

void partial_sort_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array, 
                                           int begin, 
                                           int middle, 
                                           int end)
{   
    make_heap_amp(array, begin, middle); 

    for(int i = middle; i <end; i++ )
        if(abs(array.amp[i]) < abs(array.amp[begin]))
            pop_heap_amp(array, begin, middle, i, array.basis[i],array.amp[i]);

    sort_heap_amp(array, begin, middle); 
}

void introsort_loop_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array, int begin,int end,int depthLimit)
{
    while(end - begin > threshold) 
    {   
        // threshold = 16
        if(depthLimit==0)    
        {
            //heapSortBasis(array,begin,end);
            partial_sort_amp(array,begin,end,end);  
            return ;
        }   
        
        --depthLimit; 
        
        //int cut=partitionBasis(array,begin,end,median3Basis(array,begin,begin+(end-begin)/2,end)); 
        /*
            pivot position is wrong: should be end-1! 
        */
        //int cut = unguarded_partition_amp(array,begin,end,median3_amp(array,begin,begin+(end-begin)/2,end)); 
        int cut = unguarded_partition_amp(array,begin,end,median3_amp(array,begin,begin+(end-begin)/2,end-1)); 
        introsort_loop_amp(array,cut,end,depthLimit);
        end=cut;    
    }
}   

void unguarded_linear_insertion_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array,
                                                         int end,
                                                         DET_TYPE value_basis,
                                                       FLOAT_TYPE value_amp)
{
    int next = end; 
    next--; 
	FLOAT_TYPE v_abs = abs(value_amp);// < abs( array.amp[next]))
    
	while(v_abs < abs( array.amp[next]))
    {
        array.pasteFromIdx(end,next);
        end = next; 
        next--;
    }

    array.assignValue(end,value_basis, value_amp);
}

void insertion_sort_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array,
                                             int begin,
                                             int end)
{   
    if(begin == end)
        return; 

    for (int i = begin + 1; i != end; i++)
    {
        auto value_basis = array.basis[i];
        auto value_amp = array.amp[i];

        if (abs(value_amp) < abs(array.amp[begin]))
        {
            //copy_backward(array.basis.begin()+begin, array.basis.begin()+i, array.basis.begin()+i+1);
            //copy_backward(array.amp.begin()+begin, array.amp.begin()+i, array.amp.begin()+i+1);
            copy_bckwd(array,begin,i);
            
            array.assignValue(begin, value_basis, value_amp);
        }
        else 
        {
            unguarded_linear_insertion_amp(array, i, value_basis,value_amp); 
        }
    }
}   

void final_insertion_sort_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array,int begin, int end)
{
    if(end - begin > threshold)
    {
        insertion_sort_amp(array, begin, begin + threshold);
        
        // unguarded insertion sort 
		//cout << "beep" <<endl;
	    //for (int i=0; i< threshold; i++)
		//	cout << array.amp[i] <<endl; 	
        for (int i = begin + threshold ; i< end ; i++)
		{
            //if (i%1000==0)
			//	cout << i <<endl; 
			unguarded_linear_insertion_amp(array, i, array.basis[i], array.amp[i]);     
		}
		//unguarded_insertion_sort(array, begin + threshold,end); 
    }   
    else
    {
        insertion_sort_amp(array, begin, end);
    }   
}

void introsort_amp(dynVec<DET_TYPE,FLOAT_TYPE>& array,int len)
{
    if(len>1) 
    {

		introsort_loop_amp(array,0,len,lg(len)*2);    
        //cout << "beep_amp\n";
		final_insertion_sort_amp(array,0,len); 
    }
} 

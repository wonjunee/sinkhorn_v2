#include <stdlib.h>
#include <string.h>
#ifndef BRANCHING_NUMBER
#define BRANCHING_NUMBER 8
#endif
#ifndef MAJ_SIMPLE_IMPLICIT_HEAP
#define MAJ_SIMPLE_IMPLICIT_HEAP

class i_heap_node{
public:
    k_type key;
    int originalIndex;

    i_heap_node(){
        key = 0;
        originalIndex = 0;
    }
};

class i_heap{
public:
    i_heap_node* root;
    int* locations;
    int count;

    i_heap(){
        root = NULL;
        locations = NULL;
    }

    i_heap(int ncount, int locationCount){
        count=0;
        root = new i_heap_node[ncount];
        locations = new int[locationCount];
        memset(locations,-1,locationCount*sizeof(int));
    }

    ~i_heap(){
        delete[] root;
        delete[] locations;
    }

    void i_heap_clear_heap(){
        count=0;
    }

    void i_heap_swap_nodes(int nodeIndex1, int nodeIndex2){
        int i1=root[nodeIndex1].originalIndex;
        int i2=root[nodeIndex2].originalIndex;
        k_type key1=root[nodeIndex1].key;
        k_type key2=root[nodeIndex2].key;
        root[nodeIndex1].originalIndex=i2;
        root[nodeIndex2].originalIndex=i1;
        root[nodeIndex1].key=key2;
        root[nodeIndex2].key=key1;
        locations[i1]=nodeIndex2;
        locations[i2]=nodeIndex1;
        
    }


    void i_heap_bubble_up(int nodeIndex){
        
        while (nodeIndex>0) {
            k_type myKey=root[nodeIndex].key;
            int parentIndex=nodeIndex/BRANCHING_NUMBER-((nodeIndex%BRANCHING_NUMBER)==0);
            k_type parentKey=root[parentIndex].key;
            if(myKey<parentKey){
                i_heap_swap_nodes(nodeIndex,parentIndex);
                nodeIndex=parentIndex;
            }else{
                break;
            }
        }
    }

    void i_heap_push_down(int nodeIndex){
        
        int i = 0;
        int childIndex=nodeIndex*BRANCHING_NUMBER+1;
        while(childIndex<count){
            int minIndex=childIndex;
            k_type myKey=root[nodeIndex].key;
            k_type min=myKey;
            for(i=0;i<BRANCHING_NUMBER;i++){
                if(childIndex+i<count){
                    k_type childKey=root[childIndex+i].key;
                    if(childKey<min){
                        min=childKey;
                        minIndex=childIndex+i;
                    }
                }
            }
            if(min<myKey){
                i_heap_swap_nodes(nodeIndex,minIndex);
                nodeIndex=minIndex;
                childIndex=nodeIndex*BRANCHING_NUMBER+1;
            }else{
                break;
            }
            
        }
    }    

    void i_heap_add_node_to_bottom_with_location(int originalIndex, k_type key, int location){
        
        root[count].originalIndex=originalIndex;
        root[count].key=key;
        locations[location]=count;
        count++;
    }


    void i_heap_delete_min(){
        int index=root[0].originalIndex;
        locations[index]=-1;
        i_heap_swap_nodes(count-1,0);
        count--;
        i_heap_push_down(0);
    }


    void i_heap_decrease_key(int nodeIndex, k_type newKey){  
        root[nodeIndex].key=newKey;
        i_heap_bubble_up(nodeIndex);
    }



    int i_heap_empty(){
        return count==0;
    }


};



class Fast_Marching{
public:
    int n1;
    int n2;
    int xgrid[4];
    int ygrid[4];

    double eps;

    Fast_Marching(int n1, int n2){
        this->n1=n1;
        this->n2=n2;

        xgrid[0] = 1; xgrid[1] = 0; xgrid[2] = -1; xgrid[3] = 0;
        ygrid[0] = 0; ygrid[1] = 1; ygrid[2] =  0; ygrid[3] =-1;

        eps = 1.0/sqrt(n1);
    }

    ~Fast_Marching(){
        
    }


    /* typical 4-stencil upwind eikonal update */
    double eikonal_update(double *u, int x, int y){
        double result=0;
        int xleft=fmax(0,x-1);
        int xright=fmin(n1-1,x+1);
        int ydown=fmax(0,y-1);
        int yup=fmin(n2-1,y+1);
        
        int left=y*n1+xleft;
        int right=y*n1+xright;
        int down=ydown*n1+x;
        int up=yup*n1+x;
        
       
        double finverse=1;
        
        double dx=1/(n1*1.0);
        double dy=1/(n2*1.0);
        
        double alpha=pow(dx/dy,2);
        
        double minx=fmin(u[right],u[left]);
       
        double miny=fmin(u[up],u[down]);
       
        double a=1+alpha;
        double b=-2*(minx+alpha*miny);
        double c=minx*minx+alpha*miny*miny-pow(dx*finverse,2);
        double discrim=b*b-4*a*c;
        
        if(discrim>0){
            result=(sqrt(b*b-4*a*c)-b)/(2*a);
        }else if(miny<minx){
            result=miny+dy*finverse;
        }else{
            result=minx+dx*finverse;
        }
        
        return result;
    }


    static void i_heap_insert_node_with_location(i_heap& heap, int originalIndex, k_type key, int location){
        
        int count=heap.count;
        heap.root[count].originalIndex=originalIndex;
        heap.root[count].key=key;
        heap.locations[location]=count;
        heap.count++;
        heap.i_heap_bubble_up( count);
    }


    double single_fast_marching_iteration(i_heap& heap, double *u){
        int currentIndex=heap.root[0].originalIndex;
        double min=heap.root[0].key;
        int cx=currentIndex%n1;
        int cy=currentIndex/n1;
        for(int k=0;k<4;k++){
            int x=fmin(n1-1,fmax(0,cx+xgrid[k]));
            int y=fmin(n2-1, fmax(0,cy+ygrid[k]));
            int index=y*n1+x;
            /* check if the point is in the heap */
            if(heap.locations[index]>=0){
                
                /* check to see if a shorter distance is possible */
                double possible=fmax(eikonal_update(u, x,y),min);
                
                /*if the distance is shorter then update the data accordingly*/
                if(possible<u[index]){
                    u[index]=possible;
                    int heapLocation=heap.locations[index];
                    heap.i_heap_decrease_key(heapLocation,u[index]);
                }
            }else{
                /* the point has not yet been added  */
                
                double possible=fmax(eikonal_update(u, x,y),min);
                u[index]=possible;
                i_heap_insert_node_with_location(heap, index, u[index], index);
            }
        }
        heap.i_heap_delete_min();
        return min;
    }


    /* calculate distances until u reaches the value max.  u must contain some initial values for the code to work*/

    void run(double *u, double max){
            
        i_heap heap = i_heap(n1*n2, n1*n2);

        double min=FLT_MAX;
        
        for(int i=0;i<n1*n2;i++){
            if(u[i]<max){
                i_heap_insert_node_with_location(heap, i, u[i], i);
            }
            min=fmin(min, u[i]);
        }
        
        while(!heap.i_heap_empty()&&min<max){
            min=single_fast_marching_iteration(heap, u);
        }
    }


    /* Calculate the smoothing operator */
    double calculate_K(double s, double eps){
        return exp(-M_PI*s*s/(eps*eps));
    }

    /* calculate the perimeter */
    double calculate_perimeter(double* u){

        run(u, 1);

        double perimeter = 0;

        for(int i=0;i<n1*n2;++i){
            perimeter += calculate_K(u[i],eps);
        }

        return perimeter/eps/(n1*n2);
    }
    
};








#endif







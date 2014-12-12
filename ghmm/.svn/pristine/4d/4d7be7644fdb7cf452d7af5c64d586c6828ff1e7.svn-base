#include "block_compression.h"
#include <math.h>
#include "ghmm_internals.h"

void free_block_stats(block_stats **stats){
#define CUR_PROC "free_block_stats"
    m_free((*stats)->moment1);
    m_free((*stats)->moment2);
    m_free((*stats)->length);
    m_free((*stats));
    stats = NULL;
STOP:
    return;
#undef CUR_PROC
}

//checks if the gap is small enough to be compressed
int small_enough_gap(double* seq, int len, double size, int start, int end){
    if(end-start < 2)
        return 0; 
    double max = seq[0];
    double min = max;
    int i;
    for(i = start; i < end; i++){
        if(max < seq[i])
            max = seq[i];
        if(min > seq[i])
            min = seq[i];
    }
    //printf("max - min = %f\n", fabs(max - min));
    return fabs(max-min) < size;
}

//returns index of largest gap
int get_largest_gap(double* seq, int start, int end){
    int i, index;
    double max;
    if(end-start < 2)
        return start;
    max = fabs(seq[start] - seq[start+1]);
    index = start;
    for(i = start+1; i < end-1; i++){
        if(max < fabs(seq[i] - seq[i+1])){
            max = fabs(seq[i] - seq[i+1]);
            index = i;
        }
    }
    return index;
}

//sum of data
double get_moment1(double* seq, int start, int end){
    double moment = 0;
    int i;

    for(i = start; i < end; i++){
        moment += seq[i];
    }
    return moment;
}

//sum of sqrs of data
double get_moment2(double* seq, int start, int end){
    double moment = 0;
    int i;

    for(i = start; i < end; i++){
        moment += seq[i]*seq[i];
    }
    return moment;
}

//comparation for get_median
int comp(const void* p1, const void*p2){
    double f = *((double*)p1);
    double s = *((double*)p2);
    return (f>s)-(f<s);
}

double get_median(double *x, int start, int end){
    int len = end-start;
    double tmp[len];
    int i;
    int j = 0;
    for(i = start; i < end; i++){
        tmp[j] = x[i];
        j++;
    }
    qsort(tmp, sizeof(tmp)/sizeof(*tmp), sizeof(*tmp), comp);
    if(len % 2 == 0){
        return tmp[len/2];
    }
    else{
        return (tmp[len/2] + tmp[len/2+1])/2;
    }
}

void compress_observations_helper(double* seq, int len, double width, double delta, 
        int level, int dimension, block_stats *stats, int start, int end){

    if(end-start == 1 || small_enough_gap(seq, len, width/pow(delta, level), start, end)){
        stats->moment1[stats->total] = get_moment1(seq, start, end);
        stats->moment2[stats->total] = get_moment2(seq, start, end);
        stats->length[stats->total] = end-start;
        stats->total++;
        return;
    }
    if(dimension == 1){
        int cur_start = start, cur_end = start+1;
        double median = get_median(seq, start, end);
        while(cur_end < end){
            if(cur_end < end && seq[cur_end] <= median){
                while(cur_end < end && seq[cur_end] <= median ){
                    cur_end++;
                }
                compress_observations_helper(seq, len,  width, delta, level+1, 0,
                        stats, cur_start, cur_end);
            }
            else{
                while(cur_end < end && seq[cur_end] >= median){
                    cur_end++;
                }
                compress_observations_helper(seq, len,  width, delta, level+1, 0, 
                        stats, cur_start, cur_end);
            }
            cur_start = cur_end;
            cur_end++;
        }
        if(cur_start < end){
            if(cur_end < end && seq[cur_end] <= median){
                compress_observations_helper(seq, len, width, delta, level+1, 0,
                        stats, cur_start, cur_end);
            }
            else{
                compress_observations_helper(seq, len,  width, delta, level+1, 0, 
                                stats, cur_start, cur_end);
            }
        }
    }
    else{
        int end_index = get_largest_gap(seq, start, end);
        compress_observations_helper(seq, len, width, delta, level, 1,  
                    stats, start, end_index+1);
        compress_observations_helper(seq, len,  width, delta, level, 1,  
                    stats, end_index+1, end);
    }
}

block_stats *compress_observations(double* seq, int len, double width, double delta){
#define CUR_PROC "compress_observations"
    block_stats *stats;
    ARRAY_MALLOC(stats, 1);
    ARRAY_MALLOC(stats->moment1, len);
    ARRAY_MALLOC(stats->moment2, len);
    ARRAY_MALLOC(stats->length, len);
    stats->total = 0;

    compress_observations_helper(seq, len, width, delta, 1, 1,  
            stats, 0, len);
    return stats;

STOP:
    return NULL;
#undef CUR_PROC 
}


block_stats *merge_observations(double* seq, int len, double width,
        int max_len, block_stats *stats){
#define CUR_PROC "merge_obersvations"

    block_stats *merged_stats;

    ARRAY_MALLOC(merged_stats, stats->total);
    ARRAY_MALLOC(merged_stats->moment1, stats->total);
    ARRAY_MALLOC(merged_stats->moment2, stats->total);
    ARRAY_MALLOC(merged_stats->length, stats->total);

    int index = 0;
    int i = 0;
    int isOutlier = 0;

    double min = stats->moment1[0]/stats->length[0]+width;
    double max = stats->moment1[0]/stats->length[0]-width;

    merged_stats->moment1[0] = stats->moment1[0];
    merged_stats->moment2[0] = stats->moment2[0];
    merged_stats->length[0] = stats->length[0];
    merged_stats->total = 1;

    for(i = 1; i < stats->total; i++){
        /*
          printf("len = %d, min = %f,  max = %f, mean = %f \n", stats->length[i],
                min, max, stats->moment1[i]/stats->length[i]);
        */ 
        if(stats->length[i]  == 1 && i < stats->total -1 && stats->length[i+1] > 1 &&
                !(stats->moment1[i+1]/stats->length[i+1] < max || 
                 stats->moment1[i+1]/stats->length[i+1] > min || 
                 stats->length[i+1] + 1 + merged_stats->length[index] > max_len)){

            merged_stats->moment1[index] += stats->moment1[i];
            merged_stats->moment2[index] += stats->moment2[i];
            merged_stats->length[index] += stats->length[i];

            isOutlier = 1;
        }
        else if(isOutlier==0 && (stats->moment1[i]/stats->length[i] > min ||
                    stats->moment1[i]/stats->length[i] < max ||
                    stats->length[i] + merged_stats->length[index] > max_len) ){
            //dont merge
            index++;
            merged_stats->total++;
            merged_stats->moment1[index] = stats->moment1[i];
            merged_stats->moment2[index] = stats->moment2[i];
            merged_stats->length[index] = stats->length[i];
            max = merged_stats->moment1[index]/merged_stats->length[index] - width;
            min = merged_stats->moment1[index]/merged_stats->length[index] + width;
        }
        else{
            //merge  outlier or close consecutive block
            isOutlier = 0;
            merged_stats->length[index] += stats->length[i];
            merged_stats->moment1[index] += stats->moment1[i];
            merged_stats->moment2[index] += stats->moment2[i];
        }
    }
    free_block_stats(&stats);
    return merged_stats;
STOP:
    return NULL;
#undef CUR_PROC
}

void print_stats(block_stats *stats, int length){
    printf("total blocks = %d\n", stats->total);
    printf("obs = %d\n", length);
    printf("compression = %f\n", (double)length/(double)stats->total);
    int i;
    double mean, var;

    /*
    printf("block  length     sum                sumsqrs             mean                var \n");
    for(i = 0; i < stats->total; i++){
        mean = stats->moment1[i]/stats->length[i];
        var = stats->moment2[i]/stats->length[i] - mean*mean;
        printf("%d      %d         %e      %e      %e      %e\n",
                i, stats->length[i], 
                stats->moment1[i], stats->moment2[i], 
                mean, var);
    }
    */
}

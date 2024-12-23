// UNFINISHED

#include <algorithm>

#define INPUT_DIM 3
#define OUTPUT_DIM 1

void computeShift(int* block, int* mask, int* shift_block ,int no_of_instances){
    shift_block[0] = 0;
    for(int i=0; i<no_of_instances; i++){
        if(!mask[i]){
            // skumulowany histogram TODO: da się zrównoleglić
            shift_block[i] = shift_block[(((i-1) > (0)) ? (i-1) : (0))] + block[i+1];
        }
    }

}

int computeWeightLen(int* row_pointers, int* block, int* mask, int row_pointer_size){
    int block_idx = 0;
    int acc = 0;
    for(int i=0;i<row_pointer_size; i++){
        while(block[block_idx+1]>i) block_idx++;
        if(mask[block_idx] == 1){
            acc+=row_pointers[i+1] - row_pointers[i];
        }
    }
    return acc;
}


int testReduction(){
    int* block;
    int* mask;

    int no_of_instances;

    int* column_idx;
    int* row_pointers;
    float* weights;

    int weights_size;
    int column_idx_size;
    int row_pointer_size;

    int* shift_block; // size no_of_instances

    computeShift(block, mask, shift_block, no_of_instances);
    int final_rows_size = block[no_of_instances] - shift_block[no_of_instances-1]; // old_size - cumulated_shift
    int weight_len = computeWeightLen(row_pointers, block, mask, row_pointer_size);

    int* new_row_pointers; // size final_rows_size + 1
    float* new_weights; // size weight_len
    int* new_column_idx; // size weight_len

    new_weights = (float*) malloc(weight_len * sizeof(float));
    new_column_idx = (int*) malloc(weight_len * sizeof(int));
    new_row_pointers = (int*) malloc((final_rows_size + 1) * sizeof(int));







}

int main(){

}
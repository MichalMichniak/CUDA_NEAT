#include <stdio.h>
#include <stdlib.h>


void updateRowPointers(int *row_pointers, int *out, int *blocks_edges, int *blocks_nodes, int no_instances, bool* enable){
    int idx = 0; // granulacja na poziomie idx per kernel
    for(int i=0; i<blocks_edges[no_instances]; i++){
        while(blocks_edges[idx+1]<=i) idx++;// granica pomiędzy kernelami (warunek do pętli for)
        if (enable[i])
            row_pointers[(out[i])+blocks_nodes[idx]+1] +=1;
    }
    for(int i=0; i<blocks_nodes[no_instances]; i++){
        row_pointers[i+1] += row_pointers[i];
    }
}

void updateCol_idx_weights(int *row_pointers_t, float *weights, int *col_idx, int *in, float *w, int *out, int *blocks_edges, int *blocks_nodes, int no_instances, bool* enable){
    int idx = 0; // granulacja na poziomie idx per kernel

    for(int i=0; i<blocks_edges[no_instances]; i++){
        while(blocks_edges[idx+1]<=i) idx++;
        if (enable[i]){
            int temp = row_pointers_t[out[i] +blocks_nodes[idx]];
            row_pointers_t[out[i]+blocks_nodes[idx]] +=1; // atomicAdd
            col_idx[temp] = in[i]+blocks_nodes[idx]; // jeżeli będą podwójne to nie ma znaczenia przy mnożeniu macierzowym
            weights[temp] = w[i];
        }
    }
}
int count_col_size(bool* enable, int length){
    int acc = 0;
    for(int i=0; i<length; i++) if (enable[i]) acc++;
    return acc;
}

int testToCSR(){

    int *in;
    int *out;
    float *w;
    bool *enabled; // implemented as enable below
    int *innov; // not yet implemented
    
    int no_instances;
    int *blocks_edges;
    int *translation; // not yet implemented
    int *blocks_nodes;
    FILE *plik = fopen("COOToCSR.txt", "r"); // no_instances -> Block_nodes -> Block_edges -> in -> out -> weights
    if (plik == NULL) {
        return 1;
    }
    fscanf(plik, "%d", &no_instances);
    blocks_nodes = (int*) malloc((no_instances+1) * sizeof(int));
    for(int i = 0; i<no_instances+1; i++){
        fscanf(plik, "%d", blocks_nodes+i);
    }
    blocks_edges = (int*) malloc((no_instances+1) * sizeof(int));
    for(int i = 0; i<no_instances+1; i++){
        fscanf(plik, "%d", blocks_edges+i);
    }
    in = (int*) malloc((blocks_edges[no_instances]) * sizeof(int));
    for(int i = 0; i<blocks_edges[no_instances]; i++){
        fscanf(plik, "%d", in+i);
    }
    out = (int*) malloc((blocks_edges[no_instances]) * sizeof(int));
    for(int i = 0; i<blocks_edges[no_instances]; i++){
        fscanf(plik, "%d", out+i);
    }
    w = (float*) malloc((blocks_edges[no_instances]) * sizeof(float));
    for(int i = 0; i<blocks_edges[no_instances]; i++){
        fscanf(plik, "%f", w+i);
    }
    bool enable[8] = {false,false,true,false,true,true,true,true};
    int colsize = count_col_size(enable, 8);

    int *col_idx; // size w
    float *weights; // size w
    int *row_pointers; //size translation+1
    
    col_idx = (int*) malloc((colsize) * sizeof(int));
    weights = (float*) malloc((colsize) * sizeof(float));
    row_pointers = (int*) malloc((blocks_nodes[no_instances] + 1) * sizeof(int));
    for(int i = 0; i<(blocks_nodes[no_instances] + 1); i++){
        row_pointers[i] = 0;
    }
    updateRowPointers(row_pointers, out, blocks_edges, blocks_nodes, no_instances, enable);

    int *row_pointers_t;
    row_pointers_t = (int*) malloc((blocks_nodes[no_instances] + 1) * sizeof(int));
    for(int i=0; i<(blocks_nodes[no_instances] + 1); i++){
        row_pointers_t[i] = row_pointers[i];
    }

    updateCol_idx_weights(row_pointers_t, weights, col_idx, in, w, out, blocks_edges, blocks_nodes, no_instances,enable);

    for(int i=0; i<(blocks_nodes[no_instances] + 1); i++){
        printf("%d\t", row_pointers[i]);
    }
    printf("\n");
    for(int i=0; i<(colsize); i++){
        printf("%d\t", col_idx[i]);
    }
    printf("\n");
    for(int i=0; i<(colsize); i++){
        printf("%f\t", weights[i]);
    }

    free(blocks_nodes);
    free(blocks_edges);
    free(in);
    free(out);
    free(w);
    free(col_idx);
    free(weights);
    free(row_pointers);
    return 0;
}



int main(){
    testToCSR();
}
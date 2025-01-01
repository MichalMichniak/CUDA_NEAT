#include <stdio.h>
#include <stdlib.h>

#define BLOCK_SIZE 128
// ###### CUMULATED HISTOGRAM BEGIN #######
#define SECTION_SIZE 512

__global__ void scanKernelX(int *Y, int *S, int *X, int width)
{
    //@@ INSERT CODE HERE
    __shared__ int in[SECTION_SIZE];
    if(threadIdx.x + (blockDim.x * blockIdx.x) < width){
        in[threadIdx.x] = X[threadIdx.x + (blockDim.x * blockIdx.x)];
    }else{
        in[threadIdx.x] = 0;
    }
    __syncthreads();
    
    int result = 0;
    for(int start = 1; ((start + (blockDim.x * blockIdx.x)) < width) && (start < blockDim.x); start=start*2){
        if(threadIdx.x>=start){
            result = in[threadIdx.x] + in[threadIdx.x-start];
            // printf("it: %d\t(%d, %d)\t%f\t%f\t%f\n",it,threadIdx.x, threadIdx.x+(1<<it) ,in[threadIdx.x] ,in[threadIdx.x+(1<<it)],result);
        }
        __syncthreads();
        if(threadIdx.x>=start){
            in[threadIdx.x] = result;
        }
        //printf("%f\n",in[0]);
        __syncthreads();
    }
    
    if(threadIdx.x + (blockDim.x * blockIdx.x) < width){
        Y[threadIdx.x + (blockDim.x * blockIdx.x)] = in[threadIdx.x]; 
        if(threadIdx.x == (blockDim.x - 1)){
            S[blockIdx.x] = in[threadIdx.x];
        }
    }
    
}

__global__ void scanKernelS(int *S, int width)
{
    //@@ INSERT CODE HERE
    __shared__ int in[SECTION_SIZE];
    if(threadIdx.x + (blockDim.x * blockIdx.x) < width){
        in[threadIdx.x] = S[threadIdx.x + (blockDim.x * blockIdx.x)];
    }else{
        in[threadIdx.x] = 0;
    }
    __syncthreads();

    int result = 0;
    for(int start = 1; start < width; start=start*2){
        if(threadIdx.x>=start){
            result = in[threadIdx.x] + in[threadIdx.x-start];
            // printf("it: %d\t(%d, %d)\t%f\t%f\t%f\n",it,threadIdx.x, threadIdx.x+(1<<it) ,in[threadIdx.x] ,in[threadIdx.x+(1<<it)],result);
        }
        __syncthreads();
        if(threadIdx.x>=start){
            in[threadIdx.x] = result;
        }
        //printf("%f\n",in[0]);
        __syncthreads();
    }
    
    if(threadIdx.x < width){
        S[threadIdx.x] = in[threadIdx.x]; 
    }
}

__global__ void updateYKernel(int *Y, int *S, int widthY)
{
    //@@ INSERT CODE HERE
    if((blockIdx.x >= 1) && (threadIdx.x + (blockDim.x * blockIdx.x) < widthY)){
        Y[threadIdx.x + (blockDim.x * blockIdx.x)] += S[blockIdx.x - 1];
    }
}

void cumulatedHistogram(int *d_Y, int *d_X, int width)
{
    /*
    input:
        device vectors:
            d_X - input
            d_Y - output
        width - size of vectors
    */
    int *d_S;
    cudaMalloc(&d_S, ceil((float)width/SECTION_SIZE) * sizeof(int));

    dim3 dimGrid(ceil((float)width/SECTION_SIZE),1,1);
    dim3 dimBlock(SECTION_SIZE,1,1);
    
    scanKernelX<<<dimGrid,dimBlock>>>(d_Y, d_S, d_X, width);

    dim3 dimGrid2(1,1,1);
    dim3 dimBlock2(SECTION_SIZE,1,1);

    scanKernelS<<<dimGrid2,dimBlock2>>>(d_S, ceil((float)width/SECTION_SIZE));
    
    updateYKernel<<<dimGrid,dimBlock>>>(d_Y, d_S, width);

    cudaFree(d_S);
}
// ###### CUMULATED HISTOGRAM END #######

__global__ void countOffsprings(int *first_pair, int *second_pair, int no_offsprings, int *innov, int *blocks_edges, int no_instances, int *length_offspring, int offset)
{
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<no_offsprings){
        int first = first_pair[i];
        int second = second_pair[i];
        int idx_first = blocks_edges[first];
        int idx_second = blocks_edges[second];
        int acc = 0;
        // if(i == 0) printf("\nstart\n");
        while(idx_first != blocks_edges[first+1] && idx_second != blocks_edges[second+1]){
            // if(i == 0) printf("\ninnov: %d, %d\t %d, %d\n", innov[idx_first],innov[idx_second], idx_first, idx_second);
            if(innov[idx_first] == innov[idx_second]){
                // if(i == 0) printf("%d", innov[idx_first]);
                idx_first++;
                idx_second++;
            }else if(innov[idx_first] < innov[idx_second]){
                // if(i == 0) printf("%d\t", innov[idx_first]);
                idx_first++;
            }else{
                // if(i == 0) printf("%d\t", innov[idx_second]);
                idx_second++;
            }
            acc++;
        }
        // if(i == 0) printf("\tend\t");
        while(idx_first != blocks_edges[first+1]){
            // if(i == 0) printf("%d\t", innov[idx_second]);
            idx_first++;
            acc++;
        }
        while(idx_second != blocks_edges[second+1]){
            // if(i == 0) printf("%d\t", innov[idx_second]);
            idx_second++;
            acc++;
        }
        length_offspring[i+offset]=acc;
        // if(i == 0) printf("\n%d\n", length_offspring[i+offset]);
        // if(i == 0) printf("\nendend\n");
        i += blockDim.x * gridDim.x;
    }
}

__global__ void countOffspringsNodes(int *first_pair, int *second_pair, int no_offsprings, int *translation, int *blocks_nodes, int no_instances, int *length_offspringNodes, int offset){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<no_offsprings){
        int first = first_pair[i];
        int second = second_pair[i];
        int idx_first = blocks_nodes[first];
        int idx_second = blocks_nodes[second];
        int acc = 0;
        while(idx_first != blocks_nodes[first+1] && idx_second != blocks_nodes[second+1]){
            if(translation[idx_first] == translation[idx_second]){
                idx_first++;
                idx_second++;
            }else if(translation[idx_first] < translation[idx_second]){
                idx_first++;
            }else{
                idx_second++;
            }
            acc++;
        }
        while(idx_first != blocks_nodes[first+1]){
            idx_first++;
            acc++;
        }
        while(idx_second != blocks_nodes[second+1]){
            idx_second++;
            acc++;
        }
        length_offspringNodes[i+offset]=acc;

        i += blockDim.x * gridDim.x;
    }
}

__global__ void countSurvivors(int *blocks_edges, int *blocks_nodes, int *new_blocks_edges, int *new_blocks_nodes, int no_survivors, int *no_instance_seq){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<no_survivors){
        int instance_number = no_instance_seq[i];
        new_blocks_nodes[1+i] = blocks_nodes[instance_number+1]-blocks_nodes[instance_number];
        new_blocks_edges[1+i] = blocks_edges[instance_number+1]-blocks_edges[instance_number];
        i += blockDim.x * gridDim.x;
    }
}

__global__ void copySurvivors(int *istance_numbers_seq, int *blocks_edges, int *blocks_nodes, int *new_blocks_nodes, int *new_blocks_edges, int no_survivors,
    int *new_in, int *new_out, float *new_w, bool *new_enabled, int *new_innov, int *new_translation, 
    int *in, int *out, float *w, bool *enabled, int *innov, int *translation){

    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<no_survivors){
        // 
        for(int j=0; j<new_blocks_edges[i+1] - new_blocks_edges[i]; j++){
            new_in[j+new_blocks_edges[i]] = in[j + blocks_edges[istance_numbers_seq[i]]];
            new_out[j+new_blocks_edges[i]] = out[j + blocks_edges[istance_numbers_seq[i]]];
            new_enabled[j+new_blocks_edges[i]] = enabled[j + blocks_edges[istance_numbers_seq[i]]];
            new_innov[j+new_blocks_edges[i]] = innov[j + blocks_edges[istance_numbers_seq[i]]];
            new_w[j+new_blocks_edges[i]] = w[j + blocks_edges[istance_numbers_seq[i]]];
        }
        for(int j=0; j<new_blocks_nodes[i+1] - new_blocks_nodes[i]; j++){
            new_translation[j+new_blocks_nodes[i]] = translation[j + blocks_nodes[istance_numbers_seq[i]]];
        }
        i += blockDim.x * gridDim.x;
    }
}

__global__ void crossover(
    int *blocks_nodes,
    int *blocks_edges,
    int *new_blocks_nodes,
    int *new_blocks_edges,
    int offset,
    int *in,
    int *out,
    float *w,
    bool *enabled,
    int *innov,
    int *first_pair,
    int *second_pair,
    int no_offsprings,
    int *new_in,
    int *new_out,
    float *new_w,
    bool *new_enabled,
    int *new_innov,
    int *new_translation,
    int *translation,
    int *translation_t
    ){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<no_offsprings){
        int first = first_pair[i];
        int second = second_pair[i];
        int idx_first = blocks_nodes[first];
        int idx_second = blocks_nodes[second];
        int idx = 0;
        // nodes translations
        while(idx_first != blocks_nodes[first+1] && idx_second != blocks_nodes[second+1]){
            if(translation[idx_first] == translation[idx_second]){
                new_translation[idx + new_blocks_nodes[i+offset]] = translation[idx_first];

                translation_t[idx_first] = idx;
                translation_t[idx_second] = idx;
                idx++;
                idx_first++;
                idx_second++;
            }else if(translation[idx_first] < translation[idx_second]){
                new_translation[idx + new_blocks_nodes[i+offset]] = translation[idx_first];
                translation_t[idx_first] = idx;
                idx++;
                idx_first++;
            }else{
                new_translation[idx + new_blocks_nodes[i+offset]] = translation[idx_second];
                translation_t[idx_second] = idx;
                idx++;
                idx_second++;
            }
        }
        while(idx_first != blocks_nodes[first+1]){
            new_translation[idx + new_blocks_nodes[i+offset]] = translation[idx_first];
            translation_t[idx_first] = idx;
            idx++;
            idx_first++;
        }
        while(idx_second != blocks_nodes[second+1]){
            new_translation[idx + new_blocks_nodes[i+offset]] = translation[idx_second];
            translation_t[idx_second] = idx;
            idx++;
            idx_second++;
        }

        // edges
        idx_first = blocks_edges[first];
        idx_second = blocks_edges[second];
        idx = 0;
        // if(i==0) printf("\nSTART\n");
        while(idx_first != blocks_edges[first+1] && idx_second != blocks_edges[second+1]){
            
            if(innov[idx_first] == innov[idx_second]){ // innov - only innovation number differs edges
                new_w[idx+new_blocks_edges[i+offset]] = (w[idx_first] + w[idx_second])/2;
                new_enabled[idx+new_blocks_edges[i+offset]] = enabled[idx_first] && enabled[idx_second];
                new_innov[idx+new_blocks_edges[i+offset]] = innov[idx_first];
                // if(i==0) printf("%d, %d\t", in[idx_first]+blocks_nodes[first_pair[i]],translation_t[in[idx_first]]+blocks_nodes[first_pair[i]]);
                new_in[idx+new_blocks_edges[i+offset]] = translation_t[in[idx_first]+blocks_nodes[first_pair[i]]];
                new_out[idx+new_blocks_edges[i+offset]] = translation_t[out[idx_first]+blocks_nodes[first_pair[i]]];
                idx_first++;
                idx_second++;
                idx++;
            }else if(innov[idx_first] < innov[idx_second]){
                // if(i==0) printf("%d, %d\t", in[idx_first]+blocks_nodes[first_pair[i]],translation_t[in[idx_first]+blocks_nodes[first_pair[i]]]);
                new_w[idx+new_blocks_edges[i+offset]] = w[idx_first];
                new_enabled[idx+new_blocks_edges[i+offset]] = enabled[idx_first];
                new_innov[idx+new_blocks_edges[i+offset]] = innov[idx_first];
                new_in[idx+new_blocks_edges[i+offset]] = translation_t[in[idx_first]+blocks_nodes[first_pair[i]]];
                new_out[idx+new_blocks_edges[i+offset]] = translation_t[out[idx_first]+blocks_nodes[first_pair[i]]];
                idx_first++;
                idx++;
            }else{
                // if(i==0) printf("%d, %d\t", in[idx_second]+blocks_nodes[second_pair[i]],translation_t[in[idx_second]+blocks_nodes[second_pair[i]]]);
                new_w[idx+new_blocks_edges[i+offset]] = w[idx_second];
                new_enabled[idx+new_blocks_edges[i+offset]] = enabled[idx_second];
                new_innov[idx+new_blocks_edges[i+offset]] = innov[idx_second];
                new_in[idx+new_blocks_edges[i+offset]] = translation_t[in[idx_second]+blocks_nodes[second_pair[i]]];
                new_out[idx+new_blocks_edges[i+offset]] = translation_t[out[idx_second]+blocks_nodes[second_pair[i]]];
                idx_second++;
                idx++;
            }

        }
        while(idx_first != blocks_edges[first+1]){
            // if(i==0) printf("%d, %d\t", in[idx_first]+blocks_nodes[first_pair[i]],translation_t[in[idx_first]+blocks_nodes[first_pair[i]]]);
            new_w[idx+new_blocks_edges[i+offset]] = w[idx_first];
            new_enabled[idx+new_blocks_edges[i+offset]] = enabled[idx_first];
            new_innov[idx+new_blocks_edges[i+offset]] = innov[idx_first];
            new_in[idx+new_blocks_edges[i+offset]] = translation_t[in[idx_first]+blocks_nodes[first_pair[i]]];
            new_out[idx+new_blocks_edges[i+offset]] = translation_t[out[idx_first]+blocks_nodes[first_pair[i]]];
            idx_first++;
            idx++;
        }
        while(idx_second != blocks_edges[second+1]){
            // if(i==0) printf("%d, %d\t", in[idx_second]+blocks_nodes[second_pair[i]],translation_t[in[idx_second]+blocks_nodes[second_pair[i]]]);
            new_w[idx+new_blocks_edges[i+offset]] = w[idx_second];
            new_enabled[idx+new_blocks_edges[i+offset]] = enabled[idx_second];
            new_innov[idx+new_blocks_edges[i+offset]] = innov[idx_second];
            new_in[idx+new_blocks_edges[i+offset]] = translation_t[in[idx_second]+blocks_nodes[second_pair[i]]];
            new_out[idx+new_blocks_edges[i+offset]] = translation_t[out[idx_second]+blocks_nodes[second_pair[i]]];
            idx_second++;
            idx++;
        }
        // if(i==0) printf("\nEND\n");
        i += blockDim.x * gridDim.x;
    }
    
}


void test_crosover(){
    FILE *plik = fopen("crosoverCOO_test.txt", "r");
    if (plik == NULL) {
        return;
    }
    /*
    ########## wektory Populacji wejściowej ##########
    int no_instances - ilość instancji wejściowych
    int *blocks_nodes - hstogram skumulowany ilości wierzchołków (node) instancji (zaczynający się od 0) [0, end_1 + 1, end_1 + end_2 + 1, ...] (długość no_instances+1)
    int *translation - mapowanie odcinków liczb naturalnych [0,n] do rzeczywistych numerów wierzchołków   (zał że w instancjach posortowane rosnąco) (długość blocks_nodes[no_instances])
    int *blocks_edges - hstogram skumulowany ilości krawędzi (edges) instancji (długość no_instances+1)
    int *in - wejścia synaps (krawędzi/edges) instancji zmapowane do odcinka [0,n] (długość blocks_edge[no_instances])
    int *out - wejścia synaps (krawędzi/edges) instancji zmapowane do odcinka [0,n] (długość blocks_edge[no_instances])
    float *w - wagi synaps (krawędzi/edges) instancji
    bool *enabled - czy dane krawędzie są enabled (true jeżeli biorą udział w ewaluacji, false jeżeli nie biorą udziału w ewaluacji)
    int *innov - innovation number unikatowy numer rozróżniający krawędzie pomiędzy genami
    
    
    ########## wektory wejściowe z algorytmu genetycznego ##########
    int *mask - bitowa maska 0 znaczy że instancja nie przechodzi do następnej populacji 1 że przechodzi // mask which survives
    int no_survivors - ilość instancji która przetrwa (ilość jedynek w int *mask)
    int no_mutations - ilość instancji z mutacji
    int no_offsprings - ilość instancji z krzyżowania
    int *first_pair - numer pierwszego rodzica każdy index odpowiada jednemu potomkowi(długość no_offsprings)
    int *second_pair - numer drugiego rodzica każdy index odpowiada jednemu potomkowi(długość no_offsprings)
    
    ########## wektory Populacji wyjściowej (po definicje patrz "wektory Populacji wejściowej") ########## (To mamy zwrócić)
    int *new_blocks_nodes - długość no_survivors + no_offsprings + no_mutations
    int *new_blocks_edges - długość no_survivors + no_offsprings + no_mutations
    int *new_in
    int *new_out 
    float *new_w 
    bool *new_enabled 
    int *new_innov 
    int new_no_instances - długość survivors + offsprings + mutated
    int *new_translation 

    ######### Wektory pomocnicze ##########
    int *istance_numbers_seq - funkcja LUT do mapowania [0,długość survivors) w numery instancji które przechodzą do następnej populacji
    int *translation_t - tymczasowa funkcja mapująca stare indexy w nowe używana w krzyżowaniu (długość translation)
    */

    int *in;
    int *out;
    float *w;
    bool *enabled;
    int *innov; 

    int no_instances;
    int *blocks_edges;
    int *translation; // zał że w instancjach posortowane rosnąco
    int *blocks_nodes;
    int *mask; // mask which survives
    int no_survivors; // number of 1 in mask
    int no_mutations = 0; // TODO: mutacia

    int no_offsprings; // ze starej populacji
    int *first_pair; // first parent: size no_offsprings
    int *second_pair; // second parent: size no_offsprings
    // init
    fscanf(plik, "%d", &no_instances);
    blocks_nodes = (int*) malloc((no_instances+1) * sizeof(int));
    for(int i = 0; i<no_instances+1; i++){
        fscanf(plik, "%d", blocks_nodes+i);
    }
    blocks_edges = (int*) malloc((no_instances+1) * sizeof(int));
    for(int i = 0; i<no_instances+1; i++){
        fscanf(plik, "%d", blocks_edges+i);
    }
    translation = (int*) malloc((blocks_nodes[no_instances]) * sizeof(int));
    for(int i = 0; i<blocks_nodes[no_instances]; i++){
        fscanf(plik, "%d", translation+i);
    }
    innov = (int*) malloc((blocks_edges[no_instances]) * sizeof(int));
    for(int i = 0; i<blocks_edges[no_instances]; i++){
        fscanf(plik, "%d", innov+i);
    }
    enabled = (bool*) malloc((blocks_edges[no_instances]) * sizeof(int));
    for(int i = 0; i<blocks_edges[no_instances]; i++){
        int temp;
        fscanf(plik, "%d", &temp);
        *(enabled+i) = (bool)temp;
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
    fscanf(plik, "%d", &no_survivors);
    mask = (int*) malloc(no_instances * sizeof(int));
    for(int i = 0; i<no_instances; i++){
        fscanf(plik, "%d", mask+i);
    }
    fscanf(plik, "%d", &no_offsprings);
    first_pair = (int*) malloc(no_offsprings * sizeof(int));
    for(int i = 0; i<no_offsprings; i++){
        fscanf(plik, "%d", first_pair+i);
    }
    second_pair = (int*) malloc(no_offsprings * sizeof(int));
    for(int i = 0; i<no_offsprings; i++){
        fscanf(plik, "%d", second_pair+i);
    }
    // end of init
    // printf("\nold block nodes: ");

    // for(int i=0; i<6; i++){
    //     printf("%d\t", blocks_nodes[i]);
    // }
    // printf("\nold block edges: ");
    // for(int i=0; i<6; i++){
    //     printf("%d\t", blocks_edges[i]);
    // }
    // printf("\ninnovation numbers: ");
    // for(int i=0; i<24; i++){
    //     printf("%d\t", innov[i]);
    // }
    // printf("\n");
    // ###### inicializacja i alokacja ######
    int *d_in;
    int *d_out;
    float *d_w;
    bool *d_enabled;
    int *d_innov;
    
    int *d_blocks_edges;
    int *d_translation;
    int *d_blocks_nodes;

    int *d_mask;
    int *d_first_pair;
    int *d_second_pair;
    // alocation
    cudaMalloc(&d_in, (blocks_edges[no_instances]) * sizeof(int));
    cudaMalloc(&d_out, (blocks_edges[no_instances]) * sizeof(int));
    cudaMalloc(&d_w, (blocks_edges[no_instances]) * sizeof(float));
    cudaMalloc(&d_enabled, (blocks_edges[no_instances]) * sizeof(bool));
    cudaMalloc(&d_innov, (blocks_edges[no_instances]) * sizeof(int));

    cudaMalloc(&d_blocks_edges, (no_instances+1) * sizeof(int));
    cudaMalloc(&d_translation, (blocks_nodes[no_instances]) * sizeof(int));
    cudaMalloc(&d_blocks_nodes, (no_instances+1) * sizeof(int));

    cudaMalloc(&d_mask, no_instances * sizeof(int));
    cudaMalloc(&d_first_pair, no_offsprings * sizeof(int));
    cudaMalloc(&d_second_pair, no_offsprings * sizeof(int));

    // data copy
    cudaMemcpy(d_in, in, (blocks_edges[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_out, out, (blocks_edges[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_w, w, (blocks_edges[no_instances]) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_enabled, enabled, (blocks_edges[no_instances]) * sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(d_innov, innov, (blocks_edges[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_blocks_edges, blocks_edges, (no_instances+1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_translation, translation, (blocks_nodes[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_blocks_nodes, blocks_nodes, (no_instances+1) * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_mask, mask, no_instances * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_first_pair, first_pair, no_offsprings * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_second_pair, second_pair, no_offsprings * sizeof(int), cudaMemcpyHostToDevice);

    //###### inicializacja i alokacja END ######
    // ##### to do kopiowania #####
    int *h_mask;
    h_mask = (int*) malloc(no_instances * sizeof(int));
    cudaMemcpy(h_mask, d_mask, no_instances * sizeof(int), cudaMemcpyDeviceToHost);

    int *istance_numbers_seq; // tablica numeru instancji mapowanej ze starej do nowej tablicy szie: no_survivors
    istance_numbers_seq = (int*)malloc(no_survivors*sizeof(int));
    int idx = 0;
    for(int i=0; i<no_instances; i++){ // to będzie sekwencyjne bo tylko przepisanie
        if(h_mask[i] == 1){
            istance_numbers_seq[idx] = i;
            idx++;
        }
    }
    int *d_istance_numbers_seq;
    cudaMalloc(&d_istance_numbers_seq, no_survivors * sizeof(int));
    cudaMemcpy(d_istance_numbers_seq, istance_numbers_seq, no_survivors * sizeof(int), cudaMemcpyHostToDevice);

    int *d_new_blocks_nodes;
    int *d_new_blocks_edges;
    cudaMalloc(&d_new_blocks_nodes, (1+no_survivors+no_offsprings+no_mutations) * sizeof(int));
    cudaMalloc(&d_new_blocks_edges, (1+no_survivors+no_offsprings+no_mutations) * sizeof(int));

    cudaMemset(d_new_blocks_nodes, 0, (1+no_survivors+no_offsprings+no_mutations) * sizeof(int));
    cudaMemset(d_new_blocks_edges, 0, (1+no_survivors+no_offsprings+no_mutations) * sizeof(int));

    dim3 dimGridSurvivors(ceil((float)(no_survivors)/BLOCK_SIZE),1,1);
    dim3 dimBlockSurvivors(BLOCK_SIZE,1,1);

    dim3 dimGridOffsprings(ceil((float)(no_offsprings)/BLOCK_SIZE),1,1);
    dim3 dimBlockOffsprings(BLOCK_SIZE,1,1);

    countSurvivors<<<dimGridSurvivors, dimGridOffsprings>>>(d_blocks_edges, d_blocks_nodes, d_new_blocks_edges, d_new_blocks_nodes, no_survivors, d_istance_numbers_seq); // uzupełnienie rozmiarów instancji survivors  (0 | 1 do no_survivors)

    countOffsprings<<<dimBlockOffsprings, dimBlockOffsprings>>>(d_first_pair, d_second_pair, no_offsprings, d_innov, d_blocks_edges, no_instances, d_new_blocks_edges, no_survivors + 1);
    
    countOffspringsNodes<<<dimBlockOffsprings, dimBlockOffsprings>>>(d_first_pair, d_second_pair, no_offsprings, d_translation, d_blocks_nodes, no_instances, d_new_blocks_nodes, no_survivors + 1);
    // TODO: mutations
    
    // end mutations
    dim3 dimGrid(ceil((float)((1+no_survivors+no_offsprings+no_mutations))/BLOCK_SIZE),1,1);
    dim3 dimBlock(BLOCK_SIZE,1,1);
    cumulatedHistogram(d_new_blocks_nodes, d_new_blocks_nodes, (1+no_survivors+no_offsprings+no_mutations));
    cumulatedHistogram(d_new_blocks_edges, d_new_blocks_edges, (1+no_survivors+no_offsprings+no_mutations));

    int *d_new_in;
    int *d_new_out;
    float *d_new_w;
    bool *d_new_enabled;
    int *d_new_innov;
    int new_no_instances; // survivors + offsprings + mutated
    int *d_new_translation; // 
    new_no_instances = no_survivors + no_offsprings + no_mutations;

    int no_edges;
    int no_nodes;

    cudaMemcpy(&no_edges, d_new_blocks_edges + new_no_instances, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&no_nodes, d_new_blocks_nodes + new_no_instances, sizeof(int), cudaMemcpyDeviceToHost);

    cudaMalloc(&d_new_in, no_edges * sizeof(int));
    cudaMalloc(&d_new_out, no_edges * sizeof(int));
    cudaMalloc(&d_new_w, no_edges * sizeof(float));
    cudaMalloc(&d_new_enabled, no_edges * sizeof(bool));
    cudaMalloc(&d_new_innov, no_edges * sizeof(int));
    cudaMalloc(&d_new_innov, no_edges * sizeof(int));
    cudaMalloc(&d_new_translation, no_nodes * sizeof(int));
    
    // tytaj będzie przepisywanie survivorsów
    copySurvivors<<<dimGridSurvivors, dimGridOffsprings>>>(d_istance_numbers_seq, d_blocks_edges, d_blocks_nodes, d_new_blocks_nodes, d_new_blocks_edges, no_survivors,
    d_new_in, d_new_out, d_new_w, d_new_enabled, d_new_innov, d_new_translation, d_in, d_out, d_w, d_enabled, d_innov, d_translation);
    
    // tytaj będzie tworzenie potomków
    int *d_translation_t;
    cudaMalloc(&d_translation_t, no_nodes * sizeof(int));
    
    crossover<<<dimBlockOffsprings, dimBlockOffsprings>>>(
    d_blocks_nodes,
    d_blocks_edges,
    d_new_blocks_nodes,
    d_new_blocks_edges,
    no_survivors, // offset
    d_in,
    d_out,
    d_w,
    d_enabled,
    d_innov,
    d_first_pair,
    d_second_pair,
    no_offsprings,
    d_new_in,
    d_new_out,
    d_new_w,
    d_new_enabled,
    d_new_innov,
    d_new_translation,
    d_translation,
    d_translation_t 
    );

    // podmianka
    cudaFree(d_in);
    cudaFree(d_out);
    cudaFree(d_w);
    cudaFree(d_enabled);
    cudaFree(d_innov);
    cudaFree(d_blocks_edges);
    cudaFree(d_translation);
    cudaFree(d_blocks_nodes);
    cudaFree(d_mask);
    cudaFree(d_first_pair);
    cudaFree(d_second_pair);
    d_in=d_new_in;
    d_out=d_new_out;
    d_w=d_new_w;
    d_enabled=d_new_enabled;
    d_innov=d_new_innov;

    d_blocks_edges=d_new_blocks_edges;
    d_translation=d_new_translation;
    d_blocks_nodes=d_new_blocks_nodes;

    // czyszczenie
    cudaFree(d_translation_t);
    cudaFree(d_istance_numbers_seq);
    free(istance_numbers_seq);
    free(h_mask);

    // end czyszczenie

    // sprawdzenie:
    int *new_in;
    int *new_out;
    float *new_w;
    bool *new_enabled;
    int *new_innov;

    int *new_translation; // 
    int *new_blocks_nodes;
    int *new_blocks_edges;
    new_no_instances = no_survivors + no_offsprings + no_mutations;

    new_blocks_nodes = (int*)malloc((1+no_survivors+no_offsprings+no_mutations)*sizeof(int)); // 0 | 1 do no_survivors | no_survivors + 1 do no_survivors + no_offsprings | no_survivors + no_offsprings + 1 do no_survivors + no_offsprings + no_mutation 
    new_blocks_edges = (int*)malloc((1+no_survivors+no_offsprings+no_mutations)*sizeof(int)); 
    new_in = (int*)malloc(no_edges*sizeof(int));
    new_out = (int*)malloc(no_edges*sizeof(int));
    new_w = (float*)malloc(no_edges*sizeof(float));
    new_enabled = (bool*)malloc(no_edges*sizeof(bool));
    new_innov = (int*)malloc(no_edges*sizeof(int));
    new_translation = (int*)malloc(no_nodes*sizeof(int));

    cudaMemcpy(new_in, d_in, no_edges * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(new_out, d_out, no_edges * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(new_w, d_w, no_edges * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(new_enabled, d_enabled, no_edges * sizeof(bool), cudaMemcpyDeviceToHost);
    cudaMemcpy(new_innov, d_innov, no_edges * sizeof(int), cudaMemcpyDeviceToHost);

    cudaMemcpy(new_blocks_edges, d_blocks_edges, (new_no_instances+1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(new_translation, d_translation, no_nodes * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(new_blocks_nodes, d_blocks_nodes, (new_no_instances+1) * sizeof(int), cudaMemcpyDeviceToHost);



    printf("\nnew block nodes: ");

    for(int i=0; i<new_no_instances+1; i++){
        printf("%d\t", new_blocks_nodes[i]);
    }
    printf("\nnew block edges: ");
    for(int i=0; i<new_no_instances+1; i++){
        printf("%d\t", new_blocks_edges[i]);
    }
    // printf("\ntemp translation: ");
    // for(int i=0; i<no_nodes; i++){
    //     printf("%d\t", new_translation_t[i]);
    // }
    printf("\nnew translation: ");
    for(int i=0; i<no_nodes; i++){
        printf("%d\t", new_translation[i]);
    }
    printf("\nnew in: ");
    for(int i=0; i<no_edges; i++){
        printf("%d\t", new_in[i]);
    }

    printf("\nnew out: ");
    for(int i=0; i<no_edges; i++){
        printf("%d\t", new_out[i]);
    }

    printf("\nnew innov: ");
    for(int i=0; i<no_edges; i++){
        printf("%d\t", new_innov[i]);
    }

    printf("\nnew enabled: ");
    for(int i=0; i<no_edges; i++){
        printf("%d\t", (int)new_enabled[i]);
    }
    printf("\n");
    // tutaj free TODO:

}




int main(){
    test_crosover();
}
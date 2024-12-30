#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstdlib>
#include <cmath>

#define WIDTH 800
#define HEIGHT 600
#define GAP_R 140
#define FLOPPY_RADIUS 10

#define BLOCK_SIZE 128
// ###### CUMULATED HISTOGRAM BEGIN #######
#define SECTION_SIZE 512

__global__ void scanKernelX(int *Y, int *Y_copy, int *S, int *X, int width)
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
        Y_copy[threadIdx.x + (blockDim.x * blockIdx.x)] = in[threadIdx.x];
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

__global__ void updateYKernel(int *Y, int *Y_copy, int *S, int widthY)
{
    //@@ INSERT CODE HERE
    if((blockIdx.x >= 1) && (threadIdx.x + (blockDim.x * blockIdx.x) < widthY)){
        Y[threadIdx.x + (blockDim.x * blockIdx.x)] += S[blockIdx.x - 1];
        Y_copy[threadIdx.x + (blockDim.x * blockIdx.x)] = Y[threadIdx.x + (blockDim.x * blockIdx.x)];
    }
}

void cumulatedHistogram_with_copy(int *d_Y, int *d_Y_copy, int *d_X, int width)
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
    
    scanKernelX<<<dimGrid,dimBlock>>>(d_Y, d_Y_copy, d_S, d_X, width);

    dim3 dimGrid2(1,1,1);
    dim3 dimBlock2(SECTION_SIZE,1,1);

    scanKernelS<<<dimGrid2,dimBlock2>>>(d_S, ceil((float)width/SECTION_SIZE));
    
    updateYKernel<<<dimGrid,dimBlock>>>(d_Y, d_Y_copy, d_S, width);

    cudaFree(d_S);
}
// ###### CUMULATED HISTOGRAM END #######
// ###### Reduction BEGIN #######

float reductionSequential(float *input, int width)
{
    float sum = 0.0f;
    for (int i = 0; i < width; ++i)
    {
        sum += input[i];
    }

    return sum;
}

__global__ void reductionKernelOp(float *input, float *output, int width)
{
    //@@ INSERT CODE HERE
    __shared__ float in[BLOCK_SIZE];
    if(threadIdx.x + (blockDim.x * blockIdx.x) < width){
        in[threadIdx.x] = input[threadIdx.x + (blockDim.x * blockIdx.x)];
    }else{
        in[threadIdx.x] = 0;
    }
    __syncthreads();
    int targetlevel = 0;
    int index = BLOCK_SIZE;
    while (index >>= 1) ++targetlevel;

    float result = 0;
    for(int it = 0; it < targetlevel; it++){
        if(threadIdx.x<(BLOCK_SIZE/(2<<it))){
            result = in[threadIdx.x] + in[threadIdx.x+(BLOCK_SIZE/(2<<it))];
            // printf("it: %d\t(%d, %d)\t%f\t%f\t%f\n",it,threadIdx.x, threadIdx.x+(1<<it) ,in[threadIdx.x] ,in[threadIdx.x+(1<<it)],result);
        }
        __syncthreads();
        if((threadIdx.x<(BLOCK_SIZE/(2<<it)))){
            in[threadIdx.x] = result;
        }
        //printf("%f\n",in[0]);
        __syncthreads();
    }
    
    if(threadIdx.x == 0){
        output[blockIdx.x] = in[0]; 
        // printf("%f",in[0]);
        // printf("\n%d\n",targetlevel);
    }
    
}

float launchReductionKernelOp(float *h_input, int width)
{
    //@@ INSERT CODE HERE
    float *d_input, *d_output, *h_output;
    cudaMalloc(&d_input, width * sizeof(float));
    cudaMalloc(&d_output, ceil((float)width/BLOCK_SIZE) * sizeof(float));
    h_output = (float*)malloc(ceil((float)width/BLOCK_SIZE) * sizeof(float));

    cudaMemcpy(d_input, h_input, width*sizeof(float), cudaMemcpyHostToDevice);

    dim3 dimGrid(ceil((float)width/BLOCK_SIZE),1,1);
    dim3 dimBlock(BLOCK_SIZE,1,1);
    
    reductionKernelOp<<<dimGrid,dimBlock>>>(d_input,d_output,width);
    
    cudaMemcpy(h_output, d_output,  ceil((float)width/BLOCK_SIZE)*sizeof(float), cudaMemcpyDeviceToHost);

    float result = reductionSequential(h_output,ceil((float)width/BLOCK_SIZE));
    cudaFree(d_output);
    cudaFree(d_input);
    
    free(h_output);
    return result;
}
// ###### Reduction END #######

// ###### Reduction boolean count BEGIN #######

int reductionSequential_count_col_size(int *input, int width)
{
    int sum = 0;
    for (int i = 0; i < width; ++i)
    {
        sum += input[i];
    }

    return sum;
}

__global__ void reductionKernelOp_count_col_size(bool *input, int *output, int width)
{
    //@@ INSERT CODE HERE
    __shared__ int in[BLOCK_SIZE];
    if(threadIdx.x + (blockDim.x * blockIdx.x) < width){
        in[threadIdx.x] = (int)(input[threadIdx.x + (blockDim.x * blockIdx.x)]);
    }else{
        in[threadIdx.x] = 0;
    }
    __syncthreads();
    int targetlevel = 0;
    int index = BLOCK_SIZE;
    while (index >>= 1) ++targetlevel;

    int result = 0;
    for(int it = 0; it < targetlevel; it++){
        if(threadIdx.x<(BLOCK_SIZE/(2<<it))){
            result = in[threadIdx.x] + in[threadIdx.x+(BLOCK_SIZE/(2<<it))];
            // printf("it: %d\t(%d, %d)\t%f\t%f\t%f\n",it,threadIdx.x, threadIdx.x+(1<<it) ,in[threadIdx.x] ,in[threadIdx.x+(1<<it)],result);
        }
        __syncthreads();
        if((threadIdx.x<(BLOCK_SIZE/(2<<it)))){
            in[threadIdx.x] = result;
        }
        //printf("%f\n",in[0]);
        __syncthreads();
    }
    
    if(threadIdx.x == 0){
        output[blockIdx.x] = in[0]; 
        // printf("%f",in[0]);
        // printf("\n%d\n",targetlevel);
    }
    
}

int count_col_size(bool* d_enable, int length)
{
    //@@ INSERT CODE HERE
    int *d_output, *h_output;
    cudaMalloc(&d_output, ceil((float)length/BLOCK_SIZE) * sizeof(int));
    h_output = (int*)malloc(ceil((float)length/BLOCK_SIZE) * sizeof(int));

    dim3 dimGrid(ceil((float)length/BLOCK_SIZE),1,1);
    dim3 dimBlock(BLOCK_SIZE,1,1);
    
    reductionKernelOp_count_col_size<<<dimGrid,dimBlock>>>(d_enable,d_output,length);
    
    cudaMemcpy(h_output, d_output,  ceil((float)length/BLOCK_SIZE)*sizeof(int), cudaMemcpyDeviceToHost);

    int result = reductionSequential_count_col_size(h_output,ceil((float)length/BLOCK_SIZE));
    cudaFree(d_output);
    
    free(h_output);
    return result;
}
// ###### Reduction END #######

// ###### Reduction boolean count BEGIN #######

int reductionSequential_STOP_cryterion(bool *input, int width)
{
    bool sum = false;
    for (int i = 0; i < width; ++i)
    {
        sum = sum || input[i];
    }

    return sum;
}

__global__ void reductionKernelOp_STOP_cryterion(bool *input, bool *output, int width)
{
    //@@ INSERT CODE HERE
    __shared__ bool in[BLOCK_SIZE];
    if(threadIdx.x + (blockDim.x * blockIdx.x) < width){
        in[threadIdx.x] = !(input[threadIdx.x + (blockDim.x * blockIdx.x)]);
    }else{
        in[threadIdx.x] = false;
    }
    __syncthreads();
    int targetlevel = 0;
    int index = BLOCK_SIZE;
    while (index >>= 1) ++targetlevel;

    bool result = false;
    for(int it = 0; it < targetlevel; it++){
        if(threadIdx.x<(BLOCK_SIZE/(2<<it))){
            result = in[threadIdx.x] || in[threadIdx.x+(BLOCK_SIZE/(2<<it))];
            // printf("it: %d\t(%d, %d)\t%f\t%f\t%f\n",it,threadIdx.x, threadIdx.x+(1<<it) ,in[threadIdx.x] ,in[threadIdx.x+(1<<it)],result);
        }
        __syncthreads();
        if((threadIdx.x<(BLOCK_SIZE/(2<<it)))){
            in[threadIdx.x] = result;
        }
        //printf("%f\n",in[0]);
        __syncthreads();
    }
    
    if(threadIdx.x == 0){
        output[blockIdx.x] = in[0]; 
        // printf("%f",in[0]);
        // printf("\n%d\n",targetlevel);
    }
    
}

bool STOP_cryterion(bool* d_collision, int length)
{
    //@@ INSERT CODE HERE
    bool *d_output, *h_output;
    cudaMalloc(&d_output, ceil((float)length/BLOCK_SIZE) * sizeof(bool));
    h_output = (bool*)malloc(ceil((float)length/BLOCK_SIZE) * sizeof(bool));

    dim3 dimGrid(ceil((float)length/BLOCK_SIZE),1,1);
    dim3 dimBlock(BLOCK_SIZE,1,1);
    
    reductionKernelOp_STOP_cryterion<<<dimGrid,dimBlock>>>(d_collision,d_output,length);
    
    cudaMemcpy(h_output, d_output,  ceil((float)length/BLOCK_SIZE)*sizeof(bool), cudaMemcpyDeviceToHost);

    bool result = reductionSequential_STOP_cryterion(h_output,ceil((float)length/BLOCK_SIZE));
    cudaFree(d_output);
    
    free(h_output);
    return result;
}
// ###### Reduction END #######

// Mnożenie macierzy CSR przez wektor 
__global__ void SparseMUL(int* column_idx, int* row_pointers, float* weights, float* input_vector, int vector_size, float* output_vector, int output_vector_size){
    int row = threadIdx.x + (blockDim.x * blockIdx.x);
    while(row<output_vector_size){
        float acc = 0;
        for(int col_idx=row_pointers[row]; col_idx<row_pointers[row+1]; col_idx++){
            acc += input_vector[column_idx[col_idx]] * weights[col_idx];
        }
        output_vector[row] = 1/(1+exp(-acc)); // sigmoida
        
        row += blockDim.x * gridDim.x;
    }
}

__global__ void updateRowPointers(int *row_pointers, int *out, int *blocks_edges, int *blocks_nodes, int no_instances, bool* enable){
    int inst = threadIdx.x + (blockDim.x * blockIdx.x);
    while(inst<no_instances){
        int idx = inst; // granulacja na poziomie idx per kernel
        for(int i=blocks_edges[idx]; i<blocks_edges[idx+1]; i++){
            if (enable[i])
                row_pointers[(out[i])+blocks_nodes[idx]+1] +=1;
        }
        inst += blockDim.x * gridDim.x;
    }
    
}

__global__ void updateCol_idx_weights(int *row_pointers_t, float *weights, int *col_idx, int *in, float *w, int *out, int *blocks_edges, int *blocks_nodes, int no_instances, bool* enable){
    int inst = threadIdx.x + (blockDim.x * blockIdx.x);
    while(inst<no_instances){
        int idx = inst;
        for(int i=blocks_edges[inst]; i<blocks_edges[inst+1]; i++){
            if (enable[i]){
                int temp = row_pointers_t[out[i] +blocks_nodes[idx]];
                row_pointers_t[out[i]+blocks_nodes[idx]] +=1; // no need for atomic add because [blocks_nodes[idx], blocks_nodes[idx]+1] is only for idx
                col_idx[temp] = in[i]+blocks_nodes[idx]; // jeżeli będą podwójne to nie ma znaczenia przy mnożeniu macierzowym
                weights[temp] = w[i];
            }
        }
        inst += blockDim.x * gridDim.x;
    }
}



class Pipe{
private:
    int y_upper_pipe;
    int y_bottom_pipe;
    int x_pipe;
    int r_pipe;
    int vel;

public:
    Pipe(int y_upper_pipe_, int y_bottom_pipe_, int x_pipe_, int r_pipe_, int vel_) : y_bottom_pipe(y_bottom_pipe_), y_upper_pipe(y_upper_pipe_), x_pipe(x_pipe_), r_pipe(r_pipe_), vel(vel_) {};

    Pipe(const Pipe& other) 
        : y_upper_pipe(other.y_upper_pipe), y_bottom_pipe(other.y_bottom_pipe),
          x_pipe(other.x_pipe), r_pipe(other.r_pipe), vel(other.vel) {
    }
    bool update(){
        x_pipe -= vel;
        if(x_pipe > -10) return 0;
        return 1;
    };
    int getX() const { return x_pipe; }
    int getYUpper() const { return y_upper_pipe; }
    int getYBottom() const { return y_bottom_pipe; }
    int getR() const { return r_pipe; }
    int getVelocity() const { return vel; }
    ~Pipe() = default;
};


__global__ void update_colisions(bool* colision, int* rewards, float* y, float x, int y_upper_pipe, int y_bottom_pipe, int x_pipe, int r_pipe, int no_instances){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<no_instances){
        // printf("%f\n",y[i]);
        if((((x + FLOPPY_RADIUS >=(x_pipe-r_pipe)) && (x-FLOPPY_RADIUS)<=x_pipe+r_pipe) && (y[i]-FLOPPY_RADIUS<=y_upper_pipe || y[i]+FLOPPY_RADIUS>=HEIGHT-y_bottom_pipe)) || (y[i]>=HEIGHT-FLOPPY_RADIUS || y[i]<=0+FLOPPY_RADIUS)){
            colision[i] = true;
        }
        else{
            if(!colision[i]){
                rewards[i]++;
                
            }
        }
        i += blockDim.x * gridDim.x;
    }
}

__global__ void update_step(float* vel_y, float* y, float* vect_out,int *blocks_nodes, int no_instances){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<no_instances){
        float g = 0.1;
        float flop_v = 4;
        vel_y[i] += g;
        if(vect_out[blocks_nodes[i]]>0.5){
            vel_y[i] = -flop_v;
        }
        y[i]+=vel_y[i];
        i += blockDim.x * gridDim.x;
    }
}

__global__ void update_instance_inputs(int *blocks_nodes, float* vect_in, int y_upper_pipe, int y_bottom_pipe, int x_pipe, float x, float* y, int no_instances){
    int i = threadIdx.x + (blockDim.x * blockIdx.x); // zrównoleglenie na instancje
    while(i<no_instances){
        float x1 = (x_pipe - x)/400.0;
        float x2 = (y_bottom_pipe - y[i])/600.0;
        float x3 = (y[i] - y_upper_pipe)/600.0;
        vect_in[blocks_nodes[i] + 1] = x1;
        vect_in[blocks_nodes[i] + 2] = x2;
        vect_in[blocks_nodes[i] + 3] = x3;
        // printf("%d\t%d\t%f\t%f\t%f\n",i, blocks_nodes[i],vect_in[blocks_nodes[i] + 1],vect_in[blocks_nodes[i] + 2],vect_in[blocks_nodes[i] + 3]);
        i += blockDim.x * gridDim.x;
    }
    
}

__global__ void initialize(float* array, float value, int size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        array[idx] = value;
    }
}


void symulate(int *blocks_nodes, int* column_idx, int* row_pointers, float* weights, int* rewards, int vector_size, int no_instances, int max_iterations, int no_nodes, int no_edges){
    // inicjalizacja
    Pipe now = Pipe(300,HEIGHT-(300 + 2*GAP_R),440,40,1);
    Pipe next = Pipe(120,HEIGHT-(120 + 2*GAP_R),880,40,1);
    Pipe prev = Pipe(300,HEIGHT-(300 + 2*GAP_R),440,40,1);
    // bool* isFlopping;
    // isFlopping = (bool*) malloc((no_instances) * sizeof(bool));
    float x = 120; // takie same x dla populaci
    float* d_y;
    dim3 dimGrid(ceil((float)(no_instances)/BLOCK_SIZE),1,1);
    dim3 dimBlock(BLOCK_SIZE,1,1);
    cudaMalloc(&d_y, (no_instances) * sizeof(float));
    initialize<<<dimGrid, dimBlock>>>(d_y, 50.0f, no_instances);

    float* d_vel_y;
    cudaMalloc(&d_vel_y, (no_instances) * sizeof(float));
    cudaMemset(d_vel_y, 0.0, (no_instances) * sizeof(float));

    bool* d_colision;// inicialize to false
    cudaMalloc(&d_colision, (no_instances) * sizeof(bool));
    cudaMemset(d_colision, false, (no_instances) * sizeof(bool));

    // internal state vectors:
    float* d_vect_in;
    float* d_vect_out;
    float* d_vect_temp;
    

    cudaMalloc(&d_vect_in, (no_nodes) * sizeof(float));
    cudaMemset(d_vect_in, 0.0, (no_nodes) * sizeof(float));
    cudaMalloc(&d_vect_out, (no_nodes) * sizeof(float));
    cudaMemset(d_vect_out, 0.0, (no_nodes) * sizeof(float));
    // koniec inicjalizacji
    
    for(int it = 0; it < max_iterations; it++){ // główna pętla
        
        update_colisions<<<dimGrid,dimBlock>>>(d_colision, rewards, d_y, x, now.getYUpper(), now.getYBottom(), now.getX(), now.getR(), no_instances);
        bool STOP = STOP_cryterion(d_colision, no_instances);
        // printf("%d\t%d\n", STOP, it);
        if(!STOP) break;
        // #### UPDATE ENV ####
        if(x - FLOPPY_RADIUS > now.getX() + now.getR()){
            prev = now;
            now = next;
            int randomNum = rand() %(HEIGHT - 2*GAP_R) + GAP_R;
            next = Pipe(randomNum - GAP_R,HEIGHT - (randomNum + GAP_R),880,40,1);
        }
        // #### UPDATE ENV END ####
        // #### Control ####
        d_vect_temp = d_vect_in;
        d_vect_in = d_vect_out;
        d_vect_out = d_vect_temp; // swap important as not to allocate memory (everything on device side)
        cudaMemset(d_vect_out, 0.0, (no_nodes) * sizeof(float));

        update_instance_inputs<<<dimGrid,dimBlock>>>(blocks_nodes, d_vect_in, now.getYUpper(), now.getYBottom(), now.getX(), x, d_y, no_instances);

        dim3 dimGrid2(ceil((float)(no_nodes)/BLOCK_SIZE),1,1);
        dim3 dimBlock2(BLOCK_SIZE,1,1);
        SparseMUL<<<dimGrid2, dimBlock2>>>(column_idx, row_pointers, weights, d_vect_in, no_nodes, d_vect_out, no_nodes);

        // #### update step ####
        update_step<<<dimGrid,dimBlock>>>(d_vel_y, d_y, d_vect_out, blocks_nodes, no_instances);

        prev.update();
        now.update();
        next.update();
    }
    cudaFree(d_y);
    cudaFree(d_vel_y);
    cudaFree(d_colision);
    cudaFree(d_vect_in);
    cudaFree(d_vect_out);
}

void symulationCUDA_test(){
    srand(time(0));
    int *in;
    int *out;
    float *w;
    bool *enabled;
    int *innov;
    
    int no_instances;
    int *blocks_edges;
    int *translation;
    int *blocks_nodes;

    FILE *plik = fopen("crosoverCOO_test.txt", "r");
    if (plik == NULL) {
        return;
    }
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
    // end init
    // init to device
    int *d_in;
    int *d_out;
    float *d_w;
    bool *d_enabled;
    int *d_innov;
    
    int *d_blocks_edges;
    int *d_translation;
    int *d_blocks_nodes;
    // alocation
    cudaMalloc(&d_in, (blocks_edges[no_instances]) * sizeof(int));
    cudaMalloc(&d_out, (blocks_edges[no_instances]) * sizeof(int));
    cudaMalloc(&d_w, (blocks_edges[no_instances]) * sizeof(float));
    cudaMalloc(&d_enabled, (blocks_edges[no_instances]) * sizeof(bool));
    cudaMalloc(&d_innov, (blocks_edges[no_instances]) * sizeof(int));

    cudaMalloc(&d_blocks_edges, (no_instances+1) * sizeof(int));
    cudaMalloc(&d_translation, (blocks_nodes[no_instances]) * sizeof(int));
    cudaMalloc(&d_blocks_nodes, (no_instances+1) * sizeof(int));

    // data copy
    cudaMemcpy(d_in, in, (blocks_edges[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_out, out, (blocks_edges[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_w, w, (blocks_edges[no_instances]) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_enabled, enabled, (blocks_edges[no_instances]) * sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(d_innov, innov, (blocks_edges[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_blocks_edges, blocks_edges, (no_instances+1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_translation, translation, (blocks_nodes[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_blocks_nodes, blocks_nodes, (no_instances+1) * sizeof(int), cudaMemcpyHostToDevice);

    int* d_rewards;
    cudaMalloc(&d_rewards, (no_instances) * sizeof(int)); // vector of rewards

    int no_nodes = blocks_nodes[no_instances];
    int no_edges = blocks_edges[no_instances];
    // end init to device

    // ##### PREPARATION TO SYMULATION #####
    // TO CSR MATRIX
    int colsize = count_col_size(d_enabled, no_edges); // już zrównoleglone za pomocą redukcji

    int *d_col_idx; // size w
    float *d_weights; // size w
    int *d_row_pointers; //size translation+1
    
    cudaMalloc(&d_col_idx, (colsize) * sizeof(int));
    cudaMalloc(&d_weights, (colsize) * sizeof(float));
    cudaMalloc(&d_row_pointers, (no_nodes + 1) * sizeof(int));

    cudaMemset(d_row_pointers, 0, (no_nodes + 1) * sizeof(int));

    dim3 dimGrid(ceil((float)(no_instances)/BLOCK_SIZE),1,1);
    dim3 dimBlock(BLOCK_SIZE,1,1);

    updateRowPointers<<<dimGrid,dimBlock>>>(d_row_pointers, d_out, d_blocks_edges, d_blocks_nodes, no_instances, d_enabled);

    

    int *d_row_pointers_t;
    cudaMalloc(&d_row_pointers_t, (no_nodes + 1) * sizeof(int));

    cumulatedHistogram_with_copy(d_row_pointers, d_row_pointers_t, d_row_pointers, (no_nodes + 1));

    

    updateCol_idx_weights<<<dimGrid,dimBlock>>>(d_row_pointers_t, d_weights, d_col_idx, d_in, d_w, d_out, d_blocks_edges, d_blocks_nodes, no_instances, d_enabled);
    cudaFree(d_row_pointers_t); // important free
    
    // int* h_array = (int*)malloc((colsize) * sizeof(int));
    // cudaMemcpy(h_array, d_col_idx, (colsize) * sizeof(int), cudaMemcpyDeviceToHost);
    // for(int i=0; i<(colsize); i++) printf("%d\t", h_array[i]);
    // printf("\n");

    cudaMemset(d_rewards, -1, (no_instances) * sizeof(int));//set to -1

    
    // symulate
    symulate(d_blocks_nodes, d_col_idx, d_row_pointers, d_weights, d_rewards, no_nodes, no_instances, 100, no_nodes, no_edges); // max iteration to 10
    printf("Ended succesfuly!\n");
    printf("rewards: ");
    int* rewards;
    rewards = (int*) malloc((no_instances) * sizeof(int));

    cudaMemcpy(rewards, d_rewards, (no_instances) * sizeof(int), cudaMemcpyDeviceToHost);

    for(int i=0; i<no_instances; i++) printf("%d\t", rewards[i]);
    printf("\n");
    //set to -1
    // cudaMemset(d_rewards, -1, (no_instances) * sizeof(int));//set to -1
    // // symulate 2
    // symulate(d_blocks_nodes, d_col_idx, d_row_pointers, d_weights, d_rewards, no_nodes, no_instances, 100, no_nodes, no_edges); // max iteration to 100
    // printf("Ended succesfuly!\n");

}

int main(){
    symulationCUDA_test();
}
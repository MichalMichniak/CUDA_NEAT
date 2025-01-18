#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <curand_kernel.h>

#define WIDTH 800
#define HEIGHT 600
#define GAP_R 75
#define FLOPPY_RADIUS 10

#define MAX_ITERATION_GOAL 100000

#define POPULATION_COUNT 1000
#define SURVIVORS_TOURNAMENTS 600
#define K 200
#define BLOCK_SIZE 128


#define MUTATION_T1 12
#define MUTATION_T2 35
#define MUTATION_W_P 0.05
int next_node_innov = 8;
int next_edge_innov = 26;

#define BLOCK_SIZE 128
// ###### CUMULATED HISTOGRAM BEGIN #######
#define SECTION_SIZE 512

__global__ void scanKernelX_with_copy(int *Y, int *Y_copy, int *S, int *X, int width)
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

__global__ void scanKernelS_with_copy(int *S, int width)
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

__global__ void updateYKernel_with_copy(int *Y, int *Y_copy, int *S, int widthY)
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
    
    scanKernelX_with_copy<<<dimGrid,dimBlock>>>(d_Y, d_Y_copy, d_S, d_X, width);

    dim3 dimGrid2(1,1,1);
    dim3 dimBlock2(SECTION_SIZE,1,1);

    scanKernelS_with_copy<<<dimGrid2,dimBlock2>>>(d_S, ceil((float)width/SECTION_SIZE));
    
    updateYKernel_with_copy<<<dimGrid,dimBlock>>>(d_Y, d_Y_copy, d_S, width);

    cudaFree(d_S);
}
// ###### CUMULATED HISTOGRAM END #######

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
// ###### Reduction MAX BEGIN #######

int reductionSequential_MAX(int *input, int width)
{
    int maximum = -1;
    for (int i = 0; i < width; ++i)
    {
        maximum = max(maximum,input[i]);
    }

    return maximum;
}

__global__ void reductionKernelOp_MAX(int *input, int *output, int width)
{
    //@@ INSERT CODE HERE
    __shared__ int in[BLOCK_SIZE];
    if(threadIdx.x + (blockDim.x * blockIdx.x) < width){
        in[threadIdx.x] = (input[threadIdx.x + (blockDim.x * blockIdx.x)]);
    }else{
        in[threadIdx.x] = -1;
    }
    __syncthreads();
    int targetlevel = 0;
    int index = BLOCK_SIZE;
    while (index >>= 1) ++targetlevel;

    int result = -1;
    for(int it = 0; it < targetlevel; it++){
        if(threadIdx.x<(BLOCK_SIZE/(2<<it))){
            result = (int)max((float)in[threadIdx.x],(float)in[threadIdx.x+(BLOCK_SIZE/(2<<it))]);
            // printf("it: %d\t(%d, %d)\t%d\t%d\n",it,threadIdx.x, threadIdx.x+(1<<it) ,in[threadIdx.x] ,in[threadIdx.x+(1<<it)]);
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

int MAX_cryterion(int* d_rewards, int length)
{
    //@@ INSERT CODE HERE
    int *d_output, *h_output;
    cudaMalloc(&d_output, ceil((float)length/BLOCK_SIZE) * sizeof(int));
    h_output = (int*)malloc(ceil((float)length/BLOCK_SIZE) * sizeof(int));

    dim3 dimGrid(ceil((float)length/BLOCK_SIZE),1,1);
    dim3 dimBlock(BLOCK_SIZE,1,1);
    
    reductionKernelOp_MAX<<<dimGrid,dimBlock>>>(d_rewards,d_output,length);
    
    cudaMemcpy(h_output, d_output,  ceil((float)length/BLOCK_SIZE)*sizeof(int), cudaMemcpyDeviceToHost);

    int result = reductionSequential_MAX(h_output,ceil((float)length/BLOCK_SIZE));
    cudaFree(d_output);
    
    free(h_output);
    return result;
}
// ###### Reduction MAX END #######


// ###### Reduction MAX IDX BEGIN #######

int reductionSequential_MAX_IDX(int *input, int *input_idx, int width)
{
    int maximum = -1;
    int idx = -1;
    for (int i = 0; i < width; ++i)
    {   
        if(maximum<input[i]){
            maximum = input[i];
            idx = input_idx[i];
        }
    }
    printf("Maximum: %d\n", maximum);
    return idx;
}

__global__ void reductionKernelOp_MAX_IDX(int *input, int *output, int *output_idx, int width)
{
    //@@ INSERT CODE HERE
    __shared__ int in[BLOCK_SIZE];
    __shared__ int in_idx[BLOCK_SIZE];
    if(threadIdx.x + (blockDim.x * blockIdx.x) < width){
        in[threadIdx.x] = (input[threadIdx.x + (blockDim.x * blockIdx.x)]);
        in_idx[threadIdx.x] = threadIdx.x + (blockDim.x * blockIdx.x);
    }else{
        in[threadIdx.x] = -1;
        in_idx[threadIdx.x] = -1;
    }
    __syncthreads();
    int targetlevel = 0;
    int index = BLOCK_SIZE;
    while (index >>= 1) ++targetlevel;

    int result = -1;
    int result_idx = -1;
    for(int it = 0; it < targetlevel; it++){
        if(threadIdx.x<(BLOCK_SIZE/(2<<it))){
            if(in[threadIdx.x]>in[threadIdx.x+(BLOCK_SIZE/(2<<it))]){
                result = in[threadIdx.x];
                result_idx = in_idx[threadIdx.x];
            }else{
                result = in[threadIdx.x+(BLOCK_SIZE/(2<<it))];
                result_idx = in_idx[threadIdx.x+(BLOCK_SIZE/(2<<it))];
            }
            
            // printf("it: %d\t(%d, %d)\t%d\t%d\n",it,threadIdx.x, threadIdx.x+(1<<it) ,in[threadIdx.x] ,in[threadIdx.x+(1<<it)]);
        }
        __syncthreads();
        if((threadIdx.x<(BLOCK_SIZE/(2<<it)))){
            in[threadIdx.x] = result;
            in_idx[threadIdx.x] = result_idx;
        }
        //printf("%f\n",in[0]);
        __syncthreads();
    }
    
    if(threadIdx.x == 0){
        output[blockIdx.x] = in[0]; 
        output_idx[blockIdx.x] = in_idx[0];
        // printf("%f",in[0]);
        // printf("\n%d\n",targetlevel);
    }
    
}

int MAX_cryterion_IDX(int* d_rewards, int length)
{
    //@@ INSERT CODE HERE
    int *d_output, *h_output,*d_output_idx, *h_output_idx;
    cudaMalloc(&d_output, ceil((float)length/BLOCK_SIZE) * sizeof(int));
    cudaMalloc(&d_output_idx, ceil((float)length/BLOCK_SIZE) * sizeof(int));
    h_output = (int*)malloc(ceil((float)length/BLOCK_SIZE) * sizeof(int));
    h_output_idx = (int*)malloc(ceil((float)length/BLOCK_SIZE) * sizeof(int));

    dim3 dimGrid(ceil((float)length/BLOCK_SIZE),1,1);
    dim3 dimBlock(BLOCK_SIZE,1,1);
    
    reductionKernelOp_MAX_IDX<<<dimGrid,dimBlock>>>(d_rewards,d_output,d_output_idx,length);
    
    cudaMemcpy(h_output, d_output,  ceil((float)length/BLOCK_SIZE)*sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_output_idx, d_output_idx,  ceil((float)length/BLOCK_SIZE)*sizeof(int), cudaMemcpyDeviceToHost);

    int result = reductionSequential_MAX_IDX(h_output,h_output_idx,ceil((float)length/BLOCK_SIZE));
    cudaFree(d_output);
    cudaFree(d_output_idx);
    
    free(h_output);
    free(h_output_idx);
    return result;
}
// ###### Reduction MAX IDX END #######

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

void symulate(
    int *d_in,
    int *d_out,
    float *d_w,
    bool *d_enabled,
    int *d_innov,
    int *d_blocks_edges,
    int *d_translation,
    int *d_blocks_nodes,
    int* d_rewards,
    int no_instances,
    int no_nodes,
    int no_edges,
    int max_it
){
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
    symulate(d_blocks_nodes, d_col_idx, d_row_pointers, d_weights, d_rewards, no_nodes, no_instances, max_it, no_nodes, no_edges); // max iteration to 10
}

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
        // printf("%d\t%d\t%d\t%d\n", instance_number, 1+i, blocks_edges[instance_number+1]-blocks_edges[instance_number], instance_number);
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
    int *translation_t1,
    int *translation_t2
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

                translation_t1[idx_first - blocks_nodes[first] + new_blocks_nodes[i+offset]] = idx;
                translation_t2[idx_second - blocks_nodes[second] + new_blocks_nodes[i+offset]] = idx;
                idx++;
                idx_first++;
                idx_second++;
            }else if(translation[idx_first] < translation[idx_second]){
                new_translation[idx + new_blocks_nodes[i+offset]] = translation[idx_first];
                translation_t1[idx_first - blocks_nodes[first] + new_blocks_nodes[i+offset]] = idx;
                idx++;
                idx_first++;
            }else{
                new_translation[idx + new_blocks_nodes[i+offset]] = translation[idx_second];
                translation_t2[idx_second - blocks_nodes[second] + new_blocks_nodes[i+offset]] = idx;
                idx++;
                idx_second++;
            }
        }
        while(idx_first != blocks_nodes[first+1]){
            new_translation[idx + new_blocks_nodes[i+offset]] = translation[idx_first];
            translation_t1[idx_first - blocks_nodes[first] + new_blocks_nodes[i+offset]] = idx;
            idx++;
            idx_first++;
        }
        while(idx_second != blocks_nodes[second+1]){
            new_translation[idx + new_blocks_nodes[i+offset]] = translation[idx_second];
            translation_t2[idx_second - blocks_nodes[second] + new_blocks_nodes[i+offset]] = idx;
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
                new_in[idx+new_blocks_edges[i+offset]] = translation_t1[in[idx_first]+ new_blocks_nodes[i+offset]];
                new_out[idx+new_blocks_edges[i+offset]] = translation_t1[out[idx_first]+ new_blocks_nodes[i+offset]];
                if(new_in[idx+new_blocks_edges[i+offset]]<0) printf("KRZYZ1\t%d\t%d\t%d\n",idx+new_blocks_edges[i+offset], i+offset, translation_t1[in[idx_first]+ new_blocks_nodes[i+offset]]);
                idx_first++;
                idx_second++;
                idx++;
            }else if(innov[idx_first] < innov[idx_second]){
                // if(i==0) printf("%d, %d\t", in[idx_first]+blocks_nodes[first_pair[i]],translation_t[in[idx_first]+blocks_nodes[first_pair[i]]]);
                new_w[idx+new_blocks_edges[i+offset]] = w[idx_first];
                new_enabled[idx+new_blocks_edges[i+offset]] = enabled[idx_first];
                new_innov[idx+new_blocks_edges[i+offset]] = innov[idx_first];
                new_in[idx+new_blocks_edges[i+offset]] = translation_t1[in[idx_first]+ new_blocks_nodes[i+offset]];
                if(new_in[idx+new_blocks_edges[i+offset]]<0) printf("KRZYZ1\t%d\t%d\t%d\n",idx+new_blocks_edges[i+offset], i+offset, translation_t1[in[idx_first]+ new_blocks_nodes[i+offset]]);
                new_out[idx+new_blocks_edges[i+offset]] = translation_t1[out[idx_first]+ new_blocks_nodes[i+offset]];
                idx_first++;
                idx++;
            }else{
                // if(i==0) printf("%d, %d\t", in[idx_second]+blocks_nodes[second_pair[i]],translation_t[in[idx_second]+blocks_nodes[second_pair[i]]]);
                new_w[idx+new_blocks_edges[i+offset]] = w[idx_second];
                new_enabled[idx+new_blocks_edges[i+offset]] = enabled[idx_second];
                new_innov[idx+new_blocks_edges[i+offset]] = innov[idx_second];
                new_in[idx+new_blocks_edges[i+offset]] = translation_t2[in[idx_second]+ new_blocks_nodes[i+offset]];
                if(new_in[idx+new_blocks_edges[i+offset]]<0) printf("KRZYZ1\t%d\t%d\t%d\n",idx+new_blocks_edges[i+offset], i+offset, translation_t2[in[idx_second]+ new_blocks_nodes[i+offset]]);
                new_out[idx+new_blocks_edges[i+offset]] = translation_t2[out[idx_second]+ new_blocks_nodes[i+offset]];
                idx_second++;
                idx++;
            }

        }
        while(idx_first != blocks_edges[first+1]){
            // if(i==0) printf("%d, %d\t", in[idx_first]+blocks_nodes[first_pair[i]],translation_t[in[idx_first]+blocks_nodes[first_pair[i]]]);
            new_w[idx+new_blocks_edges[i+offset]] = w[idx_first];
            new_enabled[idx+new_blocks_edges[i+offset]] = enabled[idx_first];
            new_innov[idx+new_blocks_edges[i+offset]] = innov[idx_first];
            new_in[idx+new_blocks_edges[i+offset]] = translation_t1[in[idx_first]+ new_blocks_nodes[i+offset]];
            if(new_in[idx+new_blocks_edges[i+offset]]<0) printf("KRZYZ1\t%d\t%d\t%d\n",idx+new_blocks_edges[i+offset], i+offset, translation_t1[in[idx_first]+ new_blocks_nodes[i+offset]]);
            new_out[idx+new_blocks_edges[i+offset]] = translation_t1[out[idx_first]+ new_blocks_nodes[i+offset]];
            idx_first++;
            idx++;
        }
        while(idx_second != blocks_edges[second+1]){
            // if(i==0) printf("%d, %d\t", in[idx_second]+blocks_nodes[second_pair[i]],translation_t[in[idx_second]+blocks_nodes[second_pair[i]]]);
            new_w[idx+new_blocks_edges[i+offset]] = w[idx_second];
            new_enabled[idx+new_blocks_edges[i+offset]] = enabled[idx_second];
            new_innov[idx+new_blocks_edges[i+offset]] = innov[idx_second];
            new_in[idx+new_blocks_edges[i+offset]] = translation_t2[in[idx_second]+ new_blocks_nodes[i+offset]];
            if(new_in[idx+new_blocks_edges[i+offset]]<0) printf("KRZYZ1\t%d\t%d\t%d\n",idx+new_blocks_edges[i+offset], i+offset, translation_t2[in[idx_second]+ new_blocks_nodes[i+offset]]);
            new_out[idx+new_blocks_edges[i+offset]] = translation_t2[out[idx_second]+ new_blocks_nodes[i+offset]];
            idx_second++;
            idx++;
        }
        // if(i==0) printf("\nEND\n");
        i += blockDim.x * gridDim.x;
    }
    
}
__global__ void initialize_rng(curandState* state, unsigned long seed) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(clock64()+seed, idx, 0, &state[idx]);
}

__global__ void selection_step(curandState* state,int *rewards, int *mask, int no_instances){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<SURVIVORS_TOURNAMENTS){
        int best_idx = 0;
        int best_reward = -1;
        for(int k=0; k<K; k++){
            int idx = (int)((1.0-curand_uniform(&state[i])) * no_instances);
            if(rewards[idx]>best_reward){
                best_reward = rewards[idx];
                best_idx = idx;
            }
        }
        mask[best_idx] = 1;
        i += blockDim.x * gridDim.x;
    }
}

__global__ void selection_crossover_step(curandState* state,int *rewards, int *first_parent, int *second_parent, int no_offspring, int no_instances){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<no_offspring){
        int best_idx = 0;
        int best_reward = -1;
        for(int k=0; k<K; k++){
            int idx = (int)((1.0-curand_uniform(&state[i])) * no_instances);
            if(rewards[idx]>best_reward){
                best_reward = rewards[idx];
                best_idx = idx;
            }
        }
        first_parent[i] = best_idx;
        best_idx = 0;
        best_reward = -1;
        for(int k=0; k<K; k++){
            int idx = (int)((1.0-curand_uniform(&state[i])) * no_instances);
            if(rewards[idx]>best_reward){
                best_reward = rewards[idx];
                best_idx = idx;
            }
        }
        second_parent[i] = best_idx;
        i += blockDim.x * gridDim.x;
    }
}

__global__ void countMutations(curandState* state, int* mutation_parent, int* rewards ,int *blocks_edges, int* blocks_nodes, int *new_blocks_edges, int* new_blocks_nodes, int no_instances, int offset){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<MUTATION_T1 + MUTATION_T2){
        int best_idx = 0;
        int best_reward = -1;
        for(int k=0; k<K; k++){
            int idx = (int)((1.0-curand_uniform(&state[i])) * no_instances);
            if(rewards[idx]>best_reward){
                best_reward = rewards[idx];
                best_idx = idx;
            }
        }
        if(i<MUTATION_T1){
            // printf("%d,%d\n", offset,i+offset);
            new_blocks_edges[i+offset+1] = blocks_edges[best_idx+1] - blocks_edges[best_idx] + 2;
            new_blocks_nodes[i+offset+1] = blocks_nodes[best_idx+1] - blocks_nodes[best_idx] + 1;
            mutation_parent[i] = best_idx;
        }else{
            // printf("%d,%d\n", offset,i+offset);
            new_blocks_edges[i+offset+1] = blocks_edges[best_idx+1] - blocks_edges[best_idx] + 1;
            new_blocks_nodes[i+offset+1] = blocks_nodes[best_idx+1] - blocks_nodes[best_idx];
            mutation_parent[i] = best_idx;
        }

        i += blockDim.x * gridDim.x;
    }
}

__global__ void MutateT1T2(
    curandState* state,
    int* mutation_parent,
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
    int *new_in,
    int *new_out,
    float *new_w,
    bool *new_enabled,
    int *new_innov,
    int *new_translation,
    int *translation,
    int offset_edges_innov,
    int offset_nodes_innov
    ){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<MUTATION_T1 + MUTATION_T2){
        if(i<MUTATION_T1){
            // chose random edge:
            int length_edges = blocks_edges[mutation_parent[i]+1] - blocks_edges[mutation_parent[i]];
            int split_idx = (int)((1.0-curand_uniform(&state[i])) * length_edges);

            for(int j=0; j<new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 2; j++){
                new_in[j+new_blocks_edges[offset+i]] = in[j + blocks_edges[mutation_parent[i]]];
                if(in[j + blocks_edges[mutation_parent[i]]]<0) printf("MUTACJA\t");
                new_out[j+new_blocks_edges[offset+i]] = out[j + blocks_edges[mutation_parent[i]]];
                if(split_idx != j){
                    new_enabled[j+new_blocks_edges[offset+i]] = enabled[j + blocks_edges[mutation_parent[i]]];
                }else{
                    new_enabled[j+new_blocks_edges[offset+i]] = false;
                }
                new_innov[j+new_blocks_edges[offset+i]] = innov[j + blocks_edges[mutation_parent[i]]];
                new_w[j+new_blocks_edges[offset+i]] = w[j + blocks_edges[mutation_parent[i]]];
            }
            // new edges
            new_in[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 2 +new_blocks_edges[offset+i]] = in[split_idx+blocks_edges[mutation_parent[i]]];
            new_out[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 2 +new_blocks_edges[offset+i]] = new_blocks_nodes[offset+i+1] - new_blocks_nodes[offset+i] - 1;
            new_enabled[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 2 +new_blocks_edges[offset+i]] = true;
            new_innov[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 2 +new_blocks_edges[offset+i]] = offset_edges_innov + 2*i;
            new_w[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 2 +new_blocks_edges[offset+i]] = w[split_idx+blocks_edges[mutation_parent[i]]];
            if(in[split_idx+blocks_edges[mutation_parent[i]]]<0) printf("MUTACJA\t");
            new_in[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = new_blocks_nodes[offset+i+1] - new_blocks_nodes[offset+i] - 1;
            new_out[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = out[split_idx+blocks_edges[mutation_parent[i]]];
            new_enabled[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = true;
            new_innov[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = offset_edges_innov + 2*i + 1;
            new_w[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = w[split_idx+blocks_edges[mutation_parent[i]]];

            for(int j=0; j<new_blocks_nodes[offset+i+1] - new_blocks_nodes[offset+i] - 1; j++){
                // printf("\n%d, %d\n", j+new_blocks_nodes[offset+i], translation[j + blocks_nodes[mutation_parent[i]]]);
                new_translation[j+new_blocks_nodes[offset+i]] = translation[j + blocks_nodes[mutation_parent[i]]];
            }
            // new nodes
            new_translation[new_blocks_nodes[offset+i] + new_blocks_nodes[offset+i+1] - new_blocks_nodes[offset+i] - 1] = offset_nodes_innov + i;

        }else{
            int length_nodes = new_blocks_nodes[offset+i+1] - new_blocks_nodes[offset+i];
            int first_idx = (int)((1.0-curand_uniform(&state[i])) * length_nodes);
            int second_idx = (int)((1.0-curand_uniform(&state[i])) * length_nodes);

            for(int j=0; j<new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1; j++){
                // printf("%d\n", j+new_blocks_edges[offset+i+1]);
                new_in[j+new_blocks_edges[offset+i]] = in[j + blocks_edges[mutation_parent[i]]];
                new_out[j+new_blocks_edges[offset+i]] = out[j + blocks_edges[mutation_parent[i]]];
                new_enabled[j+new_blocks_edges[offset+i]] = enabled[j + blocks_edges[mutation_parent[i]]];
                new_innov[j+new_blocks_edges[offset+i]] = innov[j + blocks_edges[mutation_parent[i]]];
                new_w[j+new_blocks_edges[offset+i]] = w[j + blocks_edges[mutation_parent[i]]];
            }

            new_in[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = first_idx;
            new_out[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = second_idx;
            new_enabled[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = true;
            new_innov[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = offset_edges_innov + 2*MUTATION_T1 + (i - MUTATION_T1);
            new_w[new_blocks_edges[offset+i+1] - new_blocks_edges[offset+i] - 1 +new_blocks_edges[offset+i]] = (float)((1.0-curand_uniform(&state[i])) - 0.5);
            if(first_idx<0) printf("MUTACJA\t");
            for(int j=0; j<new_blocks_nodes[offset+i+1] - new_blocks_nodes[offset+i]; j++){
                new_translation[j+new_blocks_nodes[offset+i]] = translation[j + blocks_nodes[mutation_parent[i]]];
            }
        }
        i += blockDim.x * gridDim.x;
    }
}

__global__ void mutate_weights(curandState* state, float *w, int no_edges){
    int i = threadIdx.x + (blockDim.x * blockIdx.x);
    while(i<no_edges){
        if(curand_uniform(&state[i])<MUTATION_W_P){
            w[i] += curand_uniform(&state[i]) - 0.5;
        }
        i += blockDim.x * gridDim.x;
    }
}


void get_new_population(
    int **d_in_,
    int **d_out_,
    float **d_w_,
    bool **d_enabled_,
    int **d_innov_,
    int **d_blocks_edges_,
    int **d_translation_,
    int **d_blocks_nodes_,
    int *d_rewards,
    int no_instances,
    int *no_nodes_,
    int *no_edges_
    ){
    // przepisanie
    int *d_in;
    int *d_out;
    float *d_w;
    bool *d_enabled;
    int *d_innov;
    
    int *d_blocks_edges;
    int *d_translation;
    int *d_blocks_nodes;

    d_in = *d_in_;
    d_out = *d_out_;
    d_w = *d_w_;
    d_enabled = *d_enabled_;
    d_innov = *d_innov_;

    d_blocks_edges = *d_blocks_edges_;
    d_translation = *d_translation_;
    d_blocks_nodes = *d_blocks_nodes_;
    //


    int no_survivors; // number of 1 in mask
    int no_mutations; // TODO: mutacia

    int no_offsprings; // ze starej populacji

    int *d_mask;
    int *d_first_pair;
    int *d_second_pair;

    curandState *d_state;
    int NGrid = ceil((float)(SURVIVORS_TOURNAMENTS)/BLOCK_SIZE);

    cudaMalloc(&d_mask, no_instances * sizeof(int));
    cudaMalloc(&d_state, NGrid * BLOCK_SIZE * sizeof(curandState));

    initialize_rng<<<NGrid, BLOCK_SIZE>>>(d_state, time(0));
    selection_step<<<NGrid, BLOCK_SIZE>>>(d_state, d_rewards, d_mask, no_instances);
    cudaFree(d_state);
    
    

    int *h_mask;
    h_mask = (int*) malloc(no_instances * sizeof(int));
    cudaMemcpy(h_mask, d_mask, no_instances * sizeof(int), cudaMemcpyDeviceToHost);
    int idx = 0;
    for(int i=0; i<no_instances; i++){ // to będzie sekwencyjne bo tylko przepisanie
        if(h_mask[i] == 1){
            idx++;
        }
    }
    no_survivors = idx;
    
    int *istance_numbers_seq; // tablica numeru instancji mapowanej ze starej do nowej tablicy szie: no_survivors
    istance_numbers_seq = (int*)malloc(no_survivors*sizeof(int));
    idx = 0;
    for(int i=0; i<no_instances; i++){ // to będzie sekwencyjne bo tylko przepisanie
        if(h_mask[i] == 1){
            istance_numbers_seq[idx] = i;
            idx++;
        }
    }

    no_offsprings = POPULATION_COUNT-no_survivors - MUTATION_T1 - MUTATION_T2;
    no_mutations = MUTATION_T1 + MUTATION_T2;

    cudaMalloc(&d_first_pair, no_offsprings * sizeof(int));
    cudaMalloc(&d_second_pair, no_offsprings * sizeof(int));
    NGrid = ceil((float)(no_offsprings)/BLOCK_SIZE);
    cudaMalloc(&d_state, NGrid * BLOCK_SIZE * sizeof(curandState));
    initialize_rng<<<NGrid, BLOCK_SIZE>>>(d_state, time(0));
    selection_crossover_step<<<NGrid, BLOCK_SIZE>>>(d_state,d_rewards, d_first_pair, d_second_pair, no_offsprings, no_instances);
    cudaFree(d_state);

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
    int *d_mutation_parent;
    d_mutation_parent = (int*)malloc((no_mutations)*sizeof(int));
    cudaMalloc(&d_mutation_parent, no_mutations * sizeof(int));
    NGrid = ceil((float)(no_mutations)/BLOCK_SIZE);
    cudaMalloc(&d_state, NGrid * BLOCK_SIZE * sizeof(curandState));
    initialize_rng<<<NGrid, BLOCK_SIZE>>>(d_state, time(0));
    countMutations<<<NGrid, BLOCK_SIZE>>>(d_state, d_mutation_parent, d_rewards , d_blocks_edges, d_blocks_nodes, d_new_blocks_edges, d_new_blocks_nodes, no_instances, no_offsprings+no_survivors);


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
    int *d_translation_t1;
    
    cudaMalloc(&d_translation_t1, no_nodes * sizeof(int));

    int *d_translation_t2;
    

    cudaMalloc(&d_translation_t2, no_nodes * sizeof(int));
    // cudaMemset(d_translation_t, 0, (no_nodes) * sizeof(int));
    // printf("%d\n", no_nodes);
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
    d_translation_t1,
    d_translation_t2
    );

    MutateT1T2<<<NGrid, BLOCK_SIZE>>>(
    d_state,
    d_mutation_parent,
    d_blocks_nodes,
    d_blocks_edges,
    d_new_blocks_nodes,
    d_new_blocks_edges,
    no_offsprings+no_survivors, // bardzo ważny offset
    d_in,
    d_out,
    d_w,
    d_enabled,
    d_innov,
    d_new_in,
    d_new_out,
    d_new_w,
    d_new_enabled,
    d_new_innov,
    d_new_translation,
    d_translation,
    next_edge_innov,
    next_node_innov
    );

    next_node_innov += MUTATION_T1;
    next_edge_innov += MUTATION_T1*2 + MUTATION_T2;

    cudaFree(d_state);
    NGrid = ceil((float)(no_edges)/BLOCK_SIZE);
    cudaMalloc(&d_state, NGrid * BLOCK_SIZE * sizeof(curandState));
    initialize_rng<<<NGrid, BLOCK_SIZE>>>(d_state, time(0));
    mutate_weights<<<NGrid, BLOCK_SIZE>>>(d_state, d_new_w, no_edges);
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
    *d_in_=d_new_in;
    *d_out_=d_new_out;
    *d_w_=d_new_w;
    *d_enabled_=d_new_enabled;
    *d_innov_=d_new_innov;

    *d_blocks_edges_=d_new_blocks_edges;
    *d_translation_=d_new_translation;
    *d_blocks_nodes_=d_new_blocks_nodes;
    *no_edges_ = no_edges;
    *no_nodes_ = no_nodes;
    // czyszczenie
    cudaFree(d_state);
    cudaFree(d_translation_t1);
    cudaFree(d_translation_t2);
    cudaFree(d_istance_numbers_seq);
    free(istance_numbers_seq);
    free(h_mask);

    // end czyszczenie
}

void main_loop(){
    // initial population file
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

    // init population:
    int *in;
    int *out;
    float *w;
    bool *enabled;
    int *innov; 

    int no_instances;
    int *blocks_edges;
    int *translation; // zał że w instancjach posortowane rosnąco
    int *blocks_nodes;
    
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
    fclose(plik);

    // ###### inicializacja i alokacja na Device ######
    int *d_in;
    int *d_out;
    float *d_w;
    bool *d_enabled;
    int *d_innov;
    
    int *d_blocks_edges;
    int *d_translation;
    int *d_blocks_nodes;

    int *d_rewards_init;
    // alocation
    cudaMalloc(&d_in, (blocks_edges[no_instances]) * sizeof(int));
    cudaMalloc(&d_out, (blocks_edges[no_instances]) * sizeof(int));
    cudaMalloc(&d_w, (blocks_edges[no_instances]) * sizeof(float));
    cudaMalloc(&d_enabled, (blocks_edges[no_instances]) * sizeof(bool));
    cudaMalloc(&d_innov, (blocks_edges[no_instances]) * sizeof(int));

    cudaMalloc(&d_blocks_edges, (no_instances+1) * sizeof(int));
    cudaMalloc(&d_translation, (blocks_nodes[no_instances]) * sizeof(int));
    cudaMalloc(&d_blocks_nodes, (no_instances+1) * sizeof(int));

    cudaMalloc(&d_rewards_init, (no_instances) * sizeof(int));

    // data copy
    cudaMemcpy(d_in, in, (blocks_edges[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_out, out, (blocks_edges[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_w, w, (blocks_edges[no_instances]) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_enabled, enabled, (blocks_edges[no_instances]) * sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(d_innov, innov, (blocks_edges[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(d_blocks_edges, blocks_edges, (no_instances+1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_translation, translation, (blocks_nodes[no_instances]) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_blocks_nodes, blocks_nodes, (no_instances+1) * sizeof(int), cudaMemcpyHostToDevice);

    cudaMemset(d_rewards_init, 0, (no_instances) * sizeof(int));
    
    //###### inicializacja i alokacja END ######
    // ##### powielenie populacji #####
    int no_edges;
    int no_nodes;
    get_new_population(
    &d_in,
    &d_out,
    &d_w,
    &d_enabled,
    &d_innov,
    &d_blocks_edges,
    &d_translation,
    &d_blocks_nodes,
    d_rewards_init,
    no_instances,
    &no_nodes,
    &no_edges
    );


    // inicjalizacja zmiennych:
    // int *new_in;
    // int *new_out;
    // float *new_w;
    // bool *new_enabled;
    // int *new_innov;

    // int *new_translation; // 
    // int *new_blocks_nodes;
    // int new_no_instances = POPULATION_COUNT;

    // new_blocks_nodes = (int*)malloc((1+POPULATION_COUNT)*sizeof(int)); // 0 | 1 do no_survivors | no_survivors + 1 do no_survivors + no_offsprings | no_survivors + no_offsprings + 1 do no_survivors + no_offsprings + no_mutation 
     
    int *d_rewards;
    cudaMalloc(&d_rewards, (POPULATION_COUNT) * sizeof(int));

    for(int i=0; i<100; i++){ // here number of iterations
        cudaMemset(d_rewards, 0, (POPULATION_COUNT) * sizeof(int));
        symulate(
        d_in,
        d_out,
        d_w,
        d_enabled,
        d_innov,
        d_blocks_edges,
        d_translation,
        d_blocks_nodes,
        d_rewards,
        POPULATION_COUNT,
        no_nodes,
        no_edges,
        MAX_ITERATION_GOAL
        );

        int best_reward;
        best_reward = MAX_cryterion(d_rewards,POPULATION_COUNT);
        printf("maximum reward: %d\n", best_reward);
        if(best_reward == MAX_ITERATION_GOAL-1) break;
        get_new_population(
        &d_in,
        &d_out,
        &d_w,
        &d_enabled,
        &d_innov,
        &d_blocks_edges,
        &d_translation,
        &d_blocks_nodes,
        d_rewards,
        POPULATION_COUNT,
        &no_nodes,
        &no_edges
        );

        // new_in = (int*)malloc(no_edges*sizeof(int));
        // new_out = (int*)malloc(no_edges*sizeof(int));
        // new_w = (float*)malloc(no_edges*sizeof(float));
        // new_enabled = (bool*)malloc(no_edges*sizeof(bool));
        // new_innov = (int*)malloc(no_edges*sizeof(int));
        // new_translation = (int*)malloc(no_nodes*sizeof(int));

        // cudaMemcpy(new_in, d_in, no_edges * sizeof(int), cudaMemcpyDeviceToHost);
        // cudaMemcpy(new_out, d_out, no_edges * sizeof(int), cudaMemcpyDeviceToHost);
        // cudaMemcpy(new_w, d_w, no_edges * sizeof(float), cudaMemcpyDeviceToHost);
        // cudaMemcpy(new_enabled, d_enabled, no_edges * sizeof(bool), cudaMemcpyDeviceToHost);
        // cudaMemcpy(new_innov, d_innov, no_edges * sizeof(int), cudaMemcpyDeviceToHost);

        
        // cudaMemcpy(new_translation, d_translation, no_nodes * sizeof(int), cudaMemcpyDeviceToHost);
        // cudaMemcpy(new_blocks_nodes, d_blocks_nodes, (new_no_instances+1) * sizeof(int), cudaMemcpyDeviceToHost);
        // printf("\nnew block nodes: ");

        // for(int i=0; i<new_no_instances+1; i++){
        //     printf("%d\t", new_blocks_nodes[i]);
        // }
        // int *new_blocks_edges;
        // new_blocks_edges = (int*)malloc((1+POPULATION_COUNT)*sizeof(int));
        // cudaMemcpy(new_blocks_edges, d_blocks_edges, (new_no_instances+1) * sizeof(int), cudaMemcpyDeviceToHost);
        // printf("\nnew block edges: ");
        // for(int i=0; i<new_no_instances+1; i++){
        //     printf("%d\t", new_blocks_edges[i]);
        // }
        
        // printf("\nnew translation: ");
        // for(int i=0; i<no_nodes; i++){
        //     printf("%d\t", new_translation[i]);
        // }
        // printf("\nnew in: ");
        // for(int i=0; i<no_edges; i++){
        //     printf("%d\t", new_in[i]);
        // }

        // printf("\nnew out: ");
        // for(int i=0; i<no_edges; i++){
        //     printf("%d\t", new_out[i]);
        // }

        // printf("\nnew innov: ");
        // for(int i=0; i<no_edges; i++){
        //     printf("%d\t", new_innov[i]);
        // }
        // // printf("%d\n",i);
        // free(new_in);
        // free(new_out);
        // free(new_w);
        // free(new_enabled);
        // free(new_innov);
        // free(new_translation);
        printf("number of nodes: %d\tnumber of edges: %d\tmaximum reward: %d\n", no_nodes, no_edges, best_reward);
    }
    int copy_idx = MAX_cryterion_IDX(d_rewards, POPULATION_COUNT);
    printf("idx: %d\n", copy_idx);
    int first_node; 
    int last_node;
    int first_edge; 
    int last_edge;

    cudaMemcpy(&first_edge, d_blocks_edges+copy_idx, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&last_edge, d_blocks_edges+copy_idx+1, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&first_node, d_blocks_nodes+copy_idx, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&last_node, d_blocks_nodes+copy_idx+1, sizeof(int), cudaMemcpyDeviceToHost);

    int edge_size = last_edge - first_edge;
    int node_size = last_node - first_node;

    int *save_in;
    int *save_out;
    float *save_w;
    bool *save_enabled;

    save_in = (int*)malloc(edge_size*sizeof(int));
    save_out = (int*)malloc(edge_size*sizeof(int));
    save_w = (float*)malloc(edge_size*sizeof(float));
    save_enabled = (bool*)malloc(edge_size*sizeof(bool));


    cudaMemcpy(save_in, d_in+first_edge, edge_size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(save_out, d_out+first_edge, edge_size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(save_w, d_w+first_edge, edge_size * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(save_enabled, d_enabled+first_edge, edge_size * sizeof(bool), cudaMemcpyDeviceToHost);

    plik = fopen("best_net.txt", "w");
    if (plik != NULL) {
        fprintf(plik,"%d %d\n", node_size, edge_size);
        for (int i = 0; i < edge_size; i++) {
            fprintf(plik, "%d ", save_in[i]);
        }
        fprintf(plik,"\n");
        for (int i = 0; i < edge_size; i++) {
            fprintf(plik, "%d ", save_out[i]);
        }
        fprintf(plik,"\n");
        for (int i = 0; i < edge_size; i++) {
            fprintf(plik, "%f ", save_w[i]);
        }
        fprintf(plik,"\n");
        for (int i = 0; i < edge_size; i++) {
            fprintf(plik, "%d ", (int)save_enabled[i]);
        }

        fclose(plik);
    }
    
    


    free(save_in);
    free(save_out);
    free(save_w);
    free(save_enabled);

    // printf("\nnew block nodes: ");

    // for(int i=0; i<new_no_instances+1; i++){
    //     printf("%d\t", new_blocks_nodes[i]);
    // }
    // int *new_blocks_edges;
    // new_blocks_edges = (int*)malloc((1+POPULATION_COUNT)*sizeof(int));
    // cudaMemcpy(new_blocks_edges, d_blocks_edges, (new_no_instances+1) * sizeof(int), cudaMemcpyDeviceToHost);
    // printf("\nnew block edges: ");
    // for(int i=0; i<new_no_instances+1; i++){
    //     printf("%d\t", new_blocks_edges[i]);
    // }
    
    // printf("\nnew translation: ");
    // for(int i=0; i<no_nodes; i++){
    //     printf("%d\t", new_translation[i]);
    // }
    // printf("\nnew in: ");
    // for(int i=0; i<no_edges; i++){
    //     printf("%d\t", new_in[i]);
    // }

    // printf("\nnew out: ");
    // for(int i=0; i<no_edges; i++){
    //     printf("%d\t", new_out[i]);
    // }

    // printf("\nnew innov: ");
    // for(int i=0; i<no_edges; i++){
    //     printf("%d\t", new_innov[i]);
    // }

    // printf("\nnew enabled: ");
    // for(int i=0; i<no_edges; i++){
    //     printf("%d\t", (int)new_enabled[i]);
    // }
    // printf("\n");

    // tutaj free TODO:
    cudaFree(d_rewards);
    cudaFree(d_in);
    cudaFree(d_out);
    cudaFree(d_w);
    cudaFree(d_enabled);
    cudaFree(d_innov);
    cudaFree(d_blocks_edges);
    cudaFree(d_blocks_nodes);
    cudaFree(d_translation);
    cudaFree(d_rewards_init);

    free(blocks_nodes);
    free(blocks_edges);
    free(translation);
    free(innov);
    free(enabled);
    free(in);
    free(out);
    free(w);
}




int main(){
    main_loop();
}
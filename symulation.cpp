#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstdlib>
#include <cmath>

#define WIDTH 800
#define HEIGHT 600
#define GAP_R 140
#define FLOPPY_RADIUS 10

float sigmoid(float x){
    return 1/(1+exp(-x));
}

// Mnożenie macierzy CSR przez wektor 
void SparseMUL(int* column_idx, int* row_pointers, float* weights, float* input_vector, int vector_size, float* output_vector, int output_vector_size){
    for(int row=0; row<output_vector_size; row++){ // To do zrównoleglenia
        // co będzie w kernelu
        // kopia __shared__ części wejściowego wektora
        float acc = 0;
        for(int col_idx=row_pointers[row]; col_idx<row_pointers[row+1]; col_idx++){
            acc += input_vector[column_idx[col_idx]] * weights[col_idx];
        }
        output_vector[row] = sigmoid(acc);
    }
}

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
void update_colisions(bool* colision, int* rewards, float* y, float x, int y_upper_pipe, int y_bottom_pipe, int x_pipe, int r_pipe, int no_instances){
    for(int i=0; i<no_instances; i++){ // TODO: zrównoleglenie na i
        if(((x + FLOPPY_RADIUS >=x_pipe-r_pipe && x-FLOPPY_RADIUS<=x_pipe+r_pipe) && (y[i]-FLOPPY_RADIUS<=y_upper_pipe || y[i]+FLOPPY_RADIUS>=HEIGHT-y_bottom_pipe)) || (y[i]>=HEIGHT-FLOPPY_RADIUS || y[i]<=0+FLOPPY_RADIUS)){
            colision[i] = true;
        }
        else{
            if(!colision[i]){
                rewards[i]++;
            }
        }
        
    }
}

void update_step(float* vel_y, float* y, float* vect_out,int *blocks_nodes, int no_instances){
    for(int i=0; i<no_instances; i++){ // TODO: zrównoleglenie na i
        float g = 0.1;
        float flop_v = 4;
        vel_y[i] += g;
        if(vect_out[blocks_nodes[i]]>0.5){
            vel_y[i] = -flop_v;
        }
        y[i]+=vel_y[i];
    }
}

void update_instance_inputs(int *blocks_nodes, float* vect_in, int y_upper_pipe, int y_bottom_pipe, int x_pipe, float x, float* y, int no_instances){
    for(int i=0; i<no_instances; i++){ // TODO: zrównoleglenie na i
        float x1 = x_pipe - x;
        float x2 = y_bottom_pipe - y[i];
        float x3 = y[i] - y_upper_pipe;
        vect_in[blocks_nodes[i] + 1] = x1;
        vect_in[blocks_nodes[i] + 2] = x2;
        vect_in[blocks_nodes[i] + 3] = x3;
    }
}

void symulate(int *blocks_nodes, int* column_idx, int* row_pointers, float* weights, int* rewards, int vector_size, int no_instances, int max_iterations){
    // inicjalizacja
    Pipe now = Pipe(300,HEIGHT-(300 + 2*GAP_R),440,40,1);
    Pipe next = Pipe(120,HEIGHT-(120 + 2*GAP_R),880,40,1);
    Pipe prev = Pipe(300,HEIGHT-(300 + 2*GAP_R),440,40,1);
    // bool* isFlopping;
    // isFlopping = (bool*) malloc((no_instances) * sizeof(bool));
    float x = 120; // takie same x dla populaci
    float* y;
    y = (float*) malloc((no_instances) * sizeof(float));
    for(int i=0; i<no_instances; i++) y[i] = 50;

    float* vel_y;
    vel_y = (float*) malloc((no_instances) * sizeof(float));
    for(int i=0; i<no_instances; i++) vel_y[i] = 0;

    bool* colision;// inicialize to false
    colision = (bool*) malloc((no_instances) * sizeof(bool));
    for(int i=0; i<no_instances; i++) colision[i] = false;

    // internal state vectors:
    float* vect_in;
    float* vect_out;
    float* vect_temp;

    vect_in = (float*) malloc((blocks_nodes[no_instances]) * sizeof(float));
    vect_out = (float*) malloc((blocks_nodes[no_instances]) * sizeof(float));
    for(int i=0; i<blocks_nodes[no_instances]; i++) vect_out[i] = 0;
    for(int i=0; i<blocks_nodes[no_instances]; i++) vect_in[i] = 0;
    // koniec inicjalizacji
    
    for(int it = 0; it < max_iterations; it++){ // główna pętla

        update_colisions(colision, rewards, y, x, now.getYUpper(), now.getYBottom(), now.getX(), now.getR(), no_instances);
        bool STOP = false;
        for(int s=0; s<no_instances; s++) STOP |= ~colision[s]; // maybe zrównoleglenie? ale nie trzeba
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
        vect_temp = vect_in;
        vect_in = vect_out;
        vect_out = vect_temp; // swap important as not to allocate memory (everything on device side)
        for(int i=0; i<blocks_nodes[no_instances]; i++) vect_out[i] = 0;

        update_instance_inputs(blocks_nodes, vect_in, now.getYUpper(), now.getYBottom(), now.getX(), x, y, no_instances);
        SparseMUL(column_idx, row_pointers, weights, vect_in, blocks_nodes[no_instances], vect_out, blocks_nodes[no_instances]);

        // #### update step ####
        update_step(vel_y, y, vect_out, blocks_nodes, no_instances);

        prev.update();
        now.update();
        next.update();
    }
    free(y);
    free(vel_y);
    free(colision);
    free(vect_in);
    free(vect_out);
}

void symulation_test(){
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
    // TO CSR MATRIX
    int colsize = count_col_size(enabled, blocks_edges[no_instances]);

    int *col_idx; // size w
    float *weights; // size w
    int *row_pointers; //size translation+1

    col_idx = (int*) malloc((colsize) * sizeof(int));
    weights = (float*) malloc((colsize) * sizeof(float));
    row_pointers = (int*) malloc((blocks_nodes[no_instances] + 1) * sizeof(int));
    for(int i = 0; i<(blocks_nodes[no_instances] + 1); i++){
        row_pointers[i] = 0;
    }
    updateRowPointers(row_pointers, out, blocks_edges, blocks_nodes, no_instances, enabled);
    int *row_pointers_t;
    row_pointers_t = (int*) malloc((blocks_nodes[no_instances] + 1) * sizeof(int));
    for(int i=0; i<(blocks_nodes[no_instances] + 1); i++){
        row_pointers_t[i] = row_pointers[i];
    }
    updateCol_idx_weights(row_pointers_t, weights, col_idx, in, w, out, blocks_edges, blocks_nodes, no_instances,enabled);
    free(row_pointers_t); // important free
    int* rewards;
    rewards = (int*) malloc((no_instances) * sizeof(int));
    //set to -1
    for(int i=0; i<no_instances; i++) rewards[i] = -1;
    // symulate
    symulate(blocks_nodes, col_idx, row_pointers, weights, rewards, blocks_nodes[no_instances], no_instances, 100); // max iteration to 100
    printf("Ended succesfuly!");

}

int main(){
    symulation_test();
}
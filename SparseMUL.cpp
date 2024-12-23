#include <stdio.h>
#include <stdlib.h>

// Mnożenie macierzy CSR przez wektor 
void SparseMUL(int* column_idx, int* row_pointers, float* weights, float* input_vector, int vector_size, float* output_vector, int output_vector_size){
    for(int row=0; row<output_vector_size; row++){ // To do zrównoleglenia
        // co będzie w kernelu
        // kopia __shared__ części wejściowego wektora
        float acc = 0;
        for(int col_idx=row_pointers[row]; col_idx<row_pointers[row+1]; col_idx++){
            acc += input_vector[column_idx[col_idx]] * weights[col_idx];
        }
        output_vector[row] = acc;
    }
}

int testSparseMUL(){
    FILE *plik = fopen("CSR_matrix_mul_test.txt", "r"); // number of elements in vector -> vector -> number of elements weights -> weights ->
                                                        // number of elements in column indices -> column indices -> number of elements Row poiters -> Row pointers
    if (plik == NULL) {
        return 1;
    }
    int* column_idx;
    int* row_pointers;
    float* weights;
    float* input_vector;
    float* output_vect;
    int input_vector_size;
    int weights_size;
    int column_idx_size;
    int row_pointer_size;
    fscanf(plik, "%d", &input_vector_size);
    input_vector = (float*) malloc(input_vector_size * sizeof(float));
    for(int i = 0; i<input_vector_size; i++){
        fscanf(plik, "%f", input_vector+i);
    }
    fscanf(plik, "%d", &weights_size);
    weights = (float*) malloc(weights_size * sizeof(float));
    for(int i = 0; i<weights_size; i++){
        fscanf(plik, "%f", weights+i);
    }
    fscanf(plik, "%d", &column_idx_size);
    column_idx = (int*) malloc(column_idx_size * sizeof(int));
    for(int i = 0; i<column_idx_size; i++){
        fscanf(plik, "%d", column_idx+i);
    }
    fscanf(plik, "%d", &row_pointer_size);
    row_pointers = (int*) malloc(row_pointer_size * sizeof(int));
    for(int i = 0; i<row_pointer_size; i++){
        fscanf(plik, "%d", row_pointers+i);
    }
    output_vect = (float*) malloc(row_pointer_size * sizeof(float));
    SparseMUL(column_idx, row_pointers, weights, input_vector, input_vector_size, output_vect, row_pointer_size-1);
    for(int i=0; i<row_pointer_size-1; i++){ // poprawne 24 3 0.5 0 10 5 0
        printf("%f\t",output_vect[i]);
    }
    fclose(plik);

    free(column_idx);
    free(row_pointers);
    free(weights);
    free(input_vector);
    free(output_vect);

    return 0;
}





int main(){
    testSparseMUL();
    return 0;
}
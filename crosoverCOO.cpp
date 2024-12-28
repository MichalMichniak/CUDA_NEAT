#include <stdio.h>
#include <stdlib.h>


void countOffsprings(int *first_pair, int *second_pair, int no_offsprings, int *innov, int *blocks_edges, int no_instances, int *length_offspring, int offset)
{
    for(int i=0; i<no_offsprings; i++){ // tutaj zrównoleglenie TODO:
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
    }
}

void countOffspringsNodes(int *first_pair, int *second_pair, int no_offsprings, int *translation, int *blocks_nodes, int no_instances, int *length_offspringNodes, int offset){
    for(int i=0; i<no_offsprings; i++){ // tutaj zrównoleglenie TODO:
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
    }
}

void countSurvivors(int *blocks_edges, int *blocks_nodes, int *new_blocks_edges, int *new_blocks_nodes, int no_survivors, int *no_instance_seq){
    for(int i=0; i<no_survivors; i++){// zrównoleglenie TODO:
        int instance_number = no_instance_seq[i];
        new_blocks_nodes[1+i] = blocks_nodes[instance_number+1]-blocks_nodes[instance_number];
        new_blocks_edges[1+i] = blocks_edges[instance_number+1]-blocks_edges[instance_number];
    }
}

void copySurvivors(int *istance_numbers_seq, int *blocks_edges, int *blocks_nodes, int *new_blocks_nodes, int *new_blocks_edges, int no_survivors,
    int *new_in, int *new_out, float *new_w, bool *new_enabled, int *new_innov, int *new_translation, 
    int *in, int *out, float *w, bool *enabled, int *innov, int *translation){

    for(int i=0; i<no_survivors; i++){ // tutaj zrównoleglenie TODO:
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
    }
}

void crossover(
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
    for(int i=0; i<no_offsprings; i++){ // tutaj zrównoleglenie TODO:
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
    printf("\nold block nodes: ");

    for(int i=0; i<6; i++){
        printf("%d\t", blocks_nodes[i]);
    }
    printf("\nold block edges: ");
    for(int i=0; i<6; i++){
        printf("%d\t", blocks_edges[i]);
    }
    printf("\ninnovation numbers: ");
    for(int i=0; i<24; i++){
        printf("%d\t", innov[i]);
    }
    printf("\n");

    int *istance_numbers_seq; // tablica numeru instancji mapowanej ze starej do nowej tablicy szie: no_survivors
    istance_numbers_seq = (int*)malloc(no_survivors*sizeof(int));
    int idx = 0;
    for(int i=0; i<no_instances; i++){ // to będzie sekwencyjne bo tylko przepisanie
        if(mask[i] == 1){
            istance_numbers_seq[idx] = i;
            idx++;
        }
    }
    
    int *new_blocks_nodes;
    int *new_blocks_edges;
    new_blocks_nodes = (int*)malloc((1+no_survivors+no_offsprings+no_mutations)*sizeof(int)); // 0 | 1 do no_survivors | no_survivors + 1 do no_survivors + no_offsprings | no_survivors + no_offsprings + 1 do no_survivors + no_offsprings + no_mutation 
    new_blocks_edges = (int*)malloc((1+no_survivors+no_offsprings+no_mutations)*sizeof(int)); 
    new_blocks_edges[0] = 0;
    new_blocks_nodes[0] = 0;
    
    countSurvivors(blocks_edges, blocks_nodes, new_blocks_edges, new_blocks_nodes, no_survivors, istance_numbers_seq); // uzupełnienie rozmiarów instancji survivors  (0 | 1 do no_survivors)

    countOffsprings(first_pair, second_pair, no_offsprings, innov, blocks_edges, no_instances, new_blocks_edges, no_survivors + 1);
    
    countOffspringsNodes(first_pair, second_pair, no_offsprings, translation, blocks_nodes, no_instances, new_blocks_nodes, no_survivors + 1);
    // TODO: mutations
    
    
    for(int i = 0; i<no_survivors+no_offsprings+no_mutations; i++){// liczenie histogramu skumulowanego (można zrównoleglić)
        new_blocks_nodes[i+1] = new_blocks_nodes[i+1] + new_blocks_nodes[i];
        new_blocks_edges[i+1] = new_blocks_edges[i+1] + new_blocks_edges[i];
    }

    int *new_in;
    int *new_out;
    float *new_w;
    bool *new_enabled;
    int *new_innov;
    int new_no_instances; // survivors + offsprings + mutated
    int *new_translation; // 
    new_no_instances = no_survivors + no_offsprings + no_mutations;
    new_in = (int*)malloc(new_blocks_edges[new_no_instances]*sizeof(int));
    new_out = (int*)malloc(new_blocks_edges[new_no_instances]*sizeof(int));
    new_w = (float*)malloc(new_blocks_edges[new_no_instances]*sizeof(float));
    new_enabled = (bool*)malloc(new_blocks_edges[new_no_instances]*sizeof(bool));
    new_innov = (int*)malloc(new_blocks_edges[new_no_instances]*sizeof(int));
    new_translation = (int*)malloc(new_blocks_nodes[new_no_instances]*sizeof(int));
    // tytaj będzie przepisywanie survivorsów
    copySurvivors(istance_numbers_seq, blocks_edges, blocks_nodes, new_blocks_nodes, new_blocks_edges, no_survivors,
    new_in, new_out, new_w, new_enabled, new_innov, new_translation, in, out, w, enabled, innov, translation);
    
    // tytaj będzie tworzenie potomków
    int *translation_t;
    translation_t = (int*)malloc(blocks_nodes[no_instances]*sizeof(int)); // temporary translation table for crossover
    crossover(
    blocks_nodes,
    blocks_edges,
    new_blocks_nodes,
    new_blocks_edges,
    no_survivors, // offset
    in,
    out,
    w,
    enabled,
    innov,
    first_pair,
    second_pair,
    no_offsprings,
    new_in,
    new_out,
    new_w,
    new_enabled,
    new_innov,
    new_translation,
    translation,
    translation_t 
    );

    printf("\nnew block nodes: ");

    for(int i=0; i<7; i++){
        printf("%d\t", new_blocks_nodes[i]);
    }
    printf("\nnew block edges: ");
    for(int i=0; i<7; i++){
        printf("%d\t", new_blocks_edges[i]);
    }
    printf("\ntemp translation: ");
    for(int i=0; i<blocks_nodes[no_instances]; i++){
        printf("%d\t", translation_t[i]);
    }
    printf("\nnew translation: ");
    for(int i=0; i<new_blocks_nodes[no_survivors+no_offsprings+no_mutations]; i++){
        printf("%d\t", new_translation[i]);
    }
    printf("\nnew in: ");
    for(int i=0; i<new_blocks_edges[no_survivors+no_offsprings+no_mutations]; i++){
        printf("%d\t", new_in[i]);
    }

    printf("\nnew out: ");
    for(int i=0; i<new_blocks_edges[no_survivors+no_offsprings+no_mutations]; i++){
        printf("%d\t", new_out[i]);
    }

    printf("\nnew innov: ");
    for(int i=0; i<new_blocks_edges[no_survivors+no_offsprings+no_mutations]; i++){
        printf("%d\t", new_innov[i]);
    }

    printf("\nnew enabled: ");
    for(int i=0; i<new_blocks_edges[no_survivors+no_offsprings+no_mutations]; i++){
        printf("%d\t", (int)new_enabled[i]);
    }
    printf("\n");
    // tutaj free TODO:

}




int main(){
    test_crosover();
}
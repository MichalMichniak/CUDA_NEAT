#include <stdio.h>
#include <stdlib.h>


int countOffsprings(int *first_pair, int *second_pair, int no_offsprings, int *innov, int *blocks_edges, int no_instances, int *length_offspring, int offset)
{
    for(int i=0; i<no_offsprings; i++){ // tutaj zrównoleglenie
        int first = first_pair[i];
        int second = second_pair[i];
        int idx_first = blocks_edges[first];
        int idx_second = blocks_edges[second];
        int acc = 0;
        while(idx_first != blocks_edges[first+1] && idx_second != blocks_edges[second+1]){
            if(innov[idx_first] == innov[idx_second]){
                idx_first++;
                idx_second++;
            }else if(innov[idx_first] < innov[idx_second]){
                idx_first++;
            }else{
                idx_second++;
            }
            acc++;
        }
        while(idx_first != blocks_edges[first+1]){
            idx_first++;
            acc++;
        }
        while(idx_second != blocks_edges[second+1]){
            idx_second++;
            acc++;
        }
        length_offspring[i+offset]=acc;
    }
}

void countOffspringsNodes(int *first_pair, int *second_pair, int no_offsprings, int *translation, int *blocks_nodes, int no_instances, int *length_offspringNodes, int offset){
    for(int i=0; i<no_offsprings; i++){ // tutaj zrównoleglenie
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
    for(int i=0; i<no_survivors; i++){// zrównoleglenie
        int instance_number = no_instance_seq[i];
        new_blocks_nodes[1+i] = blocks_nodes[instance_number+1]-blocks_nodes[instance_number];
        new_blocks_edges[1+i] = blocks_edges[instance_number+1]-blocks_edges[instance_number];
    }
}

void copySurvivors(int *istance_numbers_seq, int *blocks_edges, int *blocks_nodes, int *new_blocks_nodes, int *new_blocks_edges, int no_survivors,
    int *new_in, int *new_out, float *new_w, bool *new_enabled, int *new_innov, int *new_translation, 
    int *in, int *out, float *w, bool *enabled, int *innov, int *translation){

    for(int i=0; i<no_survivors; i++){ // tutaj zrównoleglenie
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
    int *translation_t
    ){
    

}


void test_crosover(){
    int *in;
    int *out;
    float *w;
    bool *enabled; // not yet implemented
    int *innov; // not yet implemented
    
    int no_instances;
    int *blocks_edges;
    int *translation; // zał że w instancjach posortowane rosnąco
    int *blocks_nodes;
    int *mask; // mask which survives
    int no_survivors; // number of 1 in mask
    int no_mutations = 0;

    int no_offsprings; // ze starej populacji
    int *first_pair; // first parent: size no_offsprings
    int *second_pair; // second parent: size no_offsprings
    int *length_offspring; // size no_offsprings
    int *length_offspring_nodes; // size no_offsprings

    int *length_survivors; // size no_survivors
    int *length_survivors_nodes; // size no_survivors


    
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
    bool *new_enabled; // not yet implemented
    int *new_innov; // not yet implemented
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



}




int main(){
    test_crosover();
}
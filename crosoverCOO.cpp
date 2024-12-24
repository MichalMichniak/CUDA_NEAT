


int countOffsprings(int *first_pair, int *second_pair, int no_offsprings, int *innov, int *blocks_edges, int no_instances, int *length_offspring)
{
    int sum = 0;
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
        length_offspring[i]=acc;
        sum += acc;
    }
    return sum;
}

int countOffspringsNodes(int *first_pair, int *second_pair, int no_offsprings, int *translation, int *blocks_nodes, int no_instances, int *length_offspringNodes){
    int sum = 0;
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
        length_offspringNodes[i]=acc;
        sum += acc;
    }
    return sum;
}

int countSurvivors(){
    
}

int calcNewBlocks(int *mask, int *blocks_edges, int *blocks_nodes,int length){

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
    int no_survivors;
    int no_mutations = 0;

    int no_offsprings; // ze starej populacji
    int *first_pair; // first parent: size no_offsprings
    int *second_pair; // second parent: size no_offsprings
    int *length_offspring; // size no_offsprings
    int *length_offspring_nodes; // size no_offsprings

    int *length_survivors; // size no_survivors
    int *length_survivors_nodes; // size no_survivors


    int count = countOffsprings(first_pair, second_pair, no_offsprings, innov, blocks_edges, no_instances, length_offspring);
    int count_nodes = countOffspringsNodes(first_pair, second_pair, no_offsprings, translation, blocks_nodes, no_instances, length_offspring_nodes);
    
    int *new_in;
    int *new_out;
    float *new_w;
    bool *new_enabled; // not yet implemented
    int *new_innov; // not yet implemented
    int new_no_instances; // survivors + offsprings + mutated
    int *new_blocks_edges;
    int *new_translation; // 
    int *new_blocks_nodes;
    int length; // length_survivors + count (offsprings) + mutation_length
    // tytaj będzie przepisywanie survivorsów


}




int main(){
    test_crosover();
}
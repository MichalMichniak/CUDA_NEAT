#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstdlib>
#include <cmath>

#define WIDTH 800
#define HEIGHT 600
#define GAP_R 80
class FloppyBird{
private:
    float flop_v = 4;
    float g = 0.1;
    float x,y;
    float hit_box_w, hit_box_h;
    float vel_x,vel_y;
    int max_h,min_h;
    sf::CircleShape m_shape;


public:
    FloppyBird(int x_, int y_, int vel_x_, int vel_y_){ 
        x = x_; y=y_; vel_x=vel_x_; vel_y = vel_y_;
        m_shape = sf::CircleShape(10);
        m_shape.setFillColor(sf::Color(255, 255, 0));
        m_shape.setOrigin(m_shape.getRadius(), m_shape.getRadius());
        hit_box_h = m_shape.getRadius();
        hit_box_w = m_shape.getRadius();
        }
    bool get_colision(int y_upper_pipe, int y_bottom_pipe, int x_pipe, int r_pipe){
        if(((x + hit_box_w >=x_pipe-r_pipe && x-hit_box_w<=x_pipe+r_pipe) && (y-hit_box_h<=y_upper_pipe || y+hit_box_h>=HEIGHT-y_bottom_pipe)) || (y>=HEIGHT-hit_box_h || y<=0+hit_box_h)){
            return true;
        }
        return false;
    }
    void step(bool flop){
        vel_y += g;
        if(flop){
            vel_y = -flop_v;
        }
        // symulation
        x+=vel_x;
        y+=vel_y;
        m_shape.setPosition(x,y);
    }
    virtual void draw(sf::RenderWindow& win){
        win.draw(m_shape);
    }
    float get_y(){ return y;}
    float get_x(){ return x;}
    float get_hit_box_w(){ return hit_box_w;}

};

class Pipe{
private:
    int y_upper_pipe;
    int y_bottom_pipe;
    int x_pipe;
    int r_pipe;
    int vel;
    sf::RectangleShape upperPipeShape;
    sf::RectangleShape lowerPipeShape;

public:
    Pipe(int y_upper_pipe_, int y_bottom_pipe_, int x_pipe_, int r_pipe_, int vel_) : y_bottom_pipe(y_bottom_pipe_), y_upper_pipe(y_upper_pipe_), x_pipe(x_pipe_), r_pipe(r_pipe_), vel(vel_) {
        upperPipeShape = sf::RectangleShape();
        lowerPipeShape = sf::RectangleShape();
        upperPipeShape.setSize(sf::Vector2f(r_pipe*2, y_upper_pipe));
        upperPipeShape.setFillColor(sf::Color(0, 255, 0));
        upperPipeShape.setPosition(x_pipe - r_pipe, 0);

        lowerPipeShape.setSize(sf::Vector2f(r_pipe*2, y_bottom_pipe));
        lowerPipeShape.setFillColor(sf::Color(0, 255, 0));
        lowerPipeShape.setPosition(x_pipe - r_pipe, HEIGHT - y_bottom_pipe);
    };

    Pipe(const Pipe& other) 
        : y_upper_pipe(other.y_upper_pipe), y_bottom_pipe(other.y_bottom_pipe),
          x_pipe(other.x_pipe), r_pipe(other.r_pipe), vel(other.vel),
          upperPipeShape(other.upperPipeShape), lowerPipeShape(other.lowerPipeShape) {
    }
    bool update(){
        x_pipe -= vel;

        upperPipeShape.setPosition(x_pipe - r_pipe, 0);
        lowerPipeShape.setPosition(x_pipe - r_pipe, HEIGHT - y_bottom_pipe);
        
        if(x_pipe > -10) return 0;
        return 1;
    };

    virtual void draw(sf::RenderWindow& win){
        
        win.draw(upperPipeShape);
        win.draw(lowerPipeShape);
    }
    int getX() const { return x_pipe; }
    int getYUpper() const { return y_upper_pipe; }
    int getYBottom() const { return y_bottom_pipe; }
    int getR() const { return r_pipe; }
    int getVelocity() const { return vel; }
    ~Pipe() = default;
};


int count_col_size(bool* enable, int length){
    int acc = 0;
    for(int i=0; i<length; i++) if (enable[i]) acc++;
    return acc;
}

void updateRowPointers(int *row_pointers, int *out, int no_edges, int no_nodes, bool* enable){
    for(int i=0; i<no_edges; i++){
        if (enable[i])
            row_pointers[(out[i])+1] +=1;
    }
    for(int i=0; i<no_nodes; i++){
        row_pointers[i+1] += row_pointers[i];
    }
}

void updateCol_idx_weights(int *row_pointers_t, float *weights, int *col_idx, int *in, float *w, int *out, int no_edges, int no_nodes, bool* enable){

    for(int i=0; i<no_edges; i++){
        if (enable[i]){
            int temp = row_pointers_t[out[i]];
            row_pointers_t[out[i]] +=1; // atomicAdd
            col_idx[temp] = in[i]; // jeżeli będą podwójne to nie ma znaczenia przy mnożeniu macierzowym
            weights[temp] = w[i];
        }
    }
}

void update_instance_inputs(float* vect_in, int y_upper_pipe, int y_bottom_pipe, int x_pipe, float x, float y){
    float x1 = (x_pipe - x)/400.0;
    float x2 = (y_bottom_pipe - y)/600.0;
    float x3 = (y - y_upper_pipe)/600.0;
    vect_in[1] = x1;
    vect_in[2] = x2;
    vect_in[3] = x3;

}

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

int main(){
    FILE *plik = fopen("best_net.txt", "r");
    if (plik == NULL) {
        return  1;
    }
    int no_edges;
    int no_nodes;
    int *in;
    int *out;
    float *w;
    bool *enabled;

    // init CSR
    fscanf(plik, "%d", &no_nodes);
    fscanf(plik, "%d", &no_edges);
    in = (int*) malloc(no_edges * sizeof(int));
    for(int i = 0; i<no_edges; i++){
        fscanf(plik, "%d", in+i);
    }
    out = (int*) malloc(no_edges * sizeof(int));
    for(int i = 0; i<no_edges; i++){
        fscanf(plik, "%d", out+i);
    }
    w = (float*) malloc(no_edges * sizeof(float));
    for(int i = 0; i<no_edges; i++){
        fscanf(plik, "%f", w+i);
    }
    enabled = (bool*) malloc(no_edges * sizeof(int));
    for(int i = 0; i<no_edges; i++){
        int temp;
        fscanf(plik, "%d", &temp);
        *(enabled+i) = (bool)temp;
    }

    fclose(plik);
    int colsize = count_col_size(enabled, no_edges);

    int *col_idx; // size w
    float *weights; // size w
    int *row_pointers; //size translation+1

    col_idx = (int*) malloc((colsize) * sizeof(int));
    weights = (float*) malloc((colsize) * sizeof(float));
    row_pointers = (int*) malloc((no_nodes + 1) * sizeof(int));

    for(int i = 0; i<(no_nodes + 1); i++){
        row_pointers[i] = 0;
    }

    updateRowPointers(row_pointers, out, no_edges, no_nodes, enabled);

    int *row_pointers_t;
    row_pointers_t = (int*) malloc((no_nodes + 1) * sizeof(int));
    for(int i=0; i<(no_nodes + 1); i++){
        row_pointers_t[i] = row_pointers[i];
    }

    updateCol_idx_weights(row_pointers_t, weights, col_idx, in, w, out, no_edges, no_nodes, enabled);
    // end init CSR

    // internal state vectors:
    float* vect_in;
    float* vect_out;
    float* vect_temp;

    vect_in = (float*) malloc((no_nodes) * sizeof(float));
    vect_out = (float*) malloc((no_nodes) * sizeof(float));
    for(int i=0; i<no_nodes; i++) vect_out[i] = 0;
    for(int i=0; i<no_nodes; i++) vect_in[i] = 0;
    // koniec inicjalizacji sieci

    srand(time(0));
    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "My window");
    FloppyBird fb = FloppyBird(120,50,0,0);
    Pipe now = Pipe(300,HEIGHT-(300 + 2*GAP_R),440,40,1);
    Pipe next = Pipe(120,HEIGHT-(120 + 2*GAP_R),880,40,1);
    Pipe prev = Pipe(300,HEIGHT-(300 + 2*GAP_R),440,40,1);
    bool isFlopping = false;
    while (window.isOpen())
    {
        // check all the window's events that were triggered since the last iteration of the loop
        sf::Event event;
        while (window.pollEvent(event))
        {
            // "close requested" event: we close the window
            if (event.type == sf::Event::Closed)
                window.close();
            
        }
        window.clear(sf::Color::Black);
        // if(fb.get_y() > 100){
        //     fb.step(true);
        // }else{
        //     fb.step(false);
        // }
        // printf("%d\n", fb.get_colision())
        bool colision = fb.get_colision(now.getYUpper(), now.getYBottom(), now.getX(), now.getR());
        // printf("%d\n", colision);
        if(colision) return 0;
        if(fb.get_x() - fb.get_hit_box_w() > now.getX() + now.getR()){
            prev = now;
            now = next;
            int randomNum = rand() %(HEIGHT - 2*GAP_R) + GAP_R;
            next = Pipe(randomNum - GAP_R,HEIGHT - (randomNum + GAP_R),880,40,1);
        }
        vect_temp = vect_in;
        vect_in = vect_out;
        vect_out = vect_temp;
        for(int i=0; i<no_nodes; i++) vect_out[i] = 0;
        update_instance_inputs(vect_in, now.getYUpper(), now.getYBottom(), now.getX(), fb.get_x(), fb.get_y());

        SparseMUL(col_idx, row_pointers, weights, vect_in, no_nodes, vect_out, no_nodes);

        isFlopping = (vect_out[0]>0.5);

        fb.step(isFlopping);
        prev.update();
        now.update();
        next.update();
        fb.draw(window);
        prev.draw(window);
        now.draw(window);
        next.draw(window);
        window.display();
        sf::sleep(sf::seconds(0.01));
    }
    return 0;
}
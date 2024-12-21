#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>

#define WIDTH 800
#define HEIGHT 600
#define GAP_R 140
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




int main(){
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
            if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Space) { // Spacebar to make the bird flap
                    isFlopping = true;
                }
            }
            if (event.type == sf::Event::KeyReleased) {
                if (event.key.code == sf::Keyboard::Space) { // Stop flapping when key is released
                    isFlopping = false;
                }
            }
            
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
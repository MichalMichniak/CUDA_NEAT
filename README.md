# CUDA_NEAT


## SFML Enviroment:
Pod Windows:

1. Pobierz pliki SFML ze [strony](https://www.sfml-dev.org/download/)

2.  ``` shell
    g++ -c manual_control_env.cpp -I <ścierzka do SFML/include>
    ```

3.  ``` shell
    g++ manual_control_env.o -o manual_control_env -L <ścierzka do SFML\lib> -lsfml-graphics -lsfml-window -lsfml-system -lopengl32 -lwinmm -lgdi32
    ```
4. ``` shell
   ./manual_control_env.exe
   ```

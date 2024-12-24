# CUDA_NEAT


## SFML Enviroment:
Pod Windows:

1. Pobierz pliki SFML ze [strony](https://www.sfml-dev.org/download/)

2.  ``` shell
    g++ -c manual_control_env.cpp -I <ścieżka do SFML/include>
    ```

3.  ``` shell
    g++ manual_control_env.o -o manual_control_env -L <ścieżka do SFML\lib> -lsfml-graphics -lsfml-window -lsfml-system -lopengl32 -lwinmm -lgdi32
    ```
4. ``` shell
   ./manual_control_env.exe
   ```

## Algoritm:
1. inicjalizacja populacji startowej (trywialna inicjalizacja tylko wejścia i wyjscia jako wierchołki grafu a połączenia wszystkie wejścia z wyjściami)
    * inicjalizacja innovation number [Host]
    * inicjalizacja połączeń jako kodowanie COO dla tablic weight, in, out [Host] lub [Device]
    * inicjalizacja tablicy translations [Host] lub [Device]
2. Przygotowania do symulacji:
    * zamienienie macierzy COO na CSR (nie jest wymagana zachowana kolejność kolumn w wierszach) [Device]
    * alokacja pamięci dla wejścia i wyjścia (pamięć długości tablicy translations) (wyzerowanie) [Host]
    * inicjalizacja położeń instancji [Host]
    * inicjalizacja wejść [Device]
    * inicjalizacja przeszkód (jak w manual_control_env) [Host]
    * inicjalizacja maksymalnej przebytej drogi (inicjalizacja max int) [Host]
3. Symulacja:
    * przemnożenie macierzy CSR przez wektor wejściowy [Device]
    * wektor wyjściowy przekształcić funkcją aktywacji [Device]
    * wykonać step symulacji na podstawie wyników (wyjścia dla każdej instancji wziąć z wektora blocks_nodes) [Device]
    * dla każdej instancji wywołać get_colision() i jeżeli true to wpisać minimum z tego co w tablicy rewards a dotychczasowo przebytą drogą (ze sprawdzeniem czy wszystkie < max int [kryterium STOP-u]) [Device]
    * wpisać obserwacje do wektora output [Device]
    * przepiąć wskaźniki input na output [Host]
    * sprawdzić kryterium stopu [Host]
4. Alokacja Pamięci:
    * określenie liczby survivors [Host]
    * określenie liczby offsprings [Host]
    * określenie liczby mutations [Host]
    * na podstawie tablicy rewards określenie maski (jaka instancja przetrwa jaka nie) (trzeba rozwinąć) [Host]/[Device]
    * inicjalizacja nowej tablicy new_blocks_nodes, new_blocks_edges [Host]
    * określenie długości dla każdej instancji z survivors i wpisanie do blocks_nodes i blocks_edges. [Device]
    * określenie długości dla każdej instancji z offspring i wpisanie do blocks_nodes i blocks_edges (idea w crosoverCOO.cpp) [Device]
    * podzielenie mutacji na  dodania wierzchołka, mutacje krawędzi i mutacji wag. (odpowiednio (+1 node, +2 edges), (+0 node, +1 edges), (+0 node, +0 edges)) [Host]
    * policzenie i wpisanie długości każdej z mutacji do new_blocks_nodes, new_blocks_edges [Device]
    * policzenie histogramu skumulowanego dla new_blocks_nodes, new_blocks_edges [Device]
    * alokacja [Host]
    ```
    int *new_in;
    int *new_out;
    float *new_w;
    bool *new_enabled;
    int *new_innov;
    int new_no_instances;
    int *new_translation;
    int length;
    ```
> [!IMPORTANT]
> Tablice innov, translation są i mają być posortowane wewnątrz instancji. (duże zmniejszenie złożoności)
5. Nowa Populacja
    * nadanie numeru dla każdego survivor-a [Host]
    * przepisanie względem maski i new_blocks_nodes i new_blocks_edges survivorów [Device]
    * stworzenie i wpisanie do nowej populacji offspring (work in progress) [Device]
    * podobnie mutacje [Device]
    * powrót do punktu 2. jeżeli nie osiągnięto max liczby iteracji.


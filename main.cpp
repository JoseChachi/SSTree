#include <iostream>
#include <vector>
#include <random>
#include "SStree.h"
#include <fstream>

int main() {
    // Create a random number generator
    SsTree tree;
    std::vector<Point> points;

    std::ifstream file("../dataset.txt", std::ios::in); // x y
    while (!file.eof()) {
        Point point(2);
        double x, y;
        file >> x >> y;
        point[0] = x;
        point[1] = y;
        points.push_back(point);
    }

    // Insert the points into the SStree
    for (auto& point : points) {
        //std::cout << point[0] << std::endl;
        tree.insert(point);
    }
    tree.test();
    tree.print();
    /*
    int k = 10;
    std::cout << "Linear searching\n";
    std::vector<NType> linearSearching;
    for(const auto& point: points) {
        linearSearching.push_back(distance(point,points[0]));
    }
    sort(linearSearching.begin(),linearSearching.end());

    for(int i = 0; (i < k) && (i < points.size()); i++){
        std::cout << linearSearching[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "KNN Searching" << std::endl;
    std::vector<Point> knnpoints = tree.busquedaKNN(points[0], k);

    for(int i = 0; (i < k) && (i < points.size()); i++){
        std::cout << distance(knnpoints[i], points[0]) << "\t";
        std::cout << knnpoints[i] << std::endl;
    }
    */
    sf::RenderWindow window(sf::VideoMode(800, 800), "SSTree");
    while (window.isOpen()) {
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }
        tree.draw(window);
        window.display();
        window.clear();
    }
    return 0;
}

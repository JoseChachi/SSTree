#include <iostream>
#include <vector>
#include <random>
#include "Point.h"
#include "SStree.h"
#include <json.hpp>
#include <fstream>


struct ImageData {
    std::vector<Point> embeddings;
    std::vector<std::string> paths;
};

ImageData readEmbeddingsFromJson(const std::string& FILE_NAME) {
    ImageData data;

    try {
        std::ifstream file(FILE_NAME);
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open JSON file.");
        }

        nlohmann::json jsonData;
        file >> jsonData;

        std::vector<std::vector<float>> features;
        for (const auto& featureList : jsonData["features"]) {
            std::vector<float> tempFeature = featureList;
            features.push_back(tempFeature);
        }

        data.paths = jsonData.at("paths").get<std::vector<std::string>>();

        if (features.size() != data.paths.size()) {
            throw std::runtime_error("The number of features does not match the number of paths.");
        }

        for (const auto& feature : features) {
            Point embedding(feature.size());
            for (size_t j = 0; j < feature.size(); ++j) {
                embedding[j] = feature[j];
            }
            data.embeddings.push_back(embedding);
        }

        file.close();

        for (const auto& path : data.paths) {
            if (path.empty()) {
                throw std::runtime_error("Uno de los paths está vacío.");
            }
            std::cout << "Guardado path: " << path << std::endl;
        }

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    return data;
}


int main() {
    const std::string FILE_NAME("../embedding.json");
    ImageData data = readEmbeddingsFromJson(FILE_NAME);

    SsTree tree;
    tree.print();
    for (size_t i = 0; i < data.embeddings.size(); ++i) {
        tree.insert(data.embeddings[i], data.paths[i]);
        if((!(i%250)) || (i == (data.embeddings.size()-1))) {
            std::cout << i << "\tNivel0: ";
            std::cout << tree.root->radius;

            if(!(tree.root->isLeaf())){
                auto* innerNode = dynamic_cast<SsInnerNode*>(tree.root);
                std::vector<SsNode*> hejos = innerNode->children;
                NType mx = 0;
                for(const auto& child: hejos) {
                    mx = max(mx, child->radius);
                }
                std::cout << "\tNivel1: " << mx << std::endl;
            }
            else{
                std::cout << std::endl;
            }


        }
    }
    std::string filename = "embbeding.dat";
    tree.saveToFile(filename);
    int k = 10;
    std::cout << "Linear searching\n";
    std::vector<NType> linearSearching;
    for(const auto& point: data.embeddings) {
        linearSearching.push_back(distance(point,data.embeddings[0]));
    }
    sort(linearSearching.begin(),linearSearching.end());

    for(int i = 0; (i < k+1) && (i < data.embeddings.size()); i++){
        std::cout << linearSearching[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "KNN searching k =" << k << std::endl;
    for(const auto& puntito: tree.busquedaKNN(data.embeddings[0],11)){
        std::cout << distance(data.embeddings[0],puntito) << "\t";
    }
    std::cout << std::endl;
    return 0;
}
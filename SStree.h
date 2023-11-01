#ifndef SSTREE_H
#define SSTREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include "params.h"
#include "Point.h"
#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <queue>

class SsNode {
private:
    NType varianceAlongDirection(const std::vector<Point>& centroids, size_t direction) const;
    size_t minVarianceSplit(std::vector<NType> points, size_t coordinateIndex);
    
public:
    virtual ~SsNode() = default;

    Point centroid; 
    NType radius;
    SsNode* parent = nullptr;

    virtual bool isLeaf() const = 0;
    virtual std::vector<Point> getEntriesCentroids() const = 0;
    virtual void sortEntriesByCoordinate(size_t coordinateIndex) = 0;
    virtual std::pair<SsNode*, SsNode*> split() = 0;
    virtual bool intersectsPoint(const Point& point) const {
        return distance(this->centroid, point) <= this->radius;
    }
    void draw(sf::RenderWindow &window) const;
    virtual void updateBoundingEnvelope() = 0;
    size_t directionOfMaxVariance() const;
    size_t findSplitIndex();

    virtual std::pair<SsNode*, SsNode*> insert(const Point& point) = 0;

    bool test(bool status) const;
    void print(size_t indent = 0) const;

    void draw(sf::RenderWindow window);

    virtual void saveToStream(std::ostream &out) const = 0;
    virtual void loadFromStream(std::istream &in) = 0;
};

class SsInnerNode : public SsNode {
private:
    std::vector<Point> getEntriesCentroids() const override;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    std::pair<SsNode*, SsNode*> split() override;
    std::vector<SsNode*> children;

    SsNode* findClosestChild(const Point& target) const;
    bool isLeaf() const override { return false; }
    void updateBoundingEnvelope() override;

    const float norm = 1.3;
    void draw(sf::RenderWindow &window) const;

    std::pair<SsNode*, SsNode*> insert(const Point& point) override;

    virtual void saveToStream(std::ostream &out) const override;
    virtual void loadFromStream(std::istream &in) override;
};

class SsLeaf : public SsNode {
private:
    std::vector<std::string> paths;

    std::vector<Point> getEntriesCentroids() const override;

    const float norm = 1.3;
    void sortEntriesByCoordinate(size_t coordinateIndex) override;

public:
    std::pair<SsNode*, SsNode*> split() override;
    std::vector<Point> points;

    bool isLeaf() const override { return true; }
    void updateBoundingEnvelope() override;
    void draw(sf::RenderWindow &window) const;
    std::pair<SsNode*, SsNode*> insert(const Point& point) override;

    virtual void saveToStream(std::ostream &out) const override;
    virtual void loadFromStream(std::istream &in) override;
};

//bool compare(const std::pair<NType, Point> &a,
//             const std::pair<NType, Point> &b){
//    return a.first < b.first;
//}

inline auto compare = [](std::pair<NType, Point>& a, std::pair<NType, Point>& b){
    return a.first < b.first;
};

class SsTree {
private:
    SsNode* search(SsNode* node, const Point& target);
    SsNode* searchParentLeaf(SsNode* node, const Point& target);

public:

    SsNode* root;
    SsTree() : root(nullptr) {}
    ~SsTree() {
        delete root;
    }
    
    void insert(const Point& point);
    //void insert(SsNode* node, const Point& point);
    void insert(const Point& point, const std::string& path);
    void build (const std::vector<Point>& points);
    std::vector<Point> busquedaKNN(const Point& target, size_t k);
    void kNNQuery(const Point& target, size_t k, SsNode* node, NType radius, std::priority_queue<std::pair<NType, Point>, std::vector<std::pair<NType, Point>>, decltype(compare)>& result) const;
    void draw(sf::RenderWindow &window) const;

    void print() const;
    void test() const;

    void saveToFile(const std::string &filename) const;
    void loadFromFile(const std::string &filename, size_t D=448);
};

#endif // !SSTREE_H
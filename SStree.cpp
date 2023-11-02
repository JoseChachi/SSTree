#include <numeric>
#include <fstream>
#include "SStree.h"


bool SsNode::test(bool isRoot) const {
    size_t count = 0;
    if (this->isLeaf()) {
        const SsLeaf* leaf = dynamic_cast<const SsLeaf*>(this);
        count = leaf->points.size();

        for (const Point& point : leaf->points) {
            if (distance(this->centroid, point) > this->radius) {
                std::cout << "Point outside node radius detected." << std::endl;
                return false;
            }
        }
    } else {
        const SsInnerNode* inner = dynamic_cast<const SsInnerNode*>(this);
        count = inner->children.size();

        for (const SsNode* child : inner->children) {
            if (distance(this->centroid, child->centroid) > this->radius) {
                std::cout << "Child centroid outside parent radius detected." << std::endl;
                return false;
            }
            if (!child->test(false)) {
                return false;
            }
        }
    }

    if (!isRoot && (count < Settings::m || count > Settings::M)) {
        std::cout << "Invalid number of children/points detected." << std::endl;
        return false;
    }

    if (!isRoot && !parent) {
        std::cout << "Node without parent detected." << std::endl;
        return false;
    }

    return true;
}

void SsTree::test() const {
    bool result = root->test(true);

    if (root->parent) {
        std::cout << "Root node parent pointer is not null!" << std::endl;
        result = false;
    }

    if (result) {
        std::cout << "SS-Tree is valid!" << std::endl;
    } else {
        std::cout << "SS-Tree has issues!" << std::endl;
    }
}


void SsNode::print(size_t indent) const {
    for (size_t i = 0; i < indent; ++i) {
        std::cout << "  ";
    }

    // Imprime información del nodo.
    std::cout << "Centroid: " << centroid << ", Radius: " << radius;
    if (isLeaf()) {
        const SsLeaf* leaf = dynamic_cast<const SsLeaf*>(this);
        std::cout << ", Points: [ ";
        for (const Point& p : leaf->points) {
            std::cout << p << " ";
        }
        std::cout << "]";
    } else {
        std::cout << std::endl;
        const SsInnerNode* inner = dynamic_cast<const SsInnerNode*>(this);
        for (const SsNode* child : inner->children) {
            child->print(indent + 1); 
        }
    }
    std::cout << std::endl;
}

size_t SsNode::directionOfMaxVariance() const {
    NType maxVariance = 0;
    size_t directionIndex = 0;
    std::vector<Point> centroids = this->getEntriesCentroids();
    for (int i = 0; i < this->centroid.dim(); i++){
        if(varianceAlongDirection(centroids,i) > maxVariance){
            maxVariance = varianceAlongDirection(centroids, i);
            directionIndex = i;
        }
    }
    return directionIndex;
}


size_t SsNode::findSplitIndex() {
    size_t coordinateIndex = this->directionOfMaxVariance();
    this->sortEntriesByCoordinate(coordinateIndex);
    std::vector<NType> points = {};
    for(const auto& point: this->getEntriesCentroids()){
        points.push_back(point[coordinateIndex]);
    }
    return this->minVarianceSplit(points, coordinateIndex);
}

NType variance(std::vector<NType> points) {
    NType mean = std::accumulate( points.begin(), points.end(), static_cast<NType>(0.) ) / static_cast<NType>(points.size());
    NType variance = 0;
    for (const auto& point:points) {
        variance += (point - mean) * (point - mean);
    }
    return variance;
};

size_t SsNode::minVarianceSplit(std::vector<NType> points, size_t coordinateIndex) {
    NType minVariance = NType::max_value();
    size_t splitIndex = Settings::m;

    for(int i = Settings::m; i < abs(points.size() - Settings::m); i++){
        NType variance1 = variance(std::vector<NType>(points.begin(), points.begin()+i));
        NType variance2 = variance(std::vector<NType>(points.begin()+1, points.end()));

        if((variance1 + variance2) < minVariance){
            minVariance = variance1+variance2;
            splitIndex = i;
        }
    }

    return splitIndex;
}

NType SsNode::varianceAlongDirection(const std::vector<Point> &centroids, size_t direction) const {
    NType mean = std::accumulate( centroids.begin(), centroids.end(), static_cast<NType>(0.0), [direction](NType sum, const Point& point) { return sum + point[direction]; } ) / static_cast<NType>(centroids.size());
    NType variance = 0;
    for(const auto& point: centroids){
        variance += (point[direction] - mean) * (point[direction] - mean);
    }
    return variance;
}

void SsNode::draw(sf::RenderWindow window) {
    if(!this->isLeaf()){
        auto* innerNode = dynamic_cast<SsInnerNode*>(this);
        innerNode->draw(window);
    }
    else{
        auto* leafNode = dynamic_cast<SsLeaf*>(this);
        leafNode->draw(window);
    }
}



void SsTree::print() const {
    if (root) {
        root->print();
    } else {
        std::cout << "Empty tree." << std::endl;
    }
}

SsNode *SsTree::search(SsNode *node, const Point &target) {
    if (node->isLeaf()){
        const auto* leaf = dynamic_cast<const SsLeaf*>(node);
        for(const auto& point : leaf->points){
            if(point == target) return node;
        }
    }
    else{
        const auto* innerNode = dynamic_cast<const SsInnerNode*>(node);
        for(const auto& childNode: innerNode->children){
            if(childNode->intersectsPoint(target)){
                auto result = this->search(childNode, target);
                if (result != nullptr) return result;
            }
        }
    }
    return nullptr;
}

SsNode *SsTree::searchParentLeaf(SsNode *node, const Point &target) {
    if (node->isLeaf()) return node;
    else {
        const auto* innerNode = dynamic_cast<const SsInnerNode*>(node);
        auto childNode = innerNode->findClosestChild(target);
        return searchParentLeaf(childNode, target);
    }
}

void SsTree::insert(const Point &point) {
    if(this->root == nullptr){
        this->root = new SsLeaf();
        this->root->parent = nullptr;
        this->root->centroid = point;
        this->root->radius = 0;
    }
    std::pair<SsNode*, SsNode*> result = this->root->insert(point);
    if (result.first != nullptr){
        this->root = new SsInnerNode();
        this->root->parent = nullptr;
        result.first->parent = this->root;
        result.second->parent = this->root;
        dynamic_cast<SsInnerNode*>(this->root)->children.emplace_back(result.first);
        dynamic_cast<SsInnerNode*>(this->root)->children.emplace_back(result.second);
    }
}

std::pair<SsNode *, SsNode *> SsInnerNode::insert(const Point &point) {
    auto closestChild = this->findClosestChild(point);
    auto result = closestChild->insert(point);
    if (result.first == nullptr){
        this->updateBoundingEnvelope();
        return std::make_pair(nullptr, nullptr);
    }
    else{
        this->children.erase(std::remove(this->children.begin(), this->children.end(), closestChild), this->children.end());
        this->children.push_back(result.first);
        this->children.push_back(result.second);
        this->updateBoundingEnvelope();
        if (this->children.size() <= Settings::M){
            return std::make_pair(nullptr, nullptr);
        }
    }
    return this->split();
}

std::pair<SsNode *, SsNode *> SsLeaf::insert(const Point &point) {
    for(const auto pointNode: this->points){
        if(pointNode == point) return std::make_pair(nullptr, nullptr);
    }
    this->points.push_back(point);
    this->updateBoundingEnvelope();
//    std::cout << "updating\n";
    if (this->points.size() <= Settings::M){
        return std::make_pair(nullptr, nullptr);
    }
    return this->split();
}

SsNode *SsInnerNode::findClosestChild(const Point &target) const {
    if(this->isLeaf()) throw std::runtime_error("Nani");
    NType minDistance = NType::max_value();
    SsNode* result = nullptr;

    for(auto& child: this->children){
        if (distance(child->centroid, target) < minDistance){
            minDistance = distance(child->centroid, target);
            result = child;
        }
    }
    return result;
}

void SsInnerNode::updateBoundingEnvelope() {
    std::vector<SsNode*> childNodes = this->children;
    this->centroid = childNodes[0]->centroid;
    int dims = this->centroid.dim();
    NType mx = 0;
    NType d = 0;
    auto candA = childNodes[0];
    auto candB = childNodes[0];
    for(const auto& childA: childNodes){
        for(const auto& childB: childNodes) {
            if(childA == childB) continue;
            d = distance(childA->centroid, childB->centroid) + childA->radius + childB->radius;
            if(d > mx){
                mx = d;
                candA = childA;
                candB = childB;
            }
        }
    }

    this->radius = d/2;

    std::vector<NType> unitV;
    d = distance(candA->centroid, candB->centroid);
    NType res;
    for(int i = 0; i < dims; i++){
        res = candA->centroid[i] - candB->centroid[i];
        unitV.push_back(res/d);
    }

    Point lim1(dims), lim2(dims);

    for(int i = 0; i < dims; i++){
        lim1[i] = candA->centroid[i] + (unitV[i] * candA->radius);
        lim2[i] = candB->centroid[i] - (unitV[i] * candB->radius);
    }

    for(int i = 0; i < dims; i++){
        this->centroid[i] = (lim1[i] + lim2[i])/2;
    }

    for(const auto& child: childNodes){
        d = distance(this->centroid, child->centroid) + child->radius;
        if(d <= this->radius) continue;
        this->radius = d;
    }
}

void SsLeaf::updateBoundingEnvelope() {
    std::vector<Point> points = this->getEntriesCentroids();
    this->centroid = points[0];
    int dims = this->centroid.dim();
    NType mx = 0;
    NType d = 0;
    Point candA = Point(dims);
    Point candB = Point(dims);
    for(const auto& pointA: points){
      for(const auto& pointB: points){
          if(pointA == pointB) continue;
          d = distance(pointA, pointB);
          if(d > mx){
              mx = d;
              candA = pointA;
              candB = pointB;
          }
      }
    }

    Point maxDc = Point(dims);
    for(int i = 0; i < dims; i++){
      maxDc[i] = (candA[i] + candB[i])/2;
    }
    NType maxDr = d/2;
    for(const auto& point: points){
      d = distance(maxDc, point);
      if(d <= maxDr) continue;
      maxDr = d;
    }

    Point minMBB = Point(dims);
    Point maxMBB = Point(dims);
    Point midMBB = Point(dims);
    for(int i = 0; i < dims; i++) {
        minMBB[i] = NType::max_value();
        maxMBB[i] = 0;
        for(const auto& dot: points){
            minMBB[i] = min(minMBB[i], dot[i]);
            maxMBB[i] = max(maxMBB[i], dot[i]);
        }

        midMBB[i] = (minMBB[i] + maxMBB[i])/2;
    }

    NType r = 0;
    for(const auto& point: points){
        d = distance(midMBB, point);
        if(d <= r) continue;
        r = d;
    }
    if(r <= maxDr){
        this->centroid = midMBB;
        this->radius = r;
    }
    else{
        this->centroid = maxDc;
        this->radius = maxDr;
    }


    int k = 100;
    std::vector<NType> unitV;
    NType distV;
    distV = distance(midMBB, maxDc);
    if(distV == 0) return;
    NType res;
    for(int i = 0; i < dims; i++){
        res = midMBB[i] - maxDc[i];
        unitV.push_back(res/distV);
    }


    for(int i = 1; i <= k; i++){
        Point newCenter = maxDc;
        NType newRadius = 0;
        for(int j = 0; j < dims; j++){
            newCenter[j] += unitV[j]*(distV/i);
        }
        for(const auto& dot: points){
            d = distance(newCenter, dot);
            if(d > newRadius) newRadius = d;
        }

        if(newRadius <= this->radius){
            this->radius = newRadius;
            this->centroid = newCenter;
        }

    }

}

std::pair<SsNode *, SsNode *> SsInnerNode::split() {
    auto splitIndex = this->findSplitIndex();

    SsNode* newNode1 = new SsInnerNode();
    dynamic_cast<SsInnerNode*>(newNode1)->children = std::vector<SsNode*>(this->children.begin(), this->children.begin() + splitIndex);

    SsNode* newNode2 = new SsInnerNode();
    dynamic_cast<SsInnerNode*>(newNode2)->children = std::vector<SsNode*>(this->children.begin() + splitIndex, this->children.end());

    newNode1->updateBoundingEnvelope();
    newNode2->updateBoundingEnvelope();

    newNode1->parent = this->parent;
    newNode2->parent = this->parent;

    return std::make_pair(newNode1, newNode2);
}

std::vector<Point> SsInnerNode::getEntriesCentroids() const {
    std::vector<Point> childCentroids = {};
    for(const auto& child: this->children){
        childCentroids.push_back(child->centroid);
    }
    return childCentroids;
}


std::vector<Point> SsLeaf::getEntriesCentroids() const {
    return this->points;
}

#include <ctime>
#include <cstdlib>
#include <map>
std::pair<SsNode *, SsNode *> SsLeaf::split() {
    SsNode* newNode1 = new SsLeaf();
    SsNode* newNode2 = new SsLeaf();
    auto splitIndex = this->findSplitIndex();
    dynamic_cast<SsLeaf*>(newNode1)->points = std::vector<Point>(this->points.begin(), this->points.begin() + splitIndex);
    dynamic_cast<SsLeaf*>(newNode2)->points = std::vector<Point>(this->points.begin() + splitIndex, this->points.end());
    newNode1->updateBoundingEnvelope();
    newNode2->updateBoundingEnvelope();
    newNode1->parent = this->parent;
    newNode2->parent = this->parent;
    return std::make_pair(newNode1, newNode2);


//    SsNode* newNode1 = new SsLeaf();
//    SsNode* newNode2 = new SsLeaf();
//
//    int dims = this->centroid.dim();
//    Point center1 = Point(dims);
//    Point center2 = Point(dims);
//    NType d, mx;
//    for(const auto& pointA: this->points){
//        for(const auto& pointB: this->points){
//            d = distance(pointA, pointB);
//            if(d > mx){
//                mx = d;
//                center1 = pointA;
//                center2 = pointB;
//            }
//        }
//    }
//
//    int iters = 0;
//    int szG1, szG2;
//    while(iters != 10){
//////        std::cout << "ITERRRRR\t" << iters << std::endl;
//        //redo kmeans
//        szG1 = 0; szG2 = 0;
//
//        for(auto& dot: this->points){
////            if((iters != 0) ||(dot == center1) || (dot == center2)) continue;
//            NType d1 = distance(center1, dot);
//            NType d2 = distance(center2, dot);
//            if(d1 < d2) {dot.group = 1; szG1++;}
//            else {dot.group = 2; szG2++;}
//        }
//        if(szG1 == 0 || szG2 == 0) break;
//
//        //kmeans for center1 and center2
//        for(int i = 0; i < dims; i++) {
//            center1[i] = 0;
//            center2[i] = 0;
//            for (auto &dot: this->points) {
//                if (dot.group == 1) {
//                    center1[i] += dot[i];
//                } else {
//                    center2[i] += dot[i];
//                }
//            }
//            if(szG1 != 0) center1[i] /= szG1;
//            else std::cout << "?????????????" << std::endl;
//            if(szG2 != 0) center2[i] /= szG2;
//            else std::cout << "?????????????????" << std::endl;
//        }
//        iters++;
//    }
//
////    std::cout << center1 << std::endl;
////    std::cout << center2 << std::endl;
//    auto* nodito1 = dynamic_cast<SsLeaf*>(newNode1);
//    auto* nodito2 = dynamic_cast<SsLeaf*>(newNode2);
//
//    nodito1->points.clear();
//    nodito2->points.clear();
//
//    for(const auto& dot: this->points){
//        if(dot.group == 1) nodito1->points.push_back(dot);
//        else nodito2->points.push_back(dot);
//    }
//    nodito1->updateBoundingEnvelope();
//    nodito2->updateBoundingEnvelope();
//    nodito1->parent = this->parent;
//    nodito2->parent = this->parent;
//
//    return std::make_pair(nodito1, nodito2);
}

void SsInnerNode::sortEntriesByCoordinate(size_t coordinateIndex) {
    std::sort(this->children.begin(),this->children.end(), [coordinateIndex](SsNode* a, SsNode* b) {return a->centroid[coordinateIndex] < b->centroid[coordinateIndex];});
}

void SsInnerNode::draw(sf::RenderWindow &window) const {
    sf::CircleShape circle(radius.getValue()/norm);
    circle.setPosition((centroid[0].getValue() - radius.getValue())/norm, (centroid[1].getValue() - radius.getValue())/norm);
    circle.setFillColor(sf::Color::Transparent);
    circle.setOutlineThickness(2);
    circle.setOutlineColor(sf::Color::Red);
    window.draw(circle);

    for (const auto& child : children) {
        if(child->isLeaf()){
            auto* LeafNode = dynamic_cast<SsLeaf*>(child);
            LeafNode->draw(window);
        }
        else{
            auto* InnerNode = dynamic_cast<SsInnerNode*>(child);
            InnerNode->draw(window);
        }
    }
}

void SsLeaf::sortEntriesByCoordinate(size_t coordinateIndex) {
    std::sort(this->points.begin(),this->points.end(), [coordinateIndex](Point a, Point b) {return a[coordinateIndex] < b[coordinateIndex];});
}

void SsLeaf::draw(sf::RenderWindow &window) const {
    sf::CircleShape circle(radius.getValue()/norm);
    circle.setPosition((centroid[0].getValue() - radius.getValue())/norm, (centroid[1].getValue() - radius.getValue())/norm);
    circle.setFillColor(sf::Color::Transparent);
    circle.setOutlineThickness(2);
    circle.setOutlineColor(sf::Color::Blue);
    window.draw(circle);

    for (const auto& point : points) {
        sf::CircleShape circle(2);
        circle.setPosition(point[0].getValue()/norm, point[1].getValue()/norm);
        circle.setFillColor(sf::Color::Green);
        window.draw(circle);
    }
}

void SsTree::draw(sf::RenderWindow &window) const {
    if(root->isLeaf()){
        auto* LeafNode = dynamic_cast<SsLeaf*>(root);
        LeafNode->draw(window);
    }
    else{
        auto* InnerNode = dynamic_cast<SsInnerNode*>(root);
        InnerNode->draw(window);
    }
}



void SsTree::kNNQuery(const Point &target, size_t k, SsNode *node, NType radius, std::priority_queue<std::pair<NType, Point>, std::vector<std::pair<NType, Point>>, decltype(compare)> &result) const {
    //Caso 1
    bool isItInner = !node->isLeaf();
    // si es nodo interno
    if(isItInner){
        auto* innerNode = dynamic_cast<SsInnerNode*>(node);
        //Caso 1
        for(const auto& childNode: innerNode->children) {
            if (NType(radius.getValue()) + innerNode->radius.getValue() < distance(target, childNode->centroid)) {
                //Do nothing
                continue;
            }

            //Caso 3
            if (NType(radius.getValue()) + distance(target, childNode->centroid) < childNode->radius) {
                //Do nothing
                continue;
            }

            //Caso 5
            if (k == 1) {
                if (distance(target, innerNode->centroid) + innerNode->radius < radius.getValue()) {
                    radius = distance(target, innerNode->centroid) + innerNode->radius;
                }
                this->kNNQuery(target, k, childNode, radius, result);
            }

            if(distance(target, childNode->centroid) + childNode->radius <= radius)
                this->kNNQuery(target, k, childNode, radius, result);
        }
    }
    //si es hoja
    else{
        auto* leafNode = dynamic_cast<SsLeaf*>(node);
        //buscar puntos
        for(const auto &point: leafNode->points) {
            //Caso 2
            if (NType(radius.getValue()) + distance(target, point) < distance(target, leafNode->centroid)){
                // Do nothing
                continue;
            }

            //Caso 4
            if (NType(radius.getValue()) + distance(target, leafNode->centroid) < distance(point, leafNode->centroid)){
                //Do nothing
                continue;
            }
            NType d = distance(target,point);
            if(d < radius.getValue()){
                result.push(std::make_pair(d, point));
                if(result.size() > k) result.pop();
            }

        }
    }

}

std::vector<Point> SsTree::busquedaKNN(const Point& target, size_t k) {
    std::priority_queue<std::pair<NType, Point>, std::vector<std::pair<NType, Point>>, decltype(compare)> result(compare);
    this->kNNQuery(target, k, root, NType::max_value(), result);

    std::vector<Point> ans;
    int sz = result.size();
    std::cout << result.size() << std::endl;
    for(int i = 0; (i < k) && (i < sz); i++){
        ans.push_back(result.top().second);
        result.pop();
    }
    std::reverse(ans.begin(), ans.end());

    return ans;
}

void SsTree::insert(const Point &point, const std::string &path) {
    this->insert(point);
}

void SsLeaf::saveToStream(std::ostream &out) const {
    // Guardar centroid
    auto D = centroid.dim();
    centroid.saveToFile(out, D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char*>(&radius_), sizeof(radius_));

    // Guardar el numero de puntos
    size_t numPoints = points.size();
    out.write(reinterpret_cast<const char*>(&numPoints), sizeof(numPoints));

    // Guardar los puntos
    for (const auto& point : points) {
        point.saveToFile(out, D);
    }

    // Guardar las rutas (paths)
    size_t numPaths = paths.size();
    out.write(reinterpret_cast<const char*>(&numPaths), sizeof(numPaths));
    for (const auto& p : paths) {
        size_t pathLength = p.size();
        out.write(reinterpret_cast<const char*>(&pathLength), sizeof(pathLength));
        out.write(p.c_str(), (long) pathLength);
    }
}

void SsInnerNode::saveToStream(std::ostream &out) const {
    // Guardar centroid
    centroid.saveToFile(out, centroid.dim());

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char*>(&radius_), sizeof(radius_));

    // Guardar si apunta a nodos hoja
    bool pointsToLeafs = children[0]->isLeaf();
    out.write(reinterpret_cast<const char*>(&pointsToLeafs), sizeof(pointsToLeafs));

    // Guardar la cantidad de hijos para saber cuántos nodos leer después
    size_t numChildren = children.size();
    out.write(reinterpret_cast<const char*>(&numChildren), sizeof(numChildren));

    // Guardar los hijos
    for (const auto& child : children) {
        child->saveToStream(out);
    }
}





void SsInnerNode::loadFromStream(std::istream &in) {
    // Leer centroid
    auto D = centroid.dim();
    centroid.readFromFile(in, D);

    // leer el valor del radio
    float radius_ = 0;
    in.read(reinterpret_cast<char*>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // leer si apunta a hojas o nodos internos
    bool pointsToLeaf = false;
    in.read(reinterpret_cast<char*>(&pointsToLeaf), sizeof(pointsToLeaf));

    // leer cantidad de hijos
    size_t numChildren;
    in.read(reinterpret_cast<char*>(&numChildren), sizeof(numChildren));

    // leer hijos
    for (size_t i = 0; i < numChildren; ++i) {
        SsNode* child = pointsToLeaf ? static_cast<SsNode*>(new SsLeaf()) : static_cast<SsNode*>(new SsInnerNode());
        child->loadFromStream(in);
        children.push_back(child);
    }
}

void SsLeaf::loadFromStream(std::istream &in) {
    // Leer centroid
    centroid.readFromFile(in, centroid.dim());

    // Leer radio
    float radius_ = 0;
    in.read(reinterpret_cast<char*>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // Leer numero de puntos
    size_t numPoints;
    in.read(reinterpret_cast<char*>(&numPoints), sizeof(numPoints));

    // Leer puntos
    points.resize(numPoints);
    for (size_t i = 0; i < numPoints; ++i) {
        points[i].readFromFile(in, points[i].dim());
    }

    // Leer rutas (paths)
    size_t numPaths;
    in.read(reinterpret_cast<char*>(&numPaths), sizeof(numPaths));
    paths.resize(numPaths);
    for (size_t i = 0; i < numPaths; ++i) {
        size_t pathLength;
        in.read(reinterpret_cast<char*>(&pathLength), sizeof(pathLength));
        char* buffer = new char[pathLength + 1];
        in.read(buffer, (long) pathLength);
        buffer[pathLength] = '\0';
        paths[i] = std::string(buffer);
        delete[] buffer;
    }
}





void SsTree::saveToFile(const std::string &filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        throw std::runtime_error("Cannot open file for writing");
    }

    // Guardar las dimensiones de la estructura
    auto D = root->centroid.dim();
    out.write(reinterpret_cast<const char*>(&D), sizeof(D));

    // Guardar si el root es hija o nodo interno
    bool isLeaf = root->isLeaf();
    out.write(reinterpret_cast<const char*>(&isLeaf), sizeof(isLeaf));

    // Guardar el resto de la estructura
    root->saveToStream(out);
    out.close();
}

void SsTree::loadFromFile(const std::string &filename, size_t D) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open file for reading");
    }
    if (root) {
        delete root;
        root = nullptr;
    }

    // Aquí se asume que el primer valor determina las dimensiones
    in.read(reinterpret_cast<char*>(&D), sizeof(D));

    // El segundo valor determina si el root es hoja
    bool isLeaf;
    in.read(reinterpret_cast<char*>(&isLeaf), sizeof(isLeaf));
    if (isLeaf) {
        root = new SsLeaf();
    } else {
        root = new SsInnerNode();
    }
    root->loadFromStream(in);
    in.close();
}
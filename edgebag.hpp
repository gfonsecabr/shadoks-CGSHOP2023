#pragma once

#include "tools.hpp"
#include "convex.hpp"

using Cell = std::pair<int,int>;

class EdgeBag {
  tsl::sparse_map<Cell,std::vector<int>> gridEdges;
  std::vector<Segment> edges;
  int cellSize;

public:
  size_t size() {
    return edges.size();
  }

  EdgeBag() {
  }

  EdgeBag(const std::vector<Segment> &_edges) {
    Message M("edgebag", 4, "Building edge bag");
    edges = _edges;
    assert(edges.size() > 0);
    Point boxMin, boxMax;
    boxMin = boxMax = edges[0].source();

    for(const Segment &s : edges) {
      for(Point p : {s.source(), s.target()}) {
        boxMin = Point(std::min(p.x(),boxMin.x()), std::min(p.y(),boxMin.y()));
        boxMax = Point(std::max(p.x(),boxMax.x()), std::max(p.y(),boxMax.y()));
      }
    }
    int n = edges.size();
    int width = CGAL::to_double(std::min(boxMax.x() - boxMin.x(), boxMax.y() - boxMin.y()));
    double refLength = sqrt(CGAL::to_double(edges[n/2].squared_length()));
    // width / cellSize = cellSize / refLength
    cellSize = 1 + sqrt(width * refLength)/4;

    for(int ei = 0; ei < edges.size(); ei++)
      insert(ei);

    M.close(gridEdges.size());
  }

  void write(std::string fn = "bag.json") {
    auto f = std::ofstream(fn);

    f << "{ \"colors\": [" << std::endl;

    for(auto &[c,indices] : gridEdges) {
      f << " [" << std::endl;
      for(int ei : indices) {
        Segment e = edges[ei];
        f << "  [[" << CGAL::to_double(e.source().x()) << ","
                    << CGAL::to_double(e.source().y()) << "],["
                    << CGAL::to_double(e.target().x()) << ","
                    << CGAL::to_double(e.target().y()) << "]],"
                    << std::endl;
      }
      auto corners = boxCorners(c);
      for(int i = 0; i < 4; i++) {
        f << "  [[" << CGAL::to_double(corners[i].x()) << ","
                    << CGAL::to_double(corners[i].y()) << "],["
                    << CGAL::to_double(corners[(i+1)%4].x()) << ","
                    << CGAL::to_double(corners[(i+1)%4].y()) << "]]";
        if(i != 3)
          f << ",";
        else
          f << std::endl;
      }

      f << " ]," << std::endl;
    }

    f << "[]]}" << std::endl;
  }


  void showStats() {
    std::cout << "cellSize: " << cellSize << std::endl;
    std::cout << "cells: " << gridEdges.size() << std::endl;
    std::cout << "n: " << edges.size() << std::endl;
    int stored = 0, maxload = 0;
    for(auto &[_,indices] : gridEdges) {
      stored += indices.size();
      if(indices.size() > maxload)
        maxload = indices.size();
    }
    std::cout << "stored: " << stored << std::endl;
    std::cout << "maxload: " << maxload << std::endl;
  }


  Cell cell(Point p) const {
    return std::make_pair(int(CGAL::to_double((p.x() +.5)/ cellSize)), int(CGAL::to_double((p.y()+.5) / cellSize)));
  }

  std::pair<Point,Point> cellBox(Cell c) const {
    return std::make_pair(Point(c.first * cellSize - .5,     c.second * cellSize - .5),
                          Point((c.first+1) * cellSize - .5, (c.second+1) * cellSize - .5));
  }

  std::array<Point,4> boxCorners(Cell c) const {
    auto box = cellBox(c);
    return {box.first, Point(box.second.x(), box.first.y()),
            box.second, Point(box.first.x(), box.second.y())};
  }

  bool boundaryIntersection(Cell c, const Segment &e) const {
    auto corners = boxCorners(c);
    for(int i = 0; i < 4; i++) {
      Segment s = Segment(corners[i], corners[(i+1)%4]);
      if(CGAL::do_intersect(e, s))
        return true;
    }
    return false;
  }

  std::array<Cell,8> neighboringCells(Cell c) const {
    std::array<Cell,8> ret;
    int i = 0;
    for(int dx = -1; dx <= 1; dx++)
      for(int dy = -1; dy <= 1; dy++)
        if(dx != 0 || dy != 0)
          ret[i++] = std::make_pair(c.first + dx, c.second + dy);

    return ret;
  }

  void nextCell(Cell &c, const Segment &e) const {
    std::pair<Point,Point> box = cellBox(c);
    if(CGAL::do_intersect(e, Segment(Point(box.second.x(), box.first.y()), box.second)))
      c.first++;
    else if(e.target().y() > e.source().y())
      c.second++;
    else
      c.second--;
  }

  std::vector<Cell> cells(const Segment &s) const {
    Segment e(s.source() <= s.target() ? s : s.opposite());
    Cell c0 = cell(e.source());
    Cell c1 = cell(e.target());
    if(c0 == c1)
      return {c0};
    std::vector<Cell> ret = {c0};
    while(ret.back() != c1) {
      Cell c = ret.back();
      nextCell(c, e);
      ret.push_back(c);
    }
    return ret;
  }

  // This function considers cells closed everywhere
  std::vector<Cell> cellsPlus(const Segment &s) const {
    std::vector<Cell> cellsVec = cells(s);
    std::set<Cell> cellsSet(cellsVec.begin(), cellsVec.end());
    for(Cell c : cellsVec) {
      for(Cell c2 : neighboringCells(c)) {
        if(boundaryIntersection(c2, s))
          cellsSet.insert(c2);
      }
    }

    return std::vector<Cell>(cellsSet.begin(), cellsSet.end());
  }

  void insert(int ei) {
    const Segment &e = edges[ei];
    auto l = cellsPlus(e);
    for(Cell c : l)
      gridEdges[c].push_back(ei);
  }

  struct Iterator {
    Segment s;
    Cell c, lastCell;
    size_t i;
    EdgeBag *bag;
    tsl::sparse_set<int> visited;

  public:
    Iterator(EdgeBag *_bag, const Segment &_s) {
      s = _s.source() < _s.target() ? _s : _s.opposite();
      bag = _bag;
      i = -1;
      c = bag->cell(s.source());
      lastCell = bag->cell(s.target());
      operator++();
    }

    // End iterator
    Iterator() : bag(nullptr) {
    }

    Iterator &begin() {
      return *this;
    }

    Iterator end() {
      return Iterator();
    }

    bool tooFar(Cell c) {
      if(s.target().y() > s.source().y())
        return c > lastCell;
      return std::make_pair(c.first,-c.second) > std::make_pair(lastCell.first,-lastCell.second);
    }

    Iterator &operator++() {
      i++;
      auto it = bag->gridEdges.find(c);
      while(it == bag->gridEdges.end() || i >= it->second.size() || visited.contains(curIndex())) {
        // std::cout << "Cell: " << c.first << "," << c.second << std::endl;
        if(it != bag->gridEdges.end() && i < it->second.size() - 1)
          i++;
        else {
          bag->nextCell(c,s);
          if(tooFar(c)) {
            bag = nullptr;
            return *this;
          }
          i = 0;
          it = bag->gridEdges.find(c);
        }
      }
      visited.insert(curIndex());

      return *this;
    }

    int curIndex() {
      std::vector<int> &v = bag->gridEdges.at(c);
      assert(v.size() > i);
      return v[i];
    }

    Segment operator*() {
      return bag->edges[curIndex()];
    }

    bool operator!= (const Iterator& other) {
      return bag != nullptr || other.bag != nullptr;
    }
  };

  Iterator relevantEdges(const Segment &s) {
    return Iterator(this, s);
  }
};



#pragma once
#include "instance.hpp"

std::vector<Point> edgeEdge(Instance &instance) {
  const std::vector<Segment> &edges = instance.edges();
  int es = edges.size();
  std::vector<Point> ret;

  for(int i = 0; i < es; i++) {
    Line li = Line(edges[i]);
    for(int j = 0; j < es; j++) {
      const auto result = intersection(li, edges[j]);
      if(result) {
        const Point* p = boost::get<Point >(&*result);
        if(p && instance.visible(*p, edges[i].source()))
          ret.push_back(*p);
      }
    }
  }

  return ret;
}

std::vector<Point> edgeEdgePlus(Instance &instance) {
  const std::vector<Segment> &edges = instance.edges();
  int es = edges.size();
  std::vector<Point> ret;

  for(int i = 0; i < es; i++) {
    Line li = Line(edges[i]);
    for(int j = 0; j < es; j++) {
      Line lj = Line(edges[j]);
      const auto result = intersection(li, lj);
      if(result) {
        const Point* p = boost::get<Point >(&*result);
        if(p && instance.containsPoint(*p) &&
           instance.visible(*p, edges[i].source()) &&
           instance.visible(*p, edges[j].source()))
          ret.push_back(*p);
      }
    }
  }

  return ret;
}

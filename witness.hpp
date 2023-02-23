#pragma once
#include "instance.hpp"

Arrangement buildArrangement(const std::vector<Convex> &convexes, Instance &instance) {
  Arrangement arr;

  // To identify holes, we insert a point inside each hole
  Face_handle uf = arr.unbounded_face();
  for(Point p : instance.inHole)
    arr.insert_in_face_interior(p, uf);

  std::set<Segment,segLessThan> insertedEdges;

  for(Segment e : instance.edges())
    if(!insertedEdges.contains(e)) {
      CGAL::insert(arr,e);
      insertedEdges.insert(e);
    }

  for(const Convex &poly : convexes)
    for(Segment e : poly.edges())
      if(!insertedEdges.contains(e)) {
        CGAL::insert(arr,e);
        insertedEdges.insert(e);
      }

  return arr;
}


std::vector<Witness> findVertexWitnesses(const std::vector<Convex> &convexes, Instance &instance) {
  Arrangement arr = buildArrangement(convexes, instance);

  std::vector<Witness> ret;

  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    std::vector<Point> v;
    bool found = false;

    if(!fit->is_unbounded() && fit->number_of_holes() == 0 && fit->number_of_isolated_vertices() == 0) {
      auto curr = fit->outer_ccb();
      v.push_back(curr->source()->point());
      do {
        v.push_back(curr->target()->point());
        if(!found && instance.allPointsMap.contains(curr->target()->point()))
          found = true;
      } while (++curr != fit->outer_ccb());

      if(found) {
        Point wpt = Convex(v).pointInside();
        ret.push_back(Witness(wpt,wpt));
      }
    }
  }

  return ret;
}

std::vector<Witness> findWitnesses(const std::vector<Convex> &convexes, Instance &instance) {
  Arrangement arr = buildArrangement(convexes, instance);

  std::vector<Witness> ret;

  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    std::vector<Point> v;

    if(!fit->is_unbounded() && fit->number_of_holes() == 0 && fit->number_of_isolated_vertices() == 0) {
      auto curr = fit->outer_ccb();
      do {
        v.push_back(curr->target()->point());
      } while (++curr != fit->outer_ccb());

      Point wpt = Convex(v).pointInside();
      ret.push_back(Witness(wpt,wpt));
    }
  }

  return ret;
}


std::vector<Witness> findQuickVertexWitnesses(const std::vector<Convex> &convexes, Instance &instance) {
  // Organize instance edges
  std::map<Point,Point> instanceEdgeFrom;
  std::map<Point,std::set<Point>> edgesIncident;
  for(Segment e : instance.edges()) {
    instanceEdgeFrom[e.source()] = e.target();
    edgesIncident[e.source()].insert(e.target());
    edgesIncident[e.target()].insert(e.source());
  }

  // Organize collection edges
  for(const Convex &poly : convexes) {
    for(const Segment &e : poly.edges()) {
      edgesIncident[e.source()].insert(e.target());
      edgesIncident[e.target()].insert(e.source());
    }
  }

  // Build witnesses
  std::vector<Witness> ret;
  for(auto &[p, s] : edgesIncident) {
    if(instanceEdgeFrom.contains(p)) {
      Point q = instanceEdgeFrom[p];
      std::vector<Point> neighbors(s.begin(), s.end());
      std::sort(neighbors.begin(), neighbors.end(),
                [p,q](const Point &a, const Point &b){
                  if(leftOrSameDirection(p,q,a)) {
                    if(leftOrSameDirection(p,q,b)) {
                      return CGAL::orientation(p,a,b) == CGAL::LEFT_TURN;
                    }
                    else {
                      return true;
                    }
                  }
                  if(leftOrSameDirection(p,q,b)) {
                    return false;
                  }
                  return CGAL::orientation(p,a,b) == CGAL::LEFT_TURN;
                });
      std::vector<Point> neighborsTrimmed;
      for(int i = 0; i < neighbors.size(); i++) {
        if(i == 0 ||
          CGAL::orientation(p,neighbors[i-1],neighbors[i]) != CGAL::COLLINEAR
          || (neighbors[i-1] < p && p < neighbors[i])
          || (neighbors[i-1] > p && p > neighbors[i]))
          neighborsTrimmed.push_back(neighbors[i]);
      }
      for(int i = 1; i < neighborsTrimmed.size(); i++) {
        Vector vw;
        if(CGAL::orientation(p, neighborsTrimmed[i-1], neighborsTrimmed[i]) != CGAL::COLLINEAR)
          vw = Vector(p, CGAL::midpoint(neighborsTrimmed[i-1],neighborsTrimmed[i]));
        else if(p != neighborsTrimmed[i])
          vw = Vector(p,neighborsTrimmed[i]).perpendicular(CGAL::CLOCKWISE);

        if(vw.direction().counterclockwise_in_between(
            Vector(p,neighborsTrimmed[i]).direction(),
            Vector(p,neighborsTrimmed[i-1]).direction()))
          vw = - vw;
        ret.push_back(Witness(p, vw));
      }
    }
  }

  return ret;
}

std::vector<Convex> uncoveredFaces(const std::vector<Convex> &polys, Instance &instance) {
  std::vector<Convex> ret;
  Polygon_set coverage;

  coverage.join(polys.begin(), polys.end());

  Polygon_set inst(instance.poly);
  inst.difference(coverage);
  std::vector<Polygon_with_holes> pwhs;
  inst.polygons_with_holes (std::back_inserter (pwhs));
  for(Polygon_with_holes & pwh : pwhs) {
    if(pwh.has_holes() || !pwh.outer_boundary().is_convex()) {
      for(Convex &conv : triangulate(pwh))
        ret.push_back(conv);
    }
    else {
      std::vector<Point> cv(pwh.outer_boundary().begin(), pwh.outer_boundary().end());
      ret.push_back(Convex(cv));
    }
  }

  return ret;
}

void fuseConvexes(std::vector<Convex> &polys, Instance &instance) {
  for(int i = 0; i < polys.size()-1; i++) {
    if(polys[i].size()) {
      for(int j = i+1; j < polys.size(); j++) {
        if(polys[j].size()) {
          if(polys[i][0] == polys[j][0] || instance.visible(polys[i][0],polys[j][0])) {
            std::vector<Point> v(polys[i].begin(), polys[i].end());
            for(Point p : polys[j])
              v.push_back(p);
            Convex conv(v);
            if(instance.contains(conv)) {
              polys[i] = conv;
              polys[j] = Convex();
              i--; // Repeat a value of i
              break;
            }
          }
        }
      }
    }
  }
  std::erase_if(polys, [](Convex &c){return c.size() == 0;});
}


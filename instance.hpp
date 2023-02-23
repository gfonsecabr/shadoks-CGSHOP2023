#pragma once
#include "tools.hpp"
#include "convex.hpp"
#include "edgebag.hpp"

class Instance{
public:
  Polygon_with_holes poly;
  Polygon_set comp;
  tsl::sparse_set<std::pair<int,int>> visibilityMap;
  std::map<std::array<Point,3>,bool> validTriangles;
  int numberOfHoles;
  std::string name;
  std::vector<Segment> polyEdges;
  std::vector<Segment> edgesByLength;
  std::map<Segment, Segment, osegLessThan> previousEdge;
  std::vector<Point> inHole, inHole_y;
  std::vector<Point> extraPoints;
  std::map<Point,int> allPointsMap;
  std::vector<Point> allPointsVector;
  std::string filename;
  Arrangement env;
  VIS *vis = nullptr;
  std::set<Point> isolatedVertices;
  EdgeBag bag;
  Point boxMin, boxMax;
  double triangleLimit = 2e8;

  void buildComp() {
    assert(poly.outer_boundary().is_simple());

    assert(poly.outer_boundary().orientation() == CGAL::POSITIVE);
    comp = Polygon_set(poly);
    comp.complement();
  }

  void buildVisibility() {
    Message Mvis("visibility", 4, "Building visibility map");
    auto & edgs = edges();
    CGAL::insert_non_intersecting_curves(env,edgs.begin(),edgs.end());

    vis = new VIS(env);

    for(Halfedge_const_handle he : env.halfedge_handles()) {
      Point p = he->target()->point();

      // To use only halfedges inside the polygon
      if(previousEdge.contains(Segment(he->source()->point(), p))) {
        // Visibility query
        Arrangement output_arr;
        Face_handle fh = vis->compute_visibility(p, he, output_arr);
        for (auto vit = output_arr.vertices_begin();
             vit != output_arr.vertices_end();
             ++vit) {
          Point q = vit->point();
          if(p < q && allPointsMap.contains(q)) {
            int ip = allPointsMap.at(p);
            int iq = allPointsMap.at(q);
            visibilityMap.insert(std::make_pair(ip,iq));
          }
        }
      }
    }
    Mvis.close(visibilityMap.size());
  }

  void buildPreviousEdge() {
    std::vector<Polygon> polys(poly.holes_begin(), poly.holes_end());
    polys.push_back(poly.outer_boundary());

    for(Polygon &simple : polys) {
      std::vector<Segment> pedges(simple.edges_begin(), simple.edges_end());
      pedges.push_back(pedges[0]);
      for(int i = 1; i < pedges.size(); i++) {
        const Segment &previous = pedges[i-1];
        const Segment &cur = pedges[i];
        previousEdge[cur] = previous;
      }
    }
  }

  void buildInHole() {
    for(auto hole : poly.holes()) {
      Point p = pointInsidePolygon(hole);
      inHole.push_back(p);
    }
    std::sort(inHole.begin(), inHole.end());
    inHole_y = inHole;
    std::sort(inHole_y.begin(), inHole_y.end(),smaller_y);
  }

  void buildAll(bool wantBag) {
    const auto vert = vertices();
    boxMin = boxMax = vert[0];

    for(Point p : vert) {
      boxMin = Point(std::min(p.x(),boxMin.x()), std::min(p.y(),boxMin.y()));
      boxMax = Point(std::max(p.x(),boxMax.x()), std::max(p.y(),boxMax.y()));
      if(!allPointsMap.contains(p)) {
        allPointsMap[p] = allPointsMap.size();
        allPointsVector.push_back(p);
      }
    }
    edges();
    buildPreviousEdge();
    buildComp();
    buildInHole();
    buildVisibility();
    if(wantBag)
      bag = EdgeBag(edgesByLength);
  }

  Instance(std::string fn, bool wantBag = true) {
    filename = fn;
    rapidjson::Document doc = readJson(fn);
    std::vector<Point> outer;
    std::vector<Polygon> holes;

    if(VERBOSE >= 3)
      std::cout << "Reading " << fn << std::endl;

    assert(doc["type"] == "CGSHOP2023_Instance");
    name = doc["name"].GetString();
    Message::name = name;

    Message Mob("outer_boundary",7,"Reading outer boundary");
    outer = jsonPointVec(doc["outer_boundary"]);

    int n = outer.size();
    Polygon outerp(outer.begin(), outer.end());
    assert(outerp.is_simple());
    if(outerp.orientation() == CGAL::NEGATIVE)
      outerp.reverse_orientation();
    Mob.close(outer.size());

    if(doc.HasMember("holes")) {
      Message Mh("holes",7,"Reading holes");
      for(auto &dochole: doc["holes"].GetArray()) {
        std::vector<Point> h = jsonPointVec(dochole);
        Polygon hp(h.begin(),h.end());
        if(hp.orientation() == CGAL::POSITIVE)
          hp.reverse_orientation();

        holes.push_back(hp);
        n += h.size();
      }

      poly = Polygon_with_holes(outerp, holes.begin(), holes.end());
      Mh.close(holes.size());
    }
    else {
      poly = Polygon_with_holes(outerp);
    }

    numberOfHoles = holes.size();
    Message Mi("n",4,"Number of vertices");
    Mi.close(n);
    buildAll(wantBag);
  }

  int registerPoints(std::vector<Point> pts) {
    int count = 0;
    std::map<Point,Vertex_const_handle> handles;

    for(Point p : pts) {
      if(!allPointsMap.contains(p)) {
        count++;

        extraPoints.push_back(p);
        handles[p] = CGAL::insert_point(env, p);
        allPointsMap[p] = allPointsMap.size();
        allPointsVector.push_back(p);
     }
    }

    delete vis;
    vis = new VIS(env);

    for(auto &[p, hp] : handles) {
      Arrangement output_arr;

      if(hp->is_isolated()) {
        isolatedVertices.insert(p);
        auto hf = hp->face();
        vis->compute_visibility(p, hf, output_arr);
      }
      else {
        auto he = hp->incident_halfedges();
        assert(he->target()->point() == p);

        // There are two half edges with taget p
        // We need to find the one that bounds the interior
        if(he->face()->is_unbounded() ||
          he->face()->number_of_holes() != inHole.size())
          ++he;

        vis->compute_visibility(p, he, output_arr);
      }

      for (auto vit = output_arr.vertices_begin();
            vit != output_arr.vertices_end();
            ++vit) {
        Point q = vit->point();
        if(allPointsMap.contains(q)) {
          int ip = allPointsMap.at(p);
          int iq = allPointsMap.at(q);
          if(p < q)
            visibilityMap.insert(std::make_pair(ip,iq));
          else if(p > q)
            visibilityMap.insert(std::make_pair(iq,ip));
        }
      }

      // Visibility between two isolated vertices is computed now
      for(Point p : isolatedVertices)
        for(Point q : isolatedVertices)
          if(p < q && slowVisible(p,q)) {
            int ip = allPointsMap.at(p);
            int iq = allPointsMap.at(q);
            visibilityMap.insert(std::make_pair(ip,iq));
          }

    }

    return count;
  }

  bool containsTriangle(std::array<Point, 3> tv) {
      std::sort(tv.begin(), tv.end());

      if (!visible(tv[0],tv[1]) ||
          !visible(tv[1],tv[2]) ||
          !visible(tv[2],tv[0]))
        return false;

      if(!validTriangles.contains(tv)) {
        if(validTriangles.size() > triangleLimit)
          return noHoleInside(Convex({tv[0],tv[1],tv[2]}));
        validTriangles[tv] = noHoleInside(Convex({tv[0],tv[1],tv[2]}));
      }

      return validTriangles.at(tv);
  }

  bool containsPoint(Point p) {
    return comp.oriented_side(p) != CGAL::ON_POSITIVE_SIDE;
  }

  bool contains(const Convex &c) {
    if(c.size() == 0)
      return false;
    if(c.size() == 1)
      return containsPoint(c[0]);
    if(c.size() == 2) {
      Segment s = *(c.edges_begin());
      return visible(s.source(),s.target());
    }
    if(c.size() == 3) {
      std::array<Point,3> tv;
      std::copy_n(c.begin(), 3, tv.begin());
      if(allPointsMap.contains(tv[0]) &&
         allPointsMap.contains(tv[1]) &&
         allPointsMap.contains(tv[2]))
        return containsTriangle(tv);
    }

    return !comp.do_intersect(c);
    // return false;
  }

  std::vector<Point> holeVertices() {
    std::vector<Point> s;
    for(auto hole : poly.holes())
      for(Point p : hole.vertices())
        s.push_back(p);
    return s;
  }

  std::vector<Point> vertices() {
    std::vector<Point> s = holeVertices();
    for(Point p : poly.outer_boundary().vertices())
      s.push_back(p);
    return s;
  }

  std::vector<Point> allPoints() {
    std::vector<Point> s = vertices();
    for(Point p : extraPoints)
      s.push_back(p);

    return s;
  }

  const std::vector<Segment> &edges() {
    if(polyEdges.empty()) {
      for(auto e : poly.outer_boundary().edges())
        polyEdges.push_back(e);
      for(auto hole : poly.holes())
        for(auto e : hole.edges())
          polyEdges.push_back(e);

      edgesByLength = polyEdges;
      std::sort(edgesByLength.begin(), edgesByLength.end(),
                [](const Segment &a, const Segment &b){return a.squared_length() > b.squared_length();});
    }
    return polyEdges;
  }

  bool slowVisible(Point p, Point q) {
    // std::cout<<"X"<<std::flush;

    if(p == q || p.x() < boxMin.x() || p.y() < boxMin.y()
              || q.x() < boxMin.x() || q.y() < boxMin.y()
              || p.x() > boxMax.x() || p.y() > boxMax.y()
              || q.x() > boxMax.x() || q.y() > boxMax.y())
      return false;

    if(p > q)
      std::swap(p,q);
    Segment query(p,q);

    // Polygon is oriented ccw with holes oriented cw

    std::vector<std::tuple<Segment,int,int,int,int>> relevant;

    if(bag.size()) {
      // Test if query and edge cross at a single point in the interior of both
      for(const auto &edge : bag.relevantEdges(query)) {
      // for(const auto &edge : edgesByLength) {
        int op = CGAL::orientation(edge.source(), edge.target(), p);
        int oq = CGAL::orientation(edge.source(), edge.target(), q);
        int osource = -13, otarget = -13;
        if(op != oq && op != CGAL::COLLINEAR && oq != CGAL::COLLINEAR) {
          osource = CGAL::orientation(p,q,edge.source());
          otarget = CGAL::orientation(p,q,edge.target());
          if(osource != otarget && osource != CGAL::COLLINEAR && otarget != CGAL::COLLINEAR)
            return false; // Most function calls end here
        }
        if(osource == -13 || otarget == -13 ||
          op == CGAL::COLLINEAR || oq == CGAL::COLLINEAR ||
          osource == CGAL::COLLINEAR  || otarget == CGAL::COLLINEAR)
          relevant.push_back(std::make_tuple(edge,op,oq,osource,otarget));
      }
    }
    else {
      // Test if query and edge cross at a single point in the interior of both
      for(const auto &edge : edgesByLength) {
      // for(const auto &edge : edgesByLength) {
        int op = CGAL::orientation(edge.source(), edge.target(), p);
        int oq = CGAL::orientation(edge.source(), edge.target(), q);
        int osource = -13, otarget = -13;
        if(op != oq && op != CGAL::COLLINEAR && oq != CGAL::COLLINEAR) {
          osource = CGAL::orientation(p,q,edge.source());
          otarget = CGAL::orientation(p,q,edge.target());
          if(osource != otarget && osource != CGAL::COLLINEAR && otarget != CGAL::COLLINEAR)
            return false; // Most function calls end here
        }
        if(osource == -13 || otarget == -13 ||
          op == CGAL::COLLINEAR || oq == CGAL::COLLINEAR ||
          osource == CGAL::COLLINEAR  || otarget == CGAL::COLLINEAR)
          relevant.push_back(std::make_tuple(edge,op,oq,osource,otarget));
      }
    }

    // Now there is no proper crossing
    // Test if query has an endpoint on edge interior and another enpoint is to the right
    for(auto &[edge,op,oq, osource, otarget] : relevant) {
      if(osource == -13)
        osource = CGAL::orientation(p,q,edge.source());
      if(otarget == -13)
        otarget = CGAL::orientation(p,q,edge.target());
      if(osource != otarget && osource != CGAL::COLLINEAR && otarget != CGAL::COLLINEAR) {
        // edge intersects segment at a point in the interior of edge
        if(op == CGAL::COLLINEAR && oq == CGAL::RIGHT_TURN) return false;
        if(oq == CGAL::COLLINEAR && op == CGAL::RIGHT_TURN) return false;
      }
    }

    // We now only need to handle cases where the intersection happens on an endpoint of edge

    for(auto &[edge,op,oq, osource, otarget] : relevant) {
      if(CGAL::do_intersect(edge, query)) {
        Segment previous = previousEdge.at(edge);
        assert(previous.target() == edge.source());
        if(p == edge.source()) {
          if(rightOfAngle(q, previous.source(), edge.source(), edge.target()))
            return false;
        }
        else if(q == edge.source()) {
          if(rightOfAngle(p, previous.source(), edge.source(), edge.target()))
            return false;
        }
        else if(query.has_on(edge.source())){
          if(rightOfAngle(p, previous.source(), edge.source(), edge.target()) ||
             rightOfAngle(q, previous.source(), edge.source(), edge.target()))
            return false;
        }
      }
    }

    return true;
  }

  bool noHoleInside(const Convex &c) {
    // I know all edges are visible already

    if(c.boxMax.x() - c.boxMin.x() <= c.boxMax.y() - c.boxMin.y()) {
      // Only check points with the x-coordinate in the bbox range
      for(auto it = std::lower_bound(inHole.begin(), inHole.end(), c.boxMin);
          it != inHole.end() && *it <= c.boxMax;
          ++it) {
        Point p = *it;
        if(c.contains(p))
          return false;
      }
      return true;
    }

    // Only check points with the y-coordinate in the bbox range
    for(auto it = std::lower_bound(inHole_y.begin(), inHole_y.end(), c.boxMin, smaller_y);
        it != inHole_y.end() && smallerOrEqual_y(*it, c.boxMax);
        ++it) {
      Point p = *it;
      if(c.contains(p))
        return false;
    }
    return true;

//     return !comp.do_intersect(c);
  }

  bool visible(Point p, Point q) {
    auto itp = allPointsMap.find(p);
    auto itq = allPointsMap.find(q);
    if(itp == allPointsMap.end() || itq == allPointsMap.end())
      return slowVisible(p,q);
    int ip = itp->second;
    int iq = itq->second;
    if(p > q)
      std::swap(ip,iq);

    return visibilityMap.contains(std::make_pair(ip,iq));
  }

  std::string getName() const {
    return name;
  }

  std::string getDataName() const {
    int pos = filename.find_last_of("/");
    pos = filename.find(".", pos);
    std::string datafn = filename.substr(0, pos);
    datafn += ".vis.json";
    return datafn;
  }

  ~Instance() {
    delete vis;
  }
};


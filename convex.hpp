#pragma once

#include "tools.hpp"


class Convex : public Polygon{
  std::vector<Point> sortedVertices;

public:
  Point boxMin, boxMax;

  Convex(){
  }

  Convex(const std::vector<Point> &pts) {
    CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(*this));
    for(Point p : *this)
      sortedVertices.push_back(p);

    boxMin = boxMax = sortedVertices[0];
    for(Point p : *this) {
      boxMin = Point(std::min(p.x(),boxMin.x()), std::min(p.y(),boxMin.y()));
      boxMax = Point(std::max(p.x(),boxMax.x()), std::max(p.y(),boxMax.y()));
    }

    std::sort(sortedVertices.begin(), sortedVertices.end());
  }

  // O(n) but O(log n) is possible
  std::pair<Point,Point> tangents(Point q) const {
    assert(size() >= 2);
    std::vector<Point> retv;

    if(size() == 2) {
      retv = vertices();
    }
    else {
      std::vector<Segment> vedges;
      for(Segment e : edges())
        vedges.push_back(e);

      Segment last = vedges.back();

      int previous = CGAL::orientation(q, last.source(), last.target());

      for(Segment e : vedges) {
        int cur = CGAL::orientation(q, e.source(), e.target());
        if(cur != previous && cur != 0)
          retv.push_back(e.source());

        previous = cur;
        if(retv.size() == 2)
          break;
      }
    }

    assert(retv.size() == 2);
    return std::make_pair(retv[0],retv[1]);
  }

  // O(n) but O(log n) is possible
  bool contains(Point p) const {
    if(size() <= 1)
      return false;

    if(p.x() < boxMin.x() || p.y() < boxMin.y() ||
       p.x() > boxMax.x() || p.y() > boxMax.y())
      return false;

    if(size() == 2) {
      std::vector<Point> vec(this->begin(), this->end());

      return Segment(vec[0],vec[1]).has_on(p);
    }

    for(int i = 1; i+1 < size(); i++) {
      Triangle tri((*this)[0],(*this)[i],(*this)[i+1]);
      if(tri.oriented_side(p) != CGAL::NEGATIVE)
        return true;
    }
    return false;
  }

  bool contains(Witness c) const {
    Point p = c.source();
    bool containsSource = contains(p);

    if(c.is_degenerate() || !containsSource)
      return containsSource;

    int i;
    for(i = 0; i < this->size(); i++) {
      if((*this)[(i+1) % this->size()] != p &&
         CGAL::do_intersect(p, Segment((*this)[i], (*this)[(i+1) % this->size()])))
        break;
    }

    if(i == this->size())
      return true; // p is contained in the interior

    int previousIndex = i;
    if((*this)[previousIndex] == p)
      previousIndex = i != 0 ?  i-1 : this->size() - 1;

    int nextIndex = (i+1) % this->size();
    Point previous = (*this)[previousIndex];
    Point next = (*this)[nextIndex];

    return CGAL::orientation(previous, p, c.point(1)) == CGAL::LEFT_TURN &&
           CGAL::orientation(p, next, c.point(1)) == CGAL::LEFT_TURN;
  }

  bool contains(const Convex &c) const {
    std::set<Point> sthis(begin(),end());
    std::set<Point> sc(c.begin(),c.end());
    std::vector<Point> both(begin(),end());
    for(Point p : sc)
      both.push_back(p);
    Convex c2(both);
    std::set<Point> sb(c2.begin(),c2.end());

    return sb==sthis;
  }

  Point pointInside() {
    if(size() == 1)
      return (*this)[0];
    if(size() == 2)
      return CGAL::midpoint((*this)[0], (*this)[1]);
    if(size() == 3)
      return CGAL::midpoint(CGAL::midpoint((*this)[0], (*this)[1]), (*this)[2]); // Avoid diving by 3
    return CGAL::midpoint((*this)[0],(*this)[size() / 2]);
  }

  auto operator<=>(const Convex& other) const {
    return sortedVertices <=> other.sortedVertices;
  }
};

std::vector<Convex> triangulate(const Polygon_with_holes &poly) {
  CDT cdt;
  cdt.insert_constraint(poly.outer_boundary().vertices_begin(),
                        poly.outer_boundary().vertices_end(), true);
  for(const auto &hole : poly.holes())
    cdt.insert_constraint(hole.vertices_begin(), hole.vertices_end(), true);

  std::unordered_map<CDTFace_handle, bool> in_domain_map;
  boost::associative_property_map< std::unordered_map<CDTFace_handle,bool> > in_domain(in_domain_map);
  // Mark facets that are inside the domain bounded by the polygon
  // This function requires non-stable CGAL from git
  CGAL::mark_domain_in_triangulation(cdt, in_domain);

  std::vector<Convex> ret;

  for (CDTFace_handle f : cdt.finite_face_handles())
  {
    if (get(in_domain, f)) {
      Point p0 = Point(f->vertex(0)->point().x(), f->vertex(0)->point().y());
      Point p1 = Point(f->vertex(1)->point().x(), f->vertex(1)->point().y());
      Point p2 = Point(f->vertex(2)->point().x(), f->vertex(2)->point().y());
      ret.push_back(Convex({p0,p1,p2}));
    }
  }

  return ret;
}

Point pointInsidePolygon(Polygon poly) {
  if(poly.orientation() < 0)
    poly.reverse_orientation();

  Convex convex;
  if(poly.size() <= 2 || poly.is_convex()) {
    std::vector<Point> cv;
    for(Point p : poly)
      cv.push_back(p);
    convex = Convex(cv);
  }
  else {
    convex = triangulate(Polygon_with_holes(poly))[0];
  }

  return convex.pointInside();
}



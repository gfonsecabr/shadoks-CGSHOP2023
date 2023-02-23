#pragma once
#include "instance.hpp"

bool validAddition(Instance &instance,
                   const Convex &conv,
                   const Point &v,
                   bool insideResult = true
                   ) {
  if(conv.size() == 0)
    return true;
  if(conv.size() == 1)
    return instance.visible(v, conv[0]);

  if(conv.size() == 2) {
    return instance.containsTriangle({v, conv[0], conv[1]});
  }

  if(conv.contains(v))
    return insideResult; // Valid but not useful
  auto [t1,t2] = conv.tangents(v);
  if(!instance.containsTriangle({v,t1,t2}))
    return false;

  return true;
}

void maximalConvex(Instance &instance,
                   std::vector<Convex> &ret,
                   const std::vector<Point> &cur,
                   std::set<Point> &candidates,
                   std::set<Point> &forbidden,
                   std::map<Point,int> &deg) {
  if(candidates.empty() && forbidden.empty()) {
    Convex c(cur);
    if(c.size()>=3)
      ret.push_back(c);
    if(VERBOSE > 7)
      std::cout << " " << cur.size() << std::flush;
  }

  if(cur.size() >= 3) {
    Convex c(cur);
    for(Point p : forbidden)
      if(c.contains(p))
        return; // Early termination for faster calculation with many colinear points
  }

  std::vector<Point> vert(candidates.begin(), candidates.end());
  // Sorting improves speed
  std::sort(vert.begin(), vert.end(), [&deg](Point a, Point b){return deg[a] > deg[b];});

  for(Point v : vert) {
    std::set<Point> candidates2(candidates), forbidden2(forbidden);
    std::vector<Point> cur2(cur);
    bool contaisForbidden = false;

    cur2.push_back(v);
    Convex conv2(cur2);

    for(Point p : candidates) {
      if(v==p || !validAddition(instance,conv2, p))
        candidates2.erase(p);
    }

    for(Point p : forbidden) {
      if(v==p || !validAddition(instance,conv2, p))
        forbidden2.erase(p);
    }

    maximalConvex(instance, ret, cur2, candidates2, forbidden2, deg);

    candidates.erase(v);
    forbidden.insert(v);
  }
}

// https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
std::vector<Convex> maximalConvex(Instance &instance) {
  std::set<Point> candidates;
  for(Point p : instance.allPoints())
    candidates.insert(p);

  std::set<Point> forbidden;
  std::vector<Point> cur;

  std::map<Point,int> deg;
  for(auto [i1,i2] : instance.visibilityMap) {
    Point p1 = instance.allPointsVector[i1];
    Point p2 = instance.allPointsVector[i2];
    deg[p1]++;
    deg[p2]++;
  }

  std::vector<Convex> ret;

  maximalConvex(instance, ret, cur, candidates, forbidden, deg);

  return ret;
}

std::vector<Point> rgrownConvex(Instance &instance, const Convex &convex, std::set<Point> candidates = {}) {
  std::vector<Point> cur(convex.begin(), convex.end());

  if(candidates.empty())
    for(Point p : instance.allPoints())
      if(instance.visible(cur[0], p))
        candidates.insert(p);

  for(Point p:cur)
    candidates.erase(p);

  Convex conv(cur);
  while(!candidates.empty()) {

    std::vector<Point> candidatesv(candidates.begin(), candidates.end());
    Point r = candidatesv[rand() % candidatesv.size()];
    if(validAddition(instance, conv, r, false)) {
      cur.push_back(r);
      conv = Convex(cur);
    }
    candidates.erase(r);
  }
  std::sort(cur.begin(), cur.end());

  return cur;
}

std::vector<Convex> rgrownTriangles(Instance &instance, int rep) {
  std::vector<Convex> triangulation = triangulate(instance.poly);
  std::set<std::vector<Point>> convset;

  for(const Convex &conv : triangulation) {
    if(VERBOSE > 5)
      std::cout<<'.'<<std::flush;
    for(int i = 0; i < rep; i++)
      convset.insert(rgrownConvex(instance, conv));
  }

  std::vector<Convex> ret;
  for(std::vector<Point> pts : convset)
    ret.push_back(Convex(pts));

  return ret;
}


std::set<Point> bloatPoints2(Instance &instance, Convex &conv) {
  const std::vector<Segment> &edges = instance.edges();

  std::set<Point> ret;

  std::vector<Segment> cedges(conv.edges_begin(), conv.edges_end());
  for(Segment s : cedges) {
    Line li = Line(s);
    std::vector<Point> pline = {s.source()};
    for(int j = 0; j < edges.size() + cedges.size(); j++) {
      Segment e;
      if(j < edges.size())
        e = edges[j];
      else
        e = cedges[j-edges.size()];
      const auto &result = intersection(li, Line(e));

      if(result) {
        const Point* p = boost::get<Point >(&*result);
        if(p && *p != s.source())
          pline.push_back(*p);
      }
    }
    std::sort(pline.begin(),pline.end());
    for(auto it = std::lower_bound(pline.begin(), pline.end(), s.source());
        it != pline.end();
        ++it) {
      if(*it != s.source()) {
        if(*it != s.source() && instance.visible(*it, s.source())) {
          ret.insert(*it);
        }
        else
          break;
      }
    }

    for(auto it = std::lower_bound(pline.begin(), pline.end(), s.source());
        it >= pline.begin();
        --it) {
      if(*it != s.source()) {
        if(*it != s.source() && instance.visible(*it, s.source())) {
          ret.insert(*it);
        }
        else
          break;
      }
    }
  }

  return ret;
}

std::set<Point> bloatPoints1(Instance &instance, Convex conv) {
  const std::vector<Segment> &edges = instance.edges();
  std::set<Point> ret;

  for(Segment s : conv.edges()) {
    Line li = Line(s);
    std::vector<Point> pline = {s.source()};
    for(int j = 0; j < edges.size(); j++) {
      const auto &result = intersection(li, edges[j]);

      if(result) {
        const Point* p = boost::get<Point >(&*result);
        if(p)
          pline.push_back(*p);
      }
    }
    std::sort(pline.begin(),pline.end());
    for(auto it = std::lower_bound(pline.begin(), pline.end(), s.source());
        it != pline.end();
        ++it) {
      if(*it != s.source()) {
        if(*it != s.source() && instance.visible(*it, s.source())) {
          ret.insert(*it);
        }
        else
          break;
      }
    }

    for(auto it = std::lower_bound(pline.begin(), pline.end(), s.source());
        it >= pline.begin();
        --it) {
      if(*it != s.source()) {
        if(*it != s.source() && instance.visible(*it, s.source())) {
          ret.insert(*it);
        }
        else
          break;
      }
    }
  }

  return ret;
}


std::set<Point> bloatPoints(Instance &instance, Convex &conv, int bloat) {
  std::set<Point> retset;

  if(bloat == 1)
    retset = bloatPoints1(instance, conv);
  else
    retset = bloatPoints2(instance, conv);

  return retset;
}

std::vector<Convex> bloatPolys(Instance &instance, std::vector<Convex> polys, int bloat) {
  Message Mb("bloat",4,"Bloating convex polygons");
  std::vector<Convex> ret;
  for(int i = 0 ; i < polys.size(); i++) {
    std::set<Point> bpts = bloatPoints(instance, polys[i], bloat);
    ret.push_back(rgrownConvex(instance, polys[i], bpts));
  }
  Mb.close(ret.size());

  return ret;
}

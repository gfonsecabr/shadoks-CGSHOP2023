#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <deque>
#include <utility> //pair
#include <algorithm>
#include <string>
#include <unistd.h>
#include <chrono>
#include <sstream>
#include <cmath>
#include <bits/stdc++.h>
#include "rapidjson/document.h"
#include "rapidjson/istreamwrapper.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/ostreamwrapper.h"
#include <ilcplex/ilocplex.h>
#include "tsl/sparse_set.h"
#include "tsl/sparse_map.h"
#include <boost/unordered_set.hpp> // hash_combine

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Polygon_with_holes = CGAL::Polygon_with_holes_2<K>;
using Polygon_set = CGAL::Polygon_set_2<K>;
using Point = CGAL::Point_2<K>;
using Witness = CGAL::Ray_2<K>;
using Vector = CGAL::Vector_2<K>;

using Segment = K::Segment_2;
using RT = K::RT;
using Polygon = CGAL::Polygon_2<K>;
using Line = K::Line_2;

using Triangle = CGAL::Triangle_2<K>;

using Traits = CGAL::Arr_segment_traits_2<K>;
using Arrangement = CGAL::Arrangement_2<Traits>;
using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
using Vertex_const_handle = Arrangement::Vertex_const_handle;
using Face_handle = Arrangement::Face_handle;
using VIS = CGAL::Triangular_expansion_visibility_2<Arrangement,CGAL::Tag_true>;

using Vb = CGAL::Triangulation_vertex_base_2<K>;
using Fb = CGAL::Constrained_triangulation_face_base_2<K>;
using TDS = CGAL::Triangulation_data_structure_2<Vb,Fb>;
using Itag = CGAL::Exact_predicates_tag;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag>;
using CDTFace_handle = CDT::Face_handle;


int VERBOSE = 5;
std::map<std::string, std::string> meta;
auto beginTime = std::chrono::high_resolution_clock::now();

double elapsed() {
  auto end = std::chrono::high_resolution_clock::now();
  auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(end - beginTime);

  return dur.count() / 1000.0;
}

int myrandom (int i) { return rand()%i;}

std::string doubleString(double x) {
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(3) << x;
  return ss.str();
}

struct segLessThan {
  bool operator()(const Segment &a, const Segment &b) const {
    std::pair<Point,Point> pa, pb;

    if(a.source() < a.target())
      pa = std::make_pair(a.source(), a.target());
    else
      pa = std::make_pair(a.target(), a.source());

    if(b.source() < b.target())
      pb = std::make_pair(b.source(), b.target());
    else
      pb = std::make_pair(b.target(), b.source());

    return pa < pb;
  }
};

struct osegLessThan {
  bool operator()(const Segment &a, const Segment &b) const {
    std::pair<Point,Point> pa, pb;

    pa = std::make_pair(a.source(), a.target());
    pb = std::make_pair(b.source(), b.target());

    return pa < pb;
  }
};

// Checks if a segment qb is to the right of the polyline abc
bool rightOfAngle(const Point &q, const Point &a, const Point &b, const Point &c) {
  if(q == a || q == b || q == c)
    return false;

  if(CGAL::orientation(a,b,c) == CGAL::LEFT_TURN) {
    return CGAL::orientation(a,b,q) == CGAL::RIGHT_TURN ||
           CGAL::orientation(b,c,q) == CGAL::RIGHT_TURN;
  }
  return CGAL::orientation(a,b,q) == CGAL::RIGHT_TURN &&
         CGAL::orientation(b,c,q) == CGAL::RIGHT_TURN;
}

std::vector<Point> jsonPointVec(rapidjson::Value &values) {
  std::vector<Point> ret;
  for (auto &a : values.GetArray()) {
    std::stringstream ss;

    if(a["x"].IsObject()) {
      std::string num = a["x"]["num"].GetString();
      std::string den = a["x"]["den"].GetString();
      ss << num << "/" << den;
    }
    else {
      std::string num = a["x"].GetString();
      ss << num;
    }

    ss << " ";

    if(a["y"].IsObject()) {
      std::string num = a["y"]["num"].GetString();
      std::string den = a["y"]["den"].GetString();
      ss << num << "/" << den;
    }
    else {
      std::string num = a["y"].GetString();
      ss << num;
    }

    CGAL::IO::set_ascii_mode(ss);
    Point p;
    ss >> p;

    ret.push_back(p);
  }

  return ret;
}

rapidjson::Document readJson(std::ifstream &in) {
  rapidjson::IStreamWrapper isw {in};
  rapidjson::Document doc {};
  doc.ParseStream<rapidjson::kParseNumbersAsStringsFlag>(isw);

  if (doc.HasParseError()) {
    std::cerr << "Error  : " << doc.GetParseError()  << std::endl;
    std::cerr << "Offset : " << doc.GetErrorOffset() << std::endl;
    exit(EXIT_FAILURE);
  }

  return doc;
}

rapidjson::Document readJson(std::string filename) {
  std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
  if (!in.is_open()){
    std::cerr << "Error reading " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  return readJson(in);
}

void parsePoint(std::string s,
                std::string &xnums,
                std::string &xdens,
                std::string &ynums,
                std::string &ydens) {
  std::string xs = s.substr(0, s.find(' '));
  std::string ys = s.substr(s.find(' ')+1);

  if(xs.find('/') == std::string::npos) {
    xnums = xs;
    xdens = "1";
  }
  else{
    xnums = xs.substr(0, xs.find('/'));
    xdens = xs.substr(xs.find('/')+1);
  }
  if(ys.find('/') == std::string::npos) {
    ynums = ys;
    ydens = "1";
  }
  else{
    ynums = ys.substr(0, ys.find('/'));
    ydens = ys.substr(ys.find('/')+1);
  }
}

std::string fracStr(std::string num, std::string den) {
  std::ostringstream ss;
  if(den == "1")
    ss << num;
  else
    ss << "{\"num\":" << num
       << ",\"den\":" << den
       << "}";
  return ss.str();
}

class Message {
  int vlevel;
  double t0;
  std::string id;
public:
  inline static std::string name = "";
  Message(const std::string &_id, int _vlevel = 0, const std::string &m = "") {
    t0 = elapsed();
    id = _id;
    vlevel = _vlevel;
    if(vlevel <= VERBOSE) {
      std::cout << name << "> ";
      if(m.empty())
        std::cout << id << std::flush;
      else
        std::cout << m << std::flush;
    }
  }

  template<class T>
  void close(const T &m) {
    std::ostringstream ss;
    ss << m;
    double t = elapsed() - t0;
    if(vlevel <= VERBOSE) {
      if(t > 1.0)
        std::cout << " (" << int(t) << "s)";
      std::cout << ": " << m << std::endl;
    }
    if(meta.contains(id))
      meta[id] += "; ";
    meta[id] += ss.str();
    if(meta.contains(id + "_sec"))
      meta[id + "_sec"] += "; ";
    meta[id + "_sec"] += doubleString(t);
  }
};

namespace std {
  template <> struct hash<std::pair<int,int>> {
    size_t operator()(std::pair<int,int> p) const {
      size_t seed = 0;
      boost::hash_combine(seed, p.first);
      boost::hash_combine(seed, p.second);
      return seed;
    }
  };
}

bool smaller_y(const Point &p, const Point &q) {
  return std::make_tuple(p.y(),p.x()) < std::make_tuple(q.y(),q.x());
}

bool smallerOrEqual_y(const Point &p, const Point &q) {
  return std::make_tuple(p.y(),p.x()) <= std::make_tuple(q.y(),q.x());
}

bool leftOrSameDirection(const Point &p, const Point &q, const Point &r) {
  int o = CGAL::orientation(p,q,r);
  if(o == CGAL::LEFT_TURN)
    return true;
  if(o == CGAL::RIGHT_TURN)
    return false;

  return (p < q && p < r) || (p > q && p > r);
}

bool witPointLessThan(const Witness &a, const Witness &b) {
  return a.source() < b.source();
}

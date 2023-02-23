#pragma once

#include "convex.hpp"
#include "tools.hpp"

struct Solution {
  std::vector<Convex> polys;
  std::string instance_name;
  std::string algorithm;

  Solution() {
    meta["start_time"] = timeString();
    meta["author"] = "gdf";
    char hn[80];
    gethostname(hn, 80);
    meta["host"] = std::string(hn);
  }

  Solution(std::string fn, int quiet = 0) {
    Message M("readsol",4 + quiet*1000,"Reading " + fn);
    rapidjson::Document doc = readJson(fn);
    assert(doc["type"] == "CGSHOP2023_Solution");
    instance_name = doc["instance"].GetString();

    for(auto &docconv: doc["polygons"].GetArray()) {
      std::vector<Point> v = jsonPointVec(docconv);
      if(v.size() >= 3) {
        Convex c(v);
        polys.push_back(c);
      }
      else {
        std::cerr << "Solution " << fn << " has a polygon with " << v.size() << "vertices!!!" << std::endl;
      }
    }
    M.close(polys.size());
    if(quiet == 1)
      std::cout  << " " << polys.size();
  }

  std::string writeSol(bool quiet = false, std::string filename = "", std::string type = "") const {
    meta["end_time"] = timeString();
    meta["elapsed_sec"] =  doubleString(elapsed());
    meta["algorithm"] = algorithm;
    int pos = instance_name.find(".");
    std::string instance_short = instance_name.substr(0, pos);
    filename += instance_short + "." + algorithm + "." + timeString();
    if(!type.empty())
      filename = filename + "." + type;
    filename += ".sol.json";

    if (!quiet)
      std::cout << "=" << polys.size() << "=> " << filename << std::endl;

    std::ofstream file(filename, std::fstream::out | std::ifstream::binary);

    file << "{" << std::endl;
    file << "\t\"type\": \"CGSHOP2023_Solution\"," << std::endl;
    file << "\t\"instance\": \"" << instance_name << "\"," << std::endl;
    file << "\t\"num_polygons\": " << polys.size() << "," << std::endl;
    file << "\t\"meta\": {" << std::endl;
    int i = 0;
    for(auto const& [key, val] : meta) {
      file << "\t\t\"" << key << "\": \"" << val << "\"";
      if(i < meta.size() - 1)
        file << ",";
      file << std::endl;
      i++;
    }
    file << "\t}," << std::endl;

    file << "\t\"polygons\": ["<< std::endl;
    int polyc = 0;
    for(const Convex &c : polys) {
      file << "\t\t[";
      int pointc = 0;

      for(Point p : c) {
        auto pe = p.exact();
        std::stringstream sfrac;
        sfrac << pe;
        std::string xnum, xden, ynum, yden;
        parsePoint(sfrac.str(), xnum, xden, ynum, yden);

        file << "{\"x\":"
             << fracStr(xnum,xden)
             << ", \"y\":"
             << fracStr(ynum,yden)
             << "}";

        if (++pointc < c.size())
          file << ", ";
      }

      file << "]";
      if (++polyc < polys.size())
        file << ",";
      file<<std::endl;
    }
    file << "\t]" << std::endl;
    file << "}" << std::endl;
    file.close();
    return filename;
  }


    std::string timeString(std::chrono::system_clock::time_point tp = std::chrono::system_clock::now()) const
    {
        std::time_t tt = std::chrono::system_clock::to_time_t(tp);
        struct std::tm * ptm = std::localtime(&tt);
        std::ostringstream oss;
        oss << std::put_time(ptm, "%Y%m%d-%H%M%S");
        std::string s = oss.str();
        return s;
    }
};


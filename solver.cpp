#include "cxxopts.hpp"
#include "tools.hpp"
#include "instance.hpp"
#include "maximal.hpp"
#include "witness.hpp"
#include "solution.hpp"
#include "extrapoints.hpp"
#include "setcover.hpp"

cxxopts::Options options("./solver", "CG:SHOP 2023 Cover a polygon with convex polygons");
cxxopts::ParseResult par;

void parse(int argc, char **argv) {
  options.add_options("Misc")
  ("help", "Print help.")
  ("s,seed", "Random seed, negative for automatic.", cxxopts::value<int>()->default_value("-1"))
  ("v,verbose", "Verbose level.", cxxopts::value<int>()->default_value("5"))
  ("bag", "Build edge bag.", cxxopts::value<bool>()->default_value("true"))
  ("triangle", "Maximum number of cached triangles.", cxxopts::value<double>()->default_value("2e7"))
  ("O,outpath", "Path to save solution.", cxxopts::value<std::string>())
  ("t,savetmp", "Save intermediate solutions.", cxxopts::value<bool>()->default_value("true"))
  ;
  options.add_options("Set Cover")
  ("a,algorithm", "Algorithm name. Possible values are auto, greedy, annealing and mip", cxxopts::value<std::string>()->default_value("auto"))
  ("f,fix", "Fix (at most k) uncovered regions manually.", cxxopts::value<int>()->default_value("0"))
  ("iter", "Number of iterations  (only annealing).", cxxopts::value<int>()->default_value("20000"))
  ("repeat", "Number of repetitions (only annealing).", cxxopts::value<int>()->default_value("6"))
  ("del", "Number of sets to delete at each iteration  (only annealing).", cxxopts::value<int>()->default_value("3"))
  ("to", "Time out for an iteration of mip and annealing (sec).", cxxopts::value<int>()->default_value("600"))
  ("w,witnesses", "Type of witnesses: (1) all (not to be used), (2) vertex, or (3) quick vertex.", cxxopts::value<int>()->default_value("3"))
  ;
  options.add_options("Collection")
  ("c,collection", "Set of convex polygons to choose from. Can be a solution file, 'max' for all maximal sets, 'rgrow' for randomly bloated triangles, 'merge' to merge files given as positional arguments.", cxxopts::value<std::string>()->default_value("max"))
  ("p,ipts", "Insert intermediate points from (1) edges or (2) edges including inside points.", cxxopts::value<int>()->default_value("0"))
  ("ngrow", "Number of convex polygons to produce for each triangle in a triangulation (only for rgrow).", cxxopts::value<int>()->default_value("4"))
  ("savecol", "Save collection file before running algorithm.", cxxopts::value<bool>()->default_value("false"))
  ("C,colpath", "Path to save/load collections.", cxxopts::value<std::string>())
  ("b,bloat", "Post-process the collection adding extra points. (0) none (fast), (1) edges or (2) edges including inside points (slow).", cxxopts::value<int>()->default_value("0"))
  ("collections", "Names of collections when 'merge' is given", cxxopts::value<std::vector<std::string>>());
  ;
  options.add_options("Instance")
  ("i,instance", "Instance file name (required unless a collection file is given).", cxxopts::value<std::string>())
  ("I,ipath", "Path of the instance files.", cxxopts::value<std::string>())
  ;

  options.parse_positional({"collections"});
  par = options.parse(argc, argv);
}

std::string baseName(std::string fn) {
    int i = fn.find_last_of('/');
    if(i == std::string::npos)
      i = -1;
    int j = fn.find('.', i+1);
    std::string basename = fn.substr(i+1,j-i-1);
    return basename;
}

std::string colName(std::string fn) {
    int i = fn.find_last_of('/');
    if(i == std::string::npos)
      i = -1;
    int j = fn.find('.', i+1);
    int k = fn.find('.', j+1);
    std::string colname = fn.substr(j+1,k-j-1);
    return colname;
}

int main(int argc, char **argv) {
  //-----------------------------------------------------------------
  // Read parameters
  //-----------------------------------------------------------------

  for(int i = 0; i < argc; i++) {
    meta["command"] += argv[i];
    if(i != argc-1)
      meta["command"] += " ";
  }
  std::cout << meta["command"] << std::endl;

  parse(argc, argv);
  int seed = par["seed"].as<int>();
  if(seed < 0)
    seed = static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()) % 1000000;
  srand(seed);
  meta["seed"] = std::to_string(seed);
  VERBOSE = par["verbose"].as<int>();

  if (par.count("help")) {
    std::cout << options.help() << std::endl;
    return 1;
  }

  std::string instpath;
  if (par.count("ipath")) {
    instpath = par["ipath"].as<std::string>();
    if(instpath.back() != '/')
      instpath.push_back('/');
  }

  std::string colpath;
  if (par.count("colpath")) {
    colpath = par["colpath"].as<std::string>();
    if(colpath.back() != '/')
      colpath.push_back('/');
  }

  std::string instfn;
  if (par.count("instance"))
    instfn = instpath + par["instance"].as<std::string>();

  int ipts = par["ipts"].as<int>();
  int bloat = par["bloat"].as<int>();
  int wit = par["witnesses"].as<int>();
  bool savecol = par["savecol"].as<bool>();
  bool savetmp = par["savetmp"].as<bool>();

  std::string colalg;
  std::vector<std::string> colfns;
  if (par["collection"].as<std::string>() == "max" ||
      par["collection"].as<std::string>() == "rgrow" ||
      par["collection"].as<std::string>() == "quit" ||
      par["collection"].as<std::string>() == "tri" ||
      par["collection"].as<std::string>() == "merge")
    colalg = par["collection"].as<std::string>();
  else
    colfns.push_back(par["collection"].as<std::string>());

  if (par["collection"].as<std::string>() == "merge") {
    colfns = par["collections"].as<std::vector<std::string>>();
  }

  if(!par.count("instance") && colfns.empty()) {
    std::cout << "Instance or collection parameter is required" << std::endl;
    std::cout << options.help() << std::endl;
    return 2;
  }

  bool autoAlg = false;
  std::string algorithm;
  algorithm = par["algorithm"].as<std::string>();
  if(algorithm == "auto") {
    autoAlg = true;
    algorithm = "mip";
  }

  if(instfn.empty()) {
    std::string basename = baseName(colfns[0]);
    instfn = instpath + basename + ".instance.json";
  }

  std::string outpath;
  if (par.count("outpath")) {
    outpath = par["outpath"].as<std::string>();
    if(outpath.back() != '/')
      outpath.push_back('/');
  }

  int fix = par["fix"].as<int>();

  //-----------------------------------------------------------------
  // Prepare instance, collection, and solution
  //-----------------------------------------------------------------

  // Create instance and solution
  Instance instance(instfn, par["bag"].as<bool>());

  instance.triangleLimit = par["triangle"].as<double>();
  Solution sol;
  sol.instance_name = instance.name;
  if(autoAlg)
    sol.algorithm = "auto";
  else
    sol.algorithm = algorithm;

  // Create collection
  Solution col = colfns.empty() ? Solution() : Solution(colpath + colfns[0]);

  if(colfns.size() > 1) {
    std::set<std::vector<Point>> cpolys;
    Message m("readmerge",3, "Reading collections");
    for(int i = 0; i < colfns.size(); i++) {
      Solution stmp(colpath + colfns[i], 1);
      for(Convex &c : stmp.polys) {
        std::vector v(c.begin(),c.end());
        std::sort(v.begin(), v.end());
        cpolys.insert(v);
      }
    }
    col.polys.clear();
    for(const std::vector<Point> &pol : cpolys) {
      col.polys.push_back(Convex(pol));
    }
    m.close(col.polys.size());
  }

  col.instance_name = instance.name;
  if(colfns.empty()) {
    col.algorithm = colalg;
    if(par["ngrow"].as<int>() > 1)
      col.algorithm += std::to_string(par["ngrow"].as<int>());
    if(ipts==1)
      col.algorithm += 'p';
    else if(ipts==2)
      col.algorithm += "pp";
  }
  else if(colfns.size() ==1) {
    col.algorithm = colName(colfns[0]);
  }
  else {
    col.algorithm = "merge";
  }

  //-----------------------------------------------------------------
  // Preprocess instance
  //-----------------------------------------------------------------
  // Add intermediate points
  if(ipts && !colalg.empty()) {
    Message Mep("find_extra_points", 4, "Finding extra points");
    std::vector<Point> extra;
    if(ipts == 1)
      extra = edgeEdge(instance);
    else if(ipts == 2)
      extra = edgeEdgePlus(instance);
    Mep.close(extra.size());
    Message Maep("add_extra_points", 4, "Adding extra points");
    int c = instance.registerPoints(extra);
    // instance.writeData();
    Maep.close(c);
  }

  //-----------------------------------------------------------------
  // Build collection
  //-----------------------------------------------------------------
  // Populate collection
  if(colfns.empty()) {
    if(colalg == "max") {
      Message Mms("max_sets",3,"Building maximal sets");
      col.polys = maximalConvex(instance);
      Mms.close(col.polys.size());
    }
    else if(colalg == "rgrow") {
      Message Mgt("grown",4,"Building grown triangulation");
      col.polys = rgrownTriangles(instance, par["ngrow"].as<int>());
      Mgt.close(col.polys.size());
    }
    else if(colalg == "tri") {
      Message Mt("triangulation",4,"Building triangulation");
      col.polys = triangulate(instance.poly);
      Mt.close(col.polys.size());
    }
    else if(colalg == "quit")
      return 0;
  }

  // Save collection
  if(savecol)
    col.writeSol(false, colpath);

  if(bloat) {
    if(bloat==1)
      col.algorithm += 'b';
    else
      col.algorithm += 'B';

    col.polys = bloatPolys(instance, col.polys, bloat);
    // Save collection
    if(savecol)
      col.writeSol(false, colpath);
  }

  sol.algorithm = col.algorithm + "_" + sol.algorithm;

  // Build witnesses
  std::vector<Witness> witnesses;
  Message Mwit("witnesses",4,"Finding witnesses");
  if(wit == 1)
    witnesses = findWitnesses(col.polys, instance);
  else if(wit == 2)
    witnesses = findVertexWitnesses(col.polys, instance);
  else
    witnesses = findQuickVertexWitnesses(col.polys, instance);
  Mwit.close(witnesses.size());

  //-----------------------------------------------------------------
  // Run algorithm to find a solution
  //-----------------------------------------------------------------

  // Find solution
  int bestSol = 0;
  std::vector<Convex> uncovered;
  Message Msc("set_cover", 4, "Initializing set cover");
  SetCover cover(col.polys, witnesses);
  Msc.close(std::to_string(col.polys.size()) + "," +
            std::to_string(witnesses.size()));

  int k = 0, lastFix = -1, equalFix = 0;
  do {
    int iter = par["iter"].as<int>();
    int del = par["del"].as<int>();
    Message Malg(algorithm, 2, "Running " + algorithm);
    double mipGap = 0.005;
    if(algorithm == "mip")
      mipGap = cover.mip(par["to"].as<int>());
    else if(algorithm == "greedy")
      cover.greedy();
    else if(algorithm == "annealing")
      cover.annealing(iter, del, par["repeat"].as<int>(), par["to"].as<int>());
    else {
      std::cerr << "Unknown algorithm: " << algorithm << std::endl;
      return 3;
    }
    sol.polys = cover.getSolution();
    Malg.close(sol.polys.size());

    // Test the solution found and add more witnesses or fix
    Message Mv("uncovered", 2, "Finding uncovered regions");
    uncovered = uncoveredFaces(sol.polys, instance);
    for(Convex &conv : uncovered) {
      Point inpt = conv.pointInside();
      cover.addUncoveredWitness(Witness(inpt,inpt));
    }
    if(lastFix == uncovered.size()) {
      equalFix++;
      if(equalFix > 20) {
        std::cout << "It looks like we have a bogus loop!!!";
        exit(1);
      }
    }
    else {
      lastFix = uncovered.size();
      equalFix = 0;
    }

    Mv.close(uncovered.size());

    if(uncovered.size() != 0) {
      Message Mfix("fix",2,"Fixing uncovered regions");
      if(uncovered.size() <= 100)
        fuseConvexes(uncovered, instance);
      Mfix.close(uncovered.size());
      if(bestSol == 0 || uncovered.size() + sol.polys.size() < bestSol) {
        for(Convex &conv : uncovered)
          sol.polys.push_back(conv);
        bestSol = sol.polys.size();
        if(savetmp)
          sol.writeSol(false, outpath, "tmp");
        for(int i = 0; i < uncovered.size(); i++)
          sol.polys.pop_back();
      }

      if(uncovered.size() <= fix || uncovered.size() < mipGap * sol.polys.size()) {
        for(Convex &conv : uncovered)
          sol.polys.push_back(conv);
        uncovered.clear(); // To exit the do_while
      }
    }

    if(autoAlg && mipGap > .01) {
      uncovered.push_back(Convex()); // To repeat the do_while
      algorithm = "annealing";
    }

    k++;
  } while(!uncovered.empty());
  if(VERBOSE>3 && k > 1)
    std::cout << "Total of " << k << " iterations in " << int(elapsed()) << " sec" << std::endl;

  //-----------------------------------------------------------------
  // Save solution
  //-----------------------------------------------------------------

  if(bestSol != 0 && sol.polys.size() > bestSol)
    sol.writeSol(false, outpath, "worse");
  else{
    sol.writeSol(false, outpath);
  }

  return 0;
}

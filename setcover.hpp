#pragma once
#include "convex.hpp"

using SetIndex = int;
using WitIndex = int;

class SetCover {
  std::vector<Witness> witnesses;
  std::vector<Convex> sets;
  std::vector<std::vector<WitIndex>> witnessesOf;
  std::vector<tsl::sparse_set<WitIndex>> wouldCover;
  std::set<std::pair<int,SetIndex>> wouldCoverPairs;
  std::vector<std::vector<SetIndex>> setsContaining;
  std::vector<tsl::sparse_set<SetIndex>> currentlyCoveredBy;
  tsl::sparse_set<WitIndex> uncovered;
  tsl::sparse_set<SetIndex> solution;
  std::vector<tsl::sparse_set<WitIndex>> coveringAlone;

public:
  SetCover(const std::vector<Convex> &_sets, const std::vector<Witness>  &_witnesses)
      : sets(_sets), witnesses(_witnesses) {
    std::sort(witnesses.begin(), witnesses.end(), witPointLessThan);
    for(WitIndex i = 0; i < witnesses.size(); i++)
      uncovered.insert(i);

    setsContaining = std::vector<std::vector<SetIndex>>(witnesses.size());
    currentlyCoveredBy = std::vector<tsl::sparse_set<SetIndex>>(witnesses.size());
    coveringAlone = std::vector<tsl::sparse_set<WitIndex>>(sets.size());

    for(SetIndex si = 0; si < sets.size(); si++) {
      std::vector<WitIndex> wits;
      Witness start = Witness(sets[si].boxMin, sets[si].boxMin);

      // Should only check points with the x-coordinate in the bbox range
      for(auto it = std::lower_bound(witnesses.begin(), witnesses.end(), start, witPointLessThan);
          it != witnesses.end() && it->source() <= sets[si].boxMax;
          ++it) {
        Witness w = *it;
        int wi = it - witnesses.begin();
        if(sets[si].contains(w)) {
          wits.push_back(wi);
          setsContaining[wi].push_back(si);
        }
      }
      witnessesOf.push_back(wits);
      wouldCover.push_back(tsl::sparse_set<WitIndex>(wits.begin(), wits.end()));
    }

    for(SetIndex i = 0; i < wouldCover.size(); i++)
      wouldCoverPairs.insert(std::make_pair(wouldCover[i].size(), i));

    for(WitIndex wi = 0; wi < witnesses.size(); wi++) {
      if(setsContaining[wi].empty()) {
        std::cout<< "Collection does not cover witness" << witnesses[wi] << std::endl;
      }
    }
  }

  void clear() {
    solution.clear();
    uncovered.clear();

    for(WitIndex i = 0; i < witnesses.size(); i++)
      uncovered.insert(i);
    currentlyCoveredBy = std::vector<tsl::sparse_set<SetIndex>>(witnesses.size());
    coveringAlone = std::vector<tsl::sparse_set<WitIndex>>(sets.size());
    wouldCover.clear();
    for(const std::vector<WitIndex> &wits : witnessesOf)
      wouldCover.push_back(tsl::sparse_set<WitIndex>(wits.begin(), wits.end()));
    wouldCoverPairs.clear();
    for(SetIndex i = 0; i < wouldCover.size(); i++)
      wouldCoverPairs.insert(std::make_pair(wouldCover[i].size(), i));
  }

  void addUncoveredWitness(Witness p) {
    witnesses.push_back(p);
    WitIndex pi = witnesses.size()-1;
    uncovered.insert(pi);
    setsContaining.push_back(std::vector<SetIndex>());
    currentlyCoveredBy.push_back(tsl::sparse_set<SetIndex>());

    // std::sort(witnesses.begin(), witnesses.end());
    bool covered = false;
    for(SetIndex si = 0; si < sets.size(); si++) {
      if(sets[si].contains(witnesses[pi])) {
        covered = true;
        assert(!solution.contains(si));
        setsContaining[pi].push_back(si);
        witnessesOf[si].push_back(pi);
        wouldCoverPairs.erase(std::make_pair(wouldCover[si].size(), si));
        wouldCover[si].insert(pi);
        wouldCoverPairs.insert(std::make_pair(wouldCover[si].size(), si));
      }
    }
    if(!covered) {
      std::cout << "Error, witness not covered: " << p << std::endl;
      exit(1);
    }
  }

  void solutionAdd(SetIndex s) {
    tsl::sparse_set<SetIndex> toRemove;
    for(WitIndex p : witnessesOf[s]) {
      tsl::sparse_set<SetIndex> &cover_p = currentlyCoveredBy[p];

      if(cover_p.empty()) {
        uncovered.erase(p);
        for(SetIndex i : setsContaining[p]) {
          tsl::sparse_set<WitIndex> &wci = wouldCover[i];
          wouldCoverPairs.erase(std::make_pair(wci.size(), i));
          wci.erase(p);
          wouldCoverPairs.insert(std::make_pair(wci.size(), i));
        }
        coveringAlone[s].insert(p);
      }
      else if(cover_p.size() == 1) {
        for(SetIndex i : cover_p) { // Only one
          coveringAlone[i].erase(p);
          if(coveringAlone[i].empty()) {
            toRemove.insert(i);
          }
        }
      }
      cover_p.insert(s);
    }

    for(SetIndex i : toRemove)
      if(coveringAlone[i].empty()) // Check again because a previous remove may change it
        solutionRemove(i);

    solution.insert(s);
  }

  void solutionRemove(SetIndex s) {
    if(!solution.contains(s))
      return;
    coveringAlone[s].clear();

    for(WitIndex p : witnessesOf[s]) {
      tsl::sparse_set<SetIndex> &cover_p = currentlyCoveredBy[p];
      cover_p.erase(s);

      if(cover_p.empty()) {
        uncovered.insert(p);
        for(SetIndex i : setsContaining[p]) {
          tsl::sparse_set<WitIndex> &wci = wouldCover[i];
          wouldCoverPairs.erase(std::make_pair(wci.size(), i));
          wci.insert(p);
          wouldCoverPairs.insert(std::make_pair(wci.size(), i));
        }
      }
      if(cover_p.size() ==1) {
        for(SetIndex i : cover_p) { // Only one
          coveringAlone[i].insert(p);
        }
      }
    }

    solution.erase(s);
  }

  void solutionUpdate(const tsl::sparse_set<SetIndex> &newsol) {
    for(SetIndex s : tsl::sparse_set(solution))
      if(!newsol.contains(s))
        solutionRemove(s);

    for(SetIndex s : newsol)
      if(!solution.contains(s))
        solutionAdd(s);
  }

  SetIndex bestSetIndex() const {
    std::vector<SetIndex> v;
    int bestFirst = wouldCoverPairs.rbegin()->first; // Max element
    for(auto iter = wouldCoverPairs.rbegin();
        iter != wouldCoverPairs.rend() && iter->first == bestFirst;
        iter++) {
      v.push_back(iter->second);
    }
    return v[rand()%v.size()];
  }

  void greedy() {
    while(!uncovered.empty()) {
      solutionAdd(bestSetIndex());
    }
  }

  void randomRemove(int removals) {
    std::vector<SetIndex> solv(solution.begin(), solution.end());
    for(int r = 0; r < removals; r++) {
      int rnd = rand() % solv.size();
      solutionRemove(solv[rnd]);
    }
  }

  double mip(double tlimit = 1800, bool reset = true, double mlimit = 7000) {
    if(reset)
      clear();

    IloEnv env;

    // Creating variables
    std::map<int, IloNumVar> variables;
    for(auto iter = wouldCoverPairs.rbegin();
        iter != wouldCoverPairs.rend() && iter->first != 0;
        iter++) {
      SetIndex s = iter->second;
      variables[s] = IloNumVar(env, 0, 1, ILOINT);
    }

    IloModel model(env);

    // Creating constraints
    for(WitIndex p : uncovered) {
      IloExpr expr(env);
      // std::cout << p << std::endl;
      for(int s : setsContaining.at(p)) {
        assert(variables.contains(s));
        expr += variables[s];
      }
      model.add(expr >= 1);
    }

    // Setting objective
    // Minimize the number of vertices
    IloExpr expr(env);
    for(auto &[_,var] : variables) {
      expr += var;
    }
    model.add(IloMinimize(env,expr));

    // Solving
    IloCplex cplex(model);
    if(VERBOSE < 8)
      cplex.setOut(env.getNullStream()); // Disable console output
    cplex.setParam(IloCplex::Param::TimeLimit, tlimit); // time limit in seconds
    cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, mlimit);
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, .001);

    cplex.setParam(IloCplex::Param::Threads, 1); // Single thread
    bool solved = cplex.solve();

    std::ostringstream ss;

    double ret = 1e10;

    // Saving solution
    if(solved) {
      ss << cplex.getCplexStatus() << '(' << cplex.getMIPRelativeGap() << ')';
      meta["cplex_status"] += ss.str() + ";";
      std::cout << " " << ss.str() << " ";
      ret = cplex.getMIPRelativeGap();

      // Pick sets of value 1 and add to dom
      for(auto &[s,var] : variables) {
        if(cplex.getValue(var) > .5)
          solutionAdd(s);
      }
    }

    return ret;
  }

  std::vector<Convex> getSolution() const {
    std::vector<Convex> ret;
    for(SetIndex s : solution)
      ret.push_back(sets[s]);

    return ret;
  }

  void annealing(int iter, int removals, int repeat, double timeout) {
    int goal = solution.size();
    greedy();
    tsl::sparse_set<SetIndex> best = solution, previous = solution;
    const double t0 = 100;
    double sTime = elapsed();

    for(int r = 1; r <= repeat && elapsed() - sTime <= timeout && solution.size() > goal; r++) {
      for(int i = 1; i <= iter; i++) {
        double temperature = t0 / i;
        randomRemove(removals);
        greedy();
        if(solution.size() <= best.size()) {
          best = solution;
          previous = solution;
        }
        else {
          double diff = 100.0 * (double(previous.size()) - double(solution.size())) /
                        double(solution.size());
          double changeProb = exp(diff / temperature);
          if(double(rand()) / RAND_MAX > changeProb) {
            solutionUpdate(previous); // Change back to previous
          }
          else
            previous = solution;
        }
      }
      solutionUpdate(best);
      previous = best;
      if(VERBOSE > 2)
        std::cout << " " << solution.size() << std::flush;
    }
  }

};

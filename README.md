# shadoks-CGSHOP2023
Shadoks CG:SHOP 2023 challenge code (convex covering)

https://cgshop.ibr.cs.tu-bs.de/competition/cg-shop-2023/#problem-description

To compile, download the current CGAL from https://github.com/CGAL/cgal and CPLEX
and run:

```cmake -DCMAKE_BUILD_TYPE=Release -DCGAL_DIR=~/pathtocgal/ -DCPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio221/ -DCMAKE_MODULE_PATH=.  CMakeLists.txt```
```make```

To run, simply type:

```./solver```

Which will show a help message. If you want to go with the default settings, all you need is an instance file:

```./solver -i someinstance.instance.json```

Details about the algorithms are available at

https://arxiv.org/abs/2303.07696

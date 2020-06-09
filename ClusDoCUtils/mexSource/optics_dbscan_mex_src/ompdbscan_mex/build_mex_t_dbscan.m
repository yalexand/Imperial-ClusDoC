
%mex -v -largeArrayDims CXXFLAGS="$CXXFLAGS -fPIC -fopenmp" -I/usr/include -cxx -c kdtree2.cpp utils.cpp kdtree2rnearest.cpp

%mex -v -largeArrayDims CXXFLAGS="$CXXFLAGS -fPIC -fopenmp" CXXLIBS="$CXXLIBS -fopenmp" -cxx kdtree2.o utils.o kdtree2rnearest.o -output kdtree2rnearest



mex -v -largeArrayDims CXXFLAGS="$CXXFLAGS -fPIC -fopenmp" -I/usr/include -cxx -c kdtree2.cpp dbscan.cpp utils.cpp clusters.cpp omp_dbscan.cpp pdsdbscan.cpp

mex -v -largeArrayDims CXXFLAGS="$CXXFLAGS -fPIC -fopenmp" CXXLIBS="$CXXLIBS -fopenmp" -cxx kdtree2.o dbscan.o utils.o clusters.o omp_dbscan.o pdsdbscan.o -output pdsdbscan

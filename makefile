all: val_test01_solved val_test02_solved MMult1 omp_solved2 omp_solved3 omp_solved4 omp_solved5 omp_solved6 jacobi2D-omp gs2D-omp

  
val_test01_solved: val_test01_solved.cpp
	g++ -o val_test01_solved -std=c++11 val_test01_solved.cpp  
val_test02_solved: val_test02_solved.cpp 
	g++ -o val_test02_solved -std=c++11 val_test02_solved.cpp  
MMult1: MMult1.cpp
	g++ -o MMult1 -fopenmp -O3 -march=native -std=c++11 MMult1.cpp
omp_solved2: omp_solved2.cpp
	g++ -o omp_solved2 -fopenmp -std=c++11 omp_solved2.cpp  
omp_solved3: omp_solved3.cpp 
	g++ -o omp_solved3 -fopenmp -std=c++11 omp_solved3.cpp  
omp_solved4: omp_solved4.cpp 
	g++ -o omp_solved4 -fopenmp -std=c++11 omp_solved4.cpp  
omp_solved5: omp_solved5.cpp 
	g++ -o omp_solved5 -fopenmp -std=c++11 omp_solved5.cpp  
omp_solved6: omp_solved6.cpp 
	g++ -o omp_solved6 -fopenmp -std=c++11 omp_solved6.cpp  
jacobi2D-omp: jacobi2D-omp.cpp
	g++ -o jacobi2D-omp -fopenmp -O3 -std=c++11 jacobi2D-omp.cpp
gs2D-omp: gs2D-omp.cpp
	g++ -o gs2D-omp -fopenmp -O3 -std=c++11 gs2D-omp.cpp

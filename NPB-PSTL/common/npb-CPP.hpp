/*
MIT License

Copyright (c) 2021 Parallel Applications Modelling Group - GMAP 
	GMAP website: https://gmap.pucrs.br
	
	Pontifical Catholic University of Rio Grande do Sul (PUCRS)
	Av. Ipiranga, 6681, Porto Alegre - Brazil, 90619-900

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

------------------------------------------------------------------------------

The original NPB 3.4.1 version was written in Fortran and belongs to: 
	http://www.nas.nasa.gov/Software/NPB/

------------------------------------------------------------------------------

The serial C++ version is a translation of the original NPB 3.4.1
Serial C++ version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-SER

Authors of the C++ code: 
	Dalvan Griebler <dalvangriebler@gmail.com>
	Gabriell Araujo <hexenoften@gmail.com>
 	Júnior Löff <loffjh@gmail.com>
*/ 

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <tuple>
#include <complex>
#include <algorithm>
#include <execution>
#include <span>
#include <ranges>
#include <numeric>

typedef int boolean;
typedef struct { double real; double imag; } dcomplex;

class CountIterator {
  public:
    CountIterator(const int size){
      iter_vec.resize(size);
      std::iota(iter_vec.begin(), iter_vec.end(), 0);
    };

    std::vector<int>::iterator front(){
      return iter_vec.begin();
    }

    std::vector<int>::iterator tail(){
      return iter_vec.end();
    }

    std::vector<int>::reverse_iterator rfront(){
      return iter_vec.rbegin();
    }

    std::vector<int>::reverse_iterator rtail(){
      return iter_vec.rend();
    }

  private:
    std::vector<int> iter_vec;
};

#define TRUE	1
#define FALSE	0

#define max_npb(a,b) (((a) > (b)) ? (a) : (b))
#define min_npb(a,b) (((a) < (b)) ? (a) : (b))
// using std::max;
// using std::min;
#define	pow2(a) ((a)*(a))




template <typename Matrix>
class Matrix2DIterator {
public:
    using ValueType = typename Matrix::ValueType;
    using PointerType = ValueType*;
    using ReferenceType = ValueType&;
public:
    Matrix* matrixPtr;
    PointerType currPtr;
    size_t row_size;
    size_t col_size;
    int row_;
    int col_;

    Matrix2DIterator(Matrix& matrix, size_t Max_row=0, size_t Max_col=0, size_t row = 0, size_t col = 0)
        : currPtr(nullptr), row_(row), col_(col), row_size(Max_row), col_size(Max_col) {
            matrixPtr = &matrix;
            currPtr = &(matrix[row_][col_]);
        }

    Matrix2DIterator& operator++() {
        col_++;
        currPtr++;
        if (col_ >= col_size and row_ < row_size-1) {
            col_ = 0;
            row_++;
            currPtr = &((*matrixPtr)[row_][col_]);
        }
        return *this;
    }

    Matrix2DIterator operator++(int) {
        Matrix2DIterator iterator = *this;
        ++(*this);
        return iterator;
    }

    Matrix2DIterator operator+(int i) {
        if (col_ > col_size and row_ == row_size-1)
            return Matrix2DIterator((*matrixPtr), row_size, col_size, row_, col_+i);
            
        int col_calc = (i+col_) % col_size;
        int row_calc = row_+((i+col_) / col_size);
        if (row_calc > row_size-1){
            col_calc += (row_calc-row_size-1)*(col_size-1);
            row_calc = row_size-1;
        }
        return Matrix2DIterator((*matrixPtr), row_size, col_size, row_calc, col_calc);
    }

    Matrix2DIterator& operator--() {
        if (col_<= 0 and 0 < row_) {
            col_ = col_size-1;
            row_--;
            currPtr = &((*matrixPtr)[row_][col_]);
        }
        else{
            col_--;
            currPtr--;}
        return *this;
    }

    Matrix2DIterator operator-(int i) {
        if (col_ <= 0 and row_ == 0)
            return Matrix2DIterator((*matrixPtr), row_size, col_size, row_, col_-i);

        int col_calc = (col_-i) % col_size;
        int row_calc = row_+((col_-i) / col_size);
        if (row_calc < 0){
            col_calc += (row_calc)*col_size;
            row_calc = 0;}
        return Matrix2DIterator((*matrixPtr), row_size, col_size, row_calc, col_calc);
    }


    Matrix2DIterator operator--(int) {
        Matrix2DIterator iterator = *this;
        --(*this);
        return iterator;
    }

    bool operator==(const Matrix2DIterator& other) const {
        return row_ == other.row_ and col_ == (other.col_);
    }

    bool operator!=(const Matrix2DIterator& other) const {
        return !(*this == other);
    }

    ValueType& operator*() {
        return *currPtr;
    }

    ReferenceType operator->() {
        return (*matrixPtr)[col_][row_];
    }

};




template <typename T>
class Matrix2D {
    public:
        using ValueType = T;
        using Iterator = Matrix2DIterator<Matrix2D<ValueType>>;
    public:
        int rows, columns;

        std::vector<std::vector<T>> matrix;

        Matrix2D(int n_row, int n_col, T init_value): rows(n_row), columns(n_col) {
          matrix.resize(rows, std::vector<T>(columns, init_value));
        };

        Matrix2D(): rows(0), columns(0) {};

        Matrix2D(int n_row, int n_col, int slice, std::vector<double> &base_vec, std::vector<int> dims)
          : rows(n_row), columns(n_col){
            matrix.resize(rows, std::vector<T>(columns));
            int count=slice;

            CountIterator iter(*std::max_element(dims.begin(), dims.end()));
            auto policy = std::execution::par;
            std::for_each(policy, iter.front(), iter.front()+dims[0], [&](int i3){
                int pre_comp = i3*dims[1];
                std::for_each(iter.front(), iter.front()+dims[0], [&](int i2){
                    matrix[i3][i2] = &base_vec[(pre_comp+i2)*dims[2]+slice];
                });
            });

        };

        Matrix2D(Matrix2D& base_class):rows(base_class.rows), columns(base_class.columns) {
            matrix = base_class.matrix;
        };

        std::vector<T>& operator[] (int index) {
            return matrix[index];
        };
        
        Matrix2D& operator=(const Matrix2D& base_class) {
            rows = base_class.rows;
            columns = base_class.columns;
            matrix = base_class.matrix;
            return *this;
        };

        Iterator begin() { return Iterator(*this, rows, columns); }
        Iterator end() { return Iterator(*this, rows, columns, rows-1, columns); }
};

template <typename T>
class Matrix3D {
    public:
        int rows, columns, depth;

        std::vector<std::vector<std::vector<T>>> cube;

        Matrix3D(int n_row, int n_col, int n_depth, T init_value): rows(n_row), columns(n_col), depth(n_depth) {
			cube.resize(rows, std::vector<std::vector<T>>(columns,std::vector<T>(depth)));
        };

        Matrix3D(Matrix3D& base_class):rows(base_class.rows), columns(base_class.columns), depth(base_class.depth) {
            cube = base_class.cube;
        };

        std::vector<std::vector<T>>& operator[] (int index) {
            return cube[index];
        };
};

template <typename T>
class Matrix4D {
    public:
        int x_val, y_val, z_val, k_val;

        std::vector<std::vector<std::vector<std::vector<T>>>> tensor;

        Matrix4D(int n_x_val, int n_y_val, int n_z_val, int n_k_val, T init_value):
          x_val(n_x_val), y_val(n_y_val), z_val(n_z_val), k_val(n_k_val) {
			      tensor.resize(x_val, std::vector<std::vector<std::vector<T>>>(y_val,std::vector<std::vector<T>>(z_val, std::vector<T>(k_val, init_value))));
        };

        Matrix4D(Matrix4D& base_class):
          x_val(base_class.x_val),
          y_val(base_class.y_val),
          z_val(base_class.z_val),
          k_val(base_class.k_val) {
            tensor = base_class.tensor;
        };

        std::vector<std::vector<std::vector<T>>>& operator[] (int index) {
            return tensor[index];
        };
};

/* latest version of the complex number operations */
#define dcomplex_create(r,i) (dcomplex){r, i}
#define dcomplex_add(a,b) (dcomplex){(a).real+(b).real, (a).imag+(b).imag}
#define dcomplex_sub(a,b) (dcomplex){(a).real-(b).real, (a).imag-(b).imag}
#define dcomplex_mul(a,b) (dcomplex){((a).real*(b).real)-((a).imag*(b).imag),\
	((a).real*(b).imag)+((a).imag*(b).real)}
#define dcomplex_mul2(a,b) (dcomplex){(a).real*(b), (a).imag*(b)}
static inline dcomplex dcomplex_div(dcomplex z1, dcomplex z2){
	double a = z1.real;
	double b = z1.imag;
	double c = z2.real;
	double d = z2.imag;
	double divisor = c*c + d*d;
	double real = (a*c + b*d) / divisor;
	double imag = (b*c - a*d) / divisor;
	dcomplex result = (dcomplex){real, imag};
	return result;
}
#define dcomplex_div2(a,b) (dcomplex){(a).real/(b), (a).imag/(b)}
#define dcomplex_abs(x)    sqrt(((x).real*(x).real) + ((x).imag*(x).imag))
#define dconjg(x)          (dcomplex){(x).real, -1.0*(x).imag}

extern double randlc(double &, double);
extern void vranlc(int, double &, double, double *);
extern void timer_clear(int);
extern void timer_start(int);
extern void timer_stop(int);
extern double timer_read(int);

extern void c_print_results(std::string name,
		char class_npb,
		int n1,
		int n2,
		int n3,
		int niter,
		double t,
		double mops,
		std::string optype,
		int passed_verification,
		std::string npbversion,
		std::string compiletime,
		std::string compilerversion,
		std::string cc,
		std::string clink,
		std::string c_lib,
		std::string c_inc,
		std::string cflags,
		std::string clinkflags,
		std::string rand);

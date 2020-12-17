#include <iostream>
#include <fstream>
#include "myarrays.h"
#include <cmath>

//making a uniform 1d array
MyArray1D::MyArray1D(int size) {
	mData = new double[size];
	mSize = size;
}

MyArray1D::~MyArray1D() {
	delete[] mData;
}

int My1DElement::count = 0; //initialising the static counter

//making a 1d element
My1DElement::My1DElement() {
	//count++; //counting number of objects we have created
	//mLength = length;
	double time_steps = 20;
	//mRefractive_index = refractive_index;
	mNumber_of_nodes = 2; //assume we have 2 nodes as a default
	mE = new double* [2]; //electric field
	for (int i = 0; i < 2; i++)
		mE[i] = new double[time_steps];
	mH = new double* [2]; //magnetic field
	for (int i = 0; i < 2; i++)
		mH[i] = new double[time_steps];
	mZ0H = new double* [2]; //magnetic field
	for (int i = 0; i < 2; i++)
		mZ0H[i] = new double[time_steps];
	//mTime_steps = time_steps;
	//Mid_point = length * mElementNumber - 2; //to define our midpoint at 1 for the first element
	//mElementNumber = ele_num;
	//mElementNumber = new double;
	//mZlocal = sqrt(mMu / mEpsilon);
	//mYlocal = sqrt(mEpsilon / mMu);
}





My1DElement::My1DElement(double length, double time_steps) {
	count++; //counting number of objects we have created
	mLength = length;
	//mRefractive_index = refractive_index;
	mNumber_of_nodes = 2; //assume we have 2 nodes as a default
	mE = new double* [2]; //electric field
	for (int i = 0; i < 2; i++)
		mE[i] = new double[time_steps];
	mH = new double* [2]; //magnetic field
	for (int i = 0; i < 2; i++)
		mH[i] = new double[time_steps];
	mTime_steps = time_steps;
	mMid_point = length * mElementNumber - 2; //to define our midpoint at 1 for the first element
	//mElementNumber = ele_num;
	//mElementNumber = new double;
	//mZlocal = sqrt(mMu / mEpsilon);
	//mYlocal = sqrt(mEpsilon / mMu);
}

My1DElement::My1DElement(double length, int order_of_polynomial, double time_steps) {
	count++; //counting number of objects we have created
	mLength = length;
	//mRefractive_index = refractive_index;
	mNumber_of_nodes = order_of_polynomial + 1; //use when changing the number of modes
	mE = new double* [2]; //electric field
	for (int i = 0; i < 2; i++)
		mE[i] = new double[time_steps];
	mH = new double* [2]; //magnetic field
	for (int i = 0; i < 2; i++)
		mH[i] = new double[time_steps];
	mTime_steps = time_steps;
	mMid_point = length * mElementNumber - 2; //to define our midpoint at 1 for the first element
	//mElementNumber = new double;
	//mZlocal = sqrt(mMu / mEpsilon);
	//mYlocal = sqrt(mEpsilon / mMu);
}

My1DElement::~My1DElement() {
	delete[] mData;
}




//making a uniform 2d array
MyArray2D::MyArray2D(int row_dim, int col_dim) {
	mData = new double* [row_dim];
	for (int i = 0; i < row_dim; i++)
		mData[i] = new double[col_dim];
	mrow_dim = row_dim;
	mcol_dim = col_dim;
}

MyArray2D::~MyArray2D() {
	delete[] mData;
}

//function defintions
double g25(int t) {
	//return exp(-0.5 * (25 - t) * (25 - t) / (100));
	//trying a wider Gaussian
	return exp(-0.5 * (25 - t) * (25 - t) / (6250)) * sin(2 * 3.14 * 400 * t);
}
double g50(int t) {
	//return exp(-0.5 * (50 - t) * (50 - t) / (100));
	//trying a wider Gaussian
	return exp(-0.5 * (50 - t) * (50 - t) / (6250)) * sin(2 * 3.14 * 400 * t);
}


//function to calculate the transpose of matrices
//2x2 matrix
void transpose_matrix_2x2(double** arr, double** tran_arr) {
	tran_arr[0][0] = arr[0][0];
	tran_arr[1][1] = arr[1][1];
	tran_arr[0][1] = arr[1][0];
	tran_arr[1][0] = arr[0][1];
}

//3x3 matrix
void transpose_matrix_3x3(double** arr, double** tran_arr) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			tran_arr[i][j] = arr[j][i];
		}
	}
}

//void exam_print(double impedance) {
//	std::cout << impedance;
//};







//function to invert matrix (outputs a 2d matrix), this is for a 3x3 matrix:
//double** inverting_matrix_3x3(double** arr) {
void inverting_matrix_3x3(double** arr, double** inv_arr) {
	//step 1
	//Find determinant of the matrix
	double det1 = arr[1][1] * arr[2][2] - arr[2][1] * arr[1][2];
	double det2 = arr[1][0] * arr[2][2] - arr[2][0] * arr[1][2];
	double det3 = arr[1][0] * arr[2][1] - arr[2][0] * arr[1][1];
	double det_arr = arr[0][0] * det1 - arr[0][1] * det2 + arr[0][2] * det3;
	MyArray2D adj_arr_g(3, 3); //forming a 3x3 adj matrix
	double** adj_arr = adj_arr_g.getdata();
	//MyArray2D inverse_arr_g(3, 3); //forming a 3x3 inverse matrix
	//double** inverse_arr = inverse_arr_g.getdata();
	adj_arr[0][0] = arr[1][1] * arr[2][2] - arr[2][1] * arr[1][2];
	adj_arr[0][1] = -(arr[0][1] * arr[2][2] - arr[2][1] * arr[0][2]);
	adj_arr[0][2] = arr[0][1] * arr[1][2] - arr[1][1] * arr[0][2];
	adj_arr[1][0] = -(arr[1][0] * arr[2][2] - arr[2][0] * arr[1][2]);
	adj_arr[1][1] = arr[0][0] * arr[2][2] - arr[2][0] * arr[0][2];
	adj_arr[1][2] = -(arr[0][0] * arr[1][2] - arr[1][0] * arr[0][2]);
	adj_arr[2][0] = arr[1][0] * arr[2][1] - arr[2][0] * arr[1][1];
	adj_arr[2][1] = -(arr[0][0] * arr[2][1] - arr[2][0] * arr[0][1]);
	adj_arr[2][2] = arr[0][0] * arr[1][1] - arr[1][0] * arr[0][1];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			inv_arr[i][j] = (1 / det_arr) * adj_arr[i][j];
		}
	}
	//return inv_arr; //now a void function so don't need to return anything
}

//function to invert 2x2 matrix
//double** inverting_matrix_2x2(double** arr) {
void inverting_matrix_2x2(double** arr, double** inv_arr) {
	//finding determinant of 2x2 matrix:
	double det = arr[0][0] * arr[1][1] - arr[1][0] * arr[0][1];
	MyArray2D adj_arr_g(2, 2);
	double** adj_arr = adj_arr_g.getdata();
	//MyArray2D inv_arr_g(2, 2);
	//double** inv_arr = adj_arr_g.getdata();
	adj_arr[0][0] = arr[1][1];
	adj_arr[1][1] = arr[0][0];
	adj_arr[0][1] = -arr[0][1];
	adj_arr[1][0] = -arr[1][0];
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			inv_arr[i][j] = (1 / det) * adj_arr[i][j];
		}
	}
	//return inv_arr; //this failed because of dynamic memory allocation - adj_arr would have gone out of scope.
}

//function to update the Runge Kutta coefficients, for first order polynomial approximation
/*void updating_Runge_Kutta_coefficients(double** mass_m, double** stiff_m, double** face_m_LHS, double** face_m_RHS, double** m_inv, double alpha, double element_length) { //will calculate all our updated fields!
	My1DElement elements[100]; //can change for total number of elements from 100
	//setting the time steps
	double time_steps = 20;
	//defining the time step size
	double c = 299792458; //speed of light
	double dt = element_length * 0.01 / c; //divide by c to set in correct units

	//defining free space properties
	double epsilon_freespace = 8.85e-12;
	double mu_freespace = 12.57e-7;
	double Z_freespace = sqrt(mu_freespace / epsilon_freespace);
	double Y_freespace = sqrt(epsilon_freespace / mu_freespace);
	//making freespace E-field and H-field 2D arrays
	MyArray2D E_free_arr(2, time_steps);
	double** E_free = E_free_arr.getdata();
	MyArray2D H_free_arr(2, time_steps);
	double** H_free = H_free_arr.getdata();
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < time_steps; j++) {
			E_free[i][j] = 0; //setting freespace E-field to 0
			H_free[i][j] = 0; //setting freespace H-field to 0
		}
	}
	MyArray2D HZ_free_arr(2, time_steps);
	double** HZ_free = HZ_free_arr.getdata();
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < time_steps; j++) {
			HZ_free[i][j] = H_free[i][j] * Z_freespace; //setting freespace H-field to 0
		}
	}


	for (int i = 0; i < 100; i++) {
		elements[i].mLength = element_length; //setting all their lengths to equal 2
		elements[i].mRefractive_index = 1; //setting all their refractive indices to 1
		elements[i].mEpsilon = 1* epsilon_freespace;
		elements[i].mMu = 1* mu_freespace;
		elements[i].mZlocal = sqrt(elements[i].mMu/elements[i].mEpsilon);
		elements[i].mYlocal = sqrt(elements[i].mEpsilon/elements[i].mMu);
	};
	for (int z = 0; z < 100; z++) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < time_steps; j++) {
				//setting initial E field and H field //update this so we can propagate a pulse**
				elements[z].mE[i][j] = 1; //setting initial E-field to 1
				elements[z].mZ0H[i][j] = 1; //setting initial H-field to 1 (in correct units)
			}
		}
	}
	//Running the Runge-Kutta method for each element!
	for (int T = 0; T < time_steps; T++) { //running through for each time step
		for (int i = 0; i < 100; i++) { //running through for each element
			//making arrays of size 100 for the Runge-Kutta coefficients
			//k is for node on LHS, l is for node on RHS
			double* k01;
			double* k02;
			double* k11;
			double* k12;
			double* k21;
			double* k22;
			double* k31;
			double* k32;
			double* l01;
			double* l02;
			double* l11;
			double* l12;
			double* l21;
			double* l22;
			double* l31;
			double* l32;
			MyArray1D k01_arr(100);
			MyArray1D k02_arr(100);
			MyArray1D k11_arr(100);
			MyArray1D k12_arr(100);
			MyArray1D k21_arr(100);
			MyArray1D k22_arr(100);
			MyArray1D k31_arr(100);
			MyArray1D k32_arr(100);
			MyArray1D l01_arr(100);
			MyArray1D l02_arr(100);
			MyArray1D l11_arr(100);
			MyArray1D l12_arr(100);
			MyArray1D l21_arr(100);
			MyArray1D l22_arr(100);
			MyArray1D l31_arr(100);
			MyArray1D l32_arr(100);
			k01 = k01_arr.getdata();
			k02 = k02_arr.getdata();
			k11 = k11_arr.getdata();
			k12 = k12_arr.getdata();
			k21 = k21_arr.getdata();
			k22 = k22_arr.getdata();
			k31 = k31_arr.getdata();
			k32 = k32_arr.getdata();
			l01 = l01_arr.getdata();
			l02 = l02_arr.getdata();
			l11 = l11_arr.getdata();
			l12 = l12_arr.getdata();
			l21 = l21_arr.getdata();
			l22 = l22_arr.getdata();
			l31 = l31_arr.getdata();
			l32 = l32_arr.getdata();
			//set an if condition if i=0 or i=99 as we will use periodic b.c.s!!            /free space (absorbing) b.c.s 
			
			//calculating first coefficients
			if (i == 0) { //for first element
				k01[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * (face_m_LHS[0][0] * (elements[99].mE[1][T] - elements[i].mE[0][T]) + face_m_LHS[0][1] * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[99].mZlocal * (face_m_LHS[0][0] * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_LHS[0][1] * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[99].mZlocal + elements[i].mZlocal));
				k02[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * (face_m_RHS[1][0] * (elements[99].mE[1][T] - elements[i].mE[0][T]) + face_m_RHS[1][1] * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[i+1].mZlocal * (face_m_RHS[1][0] * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_RHS[1][1] * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i+1].mZlocal + elements[i].mZlocal));
				l01[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * (face_m_LHS[0][0] * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_LHS[0][1] * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[99].mYlocal * (face_m_LHS[0][0] * (elements[99].mE[1][T] - elements[i].mE[0][T]) + face_m_LHS[0][1] * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[99].mYlocal + elements[i].mYlocal));
				l02[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * (face_m_RHS[1][0] * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_RHS[1][1] * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i + 1].mYlocal * (face_m_RHS[1][0] * elements[99].mE[1][T] - elements[i].mE[0][T] + face_m_RHS[1][1] * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[i+1].mYlocal + elements[i].mYlocal));
			}
			else if (i != 0 && i != 99) {
				k01[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * (face_m_LHS[0][0] * (elements[i-1].mE[1][T] - elements[i].mE[0][T]) + face_m_LHS[0][1] * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[i-1].mZlocal * (face_m_LHS[0][0] * (elements[i-1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_LHS[0][1] * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i-1].mZlocal + elements[i].mZlocal));
				k02[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * (face_m_RHS[1][0] * (elements[i-1].mE[1][T] - elements[i].mE[0][T]) + face_m_RHS[1][1] * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[i + 1].mZlocal * (face_m_RHS[1][0] * (elements[i-1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_RHS[1][1] * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
				l01[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * (face_m_LHS[0][0] * (elements[i-1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_LHS[0][1] * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i-1].mYlocal * (face_m_LHS[0][0] * (elements[i-1].mE[1][T] - elements[i].mE[0][T]) + face_m_LHS[0][1] * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[i-1].mYlocal + elements[i].mYlocal));
				l02[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * (face_m_RHS[1][0] * (elements[i-1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_RHS[1][1] * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i + 1].mYlocal * (face_m_RHS[1][0] * elements[i-1].mE[1][T] - elements[i].mE[0][T] + face_m_RHS[1][1] * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[i + 1].mYlocal + elements[i].mYlocal));

			}
			else { //for last element
				k01[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * (face_m_LHS[0][0] * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + face_m_LHS[0][1] * (elements[0].mE[0][T] - elements[i].mE[1][T])) + elements[i - 1].mZlocal * (face_m_LHS[0][0] * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_LHS[0][1] * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
				k02[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * (face_m_RHS[1][0] * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + face_m_RHS[1][1] * (elements[0].mE[0][T] - elements[i].mE[1][T])) + elements[0].mZlocal * (face_m_RHS[1][0] * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_RHS[1][1] * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[0].mZlocal + elements[i].mZlocal));
				l01[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * (face_m_LHS[0][0] * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_LHS[0][1] * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i - 1].mYlocal * (face_m_LHS[0][0] * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + face_m_LHS[0][1] * (elements[0].mE[0][T] - elements[i].mE[1][T]))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
				l02[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * (face_m_RHS[1][0] * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + face_m_RHS[1][1] * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[0].mYlocal * (face_m_RHS[1][0] * elements[i - 1].mE[1][T] - elements[i].mE[0][T] + face_m_RHS[1][1] * (elements[0].mE[0][T] - elements[i].mE[1][T]))) / (elements[0].mYlocal + elements[i].mYlocal));
			} 

			//calculating second coefficients
			if (i == 0) { //for first element
				k11[i] = (dt/2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i]/2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i]/2)) + (alpha * (face_m_LHS[0][0] * ((elements[99].mE[1][T] + k02[99]/2) - (elements[i].mE[0][T] + k01[i]/2)) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k01[i + 1]/2) - (elements[i].mE[1][T] + k02[i]/2))) + elements[99].mZlocal * (face_m_LHS[0][0] * ((elements[99].mZ0H[1][T] + l02[99]/2) - (elements[i].mZ0H[0][T] + l01[i]/2)) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l01[i + 1]/2) - (elements[i].mZ0H[1][T] + l02[i]/2)))) / (elements[99].mZlocal + elements[i].mZlocal));
				k12[i] = (dt/2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i]/2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i]/2)) + (alpha * (face_m_RHS[1][0] * ((elements[99].mE[1][T] + k02[99]/2) - (elements[i].mE[0][T] + k01[i]/2)) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k01[i + 1]/2) - (elements[i].mE[1][T] + k02[i]/2))) + elements[i + 1].mZlocal * (face_m_RHS[1][0] * ((elements[99].mZ0H[1][T] + l02[99]/2) - (elements[i].mZ0H[0][T] + l01[i]/2)) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l01[i + 1]/2) - (elements[i].mZ0H[1][T] + l02[i]/2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
				l11[i] = (dt/2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i]/2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i]/2)) + (alpha * (face_m_LHS[0][0] * ((elements[99].mZ0H[1][T] + l02[99]/2) - (elements[i].mZ0H[0][T] + l01[i]/2)) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l01[i + 1]/2) - (elements[i].mZ0H[1][T] + l02[i]/2))) - elements[99].mYlocal * (face_m_LHS[0][0] * ((elements[99].mE[1][T] + k02[99]/2) - (elements[i].mE[0][T] + k01[i]/2)) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k01[i + 1]/2) - (elements[i].mE[1][T] + k02[i]/2)))) / (elements[99].mYlocal + elements[i].mYlocal));
				l12[i] = (dt/2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i]/2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i]/2)) + (alpha * (face_m_RHS[1][0] * ((elements[99].mZ0H[1][T] + l02[99]/2) - (elements[i].mZ0H[0][T] + l01[i]/2)) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l01[i + 1]/2) - (elements[i].mZ0H[1][T] + l02[i]/2))) - elements[i + 1].mYlocal * (face_m_RHS[1][0] * ((elements[99].mE[1][T] + k02[99]/2) - (elements[i].mE[0][T] + k01[i]/2)) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k01[i + 1]/2) - (elements[i].mE[1][T] + k02[i]/2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
			}
			else if (i != 0 && i != 99) {
				k11[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i - 1].mZlocal * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
				k12[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i + 1].mZlocal * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
				l11[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i - 1].mYlocal * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
				l12[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i + 1].mYlocal * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
			}
			else {
				k11[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + face_m_LHS[0][1] * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i - 1].mZlocal * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + face_m_LHS[0][1] * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
				k12[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + face_m_RHS[1][1] * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[0].mZlocal * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + face_m_RHS[1][1] * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[0].mZlocal + elements[i].mZlocal));
				l11[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + face_m_LHS[0][1] * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i - 1].mYlocal * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + face_m_LHS[0][1] * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
				l12[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + face_m_RHS[1][1] * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[0].mYlocal * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + face_m_RHS[1][1] * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[0].mYlocal + elements[i].mYlocal));
			}

			//calculating third coefficients
			if (i == 0) { //for first element
				k21[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[99].mZlocal * (face_m_LHS[0][0] * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[99].mZlocal + elements[i].mZlocal));
				k22[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i + 1].mZlocal * (face_m_RHS[1][0] * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
				l21[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[99].mYlocal * (face_m_LHS[0][0] * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[99].mYlocal + elements[i].mYlocal));
				l22[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i + 1].mYlocal * (face_m_RHS[1][0] * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
			}
			else if (i != 0 && i != 99) {
				k21[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i - 1].mZlocal * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
				k22[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i + 1].mZlocal * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
				l21[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i - 1].mYlocal * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
				l22[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i + 1].mYlocal * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
			}
			else {
				k21[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_LHS[0][1] * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i - 1].mZlocal * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_LHS[0][1] * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
				k22[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_RHS[1][1] * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[0].mZlocal * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_RHS[1][1] * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[0].mZlocal + elements[i].mZlocal));
				l21[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_LHS[0][1] * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i - 1].mYlocal * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_LHS[0][1] * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
				l22[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + face_m_RHS[1][1] * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[0].mYlocal * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + face_m_RHS[1][1] * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[0].mYlocal + elements[i].mYlocal));
			}

			//calculating fourth coefficients
			if (i == 0) { //for first element
				k31[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * (face_m_LHS[0][0] * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[99].mZlocal * (face_m_LHS[0][0] * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[99].mZlocal + elements[i].mZlocal));
				k32[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * (face_m_RHS[1][0] * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[i + 1].mZlocal * (face_m_RHS[1][0] * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
				l31[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * (face_m_LHS[0][0] * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[99].mYlocal * (face_m_LHS[0][0] * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[99].mYlocal + elements[i].mYlocal));
				l32[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * (face_m_RHS[1][0] * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i + 1].mYlocal * (face_m_RHS[1][0] * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
			}
			else if (i != 0 && i != 99) {
				k31[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[i - 1].mZlocal * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
				k32[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[i + 1].mZlocal * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
				l31[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_LHS[0][1] * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i - 1].mYlocal * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + face_m_LHS[0][1] * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
				l32[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_RHS[1][1] * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i + 1].mYlocal * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + face_m_RHS[1][1] * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
			}
			else {
				k31[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + face_m_LHS[0][1] * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i]))) + elements[i - 1].mZlocal * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_LHS[0][1] * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
				k32[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + face_m_RHS[1][1] * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i]))) + elements[0].mZlocal * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_RHS[1][1] * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[0].mZlocal + elements[i].mZlocal));
				l31[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * (face_m_LHS[0][0] * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_LHS[0][1] * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i - 1].mYlocal * (face_m_LHS[0][0] * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + face_m_LHS[0][1] * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
				l32[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * (face_m_RHS[1][0] * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + face_m_RHS[1][1] * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[0].mYlocal * (face_m_RHS[1][0] * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + face_m_RHS[1][1] * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i])))) / (elements[0].mYlocal + elements[i].mYlocal));
			}


			//updating E and H fields
			elements[i].mE[0][T + 1] = elements[i].mE[0][T] + ((k01[i] + 2 * k11[i] + 2 * k21[i] + k31[i]) / 6);
			elements[i].mE[1][T + 1] = elements[i].mE[1][T] + ((k02[i] + 2 * k12[i] + 2 * k22[i] + k32[i]) / 6);
			elements[i].mZ0H[0][T + 1] = elements[i].mZ0H[0][T] + ((l01[i] + 2 * l11[i] + 2 * l21[i] + l31[i]) / 6);
			elements[i].mZ0H[1][T + 1] = elements[i].mZ0H[1][T] + ((l02[i] + 2 * l12[i] + 2 * l22[i] + l32[i]) / 6);
		
		}

	}
	



	//save output in text files so we can access it!
	//saving E fields
	std::ofstream MyFile1("E_field_RK.txt");
	for (int i = 0; i < 100; i++) {
		for (int m = 0; m < 2; m++) {
			for (int z = 0; z < time_steps; z++) {
				MyFile1 << elements[i].mE[m][z] << " ";
				std::cout << elements[i].mE[m][z];
			}
			MyFile1 << "" << "\n";
		}
	}
	MyFile1.close();
	//saving H fields
	std::ofstream MyFile2("H_field_RK.txt");
	for (int i = 0; i < 100; i++) {
		for (int m = 0; m < 2; m++) {
			for (int z = 0; z < time_steps; z++) {
				MyFile2 << elements[i].mZ0H[m][z] << " ";
			}
			MyFile2 << "" << "\n";
		}
	}
	MyFile2.close();
}*/
																																																												  
																																																												  
																																																												  
																																																												  
																																																												  
																																																												  
																																																												  //for each time loop need to run for all elements
//create elements inside this function
//create freespace parameters within loop
//}
//for each element need to calculate for LHS then RHS node
//all with refractive index of 1 (so is air)

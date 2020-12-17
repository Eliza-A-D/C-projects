#include <iostream>
#include "myarrays.h"
#include <cmath>
#include <fstream>

int main()
{
	double alpha = 1; //upwind parameter //start with alpha = 1 as is the most stable!
	//alpha = 1; //pure upwind flux
    //formulating basis vectors

    //mesh size
    double h;
    //mesh size is half of wavelength
    double lambda;
    lambda = 600;
    h = 0.5 * lambda;
    double element_size;
    element_size = 2.0;
    double element_length = 2.0;

    //double time_steps = 20;
    int poly_order = 1; //we will normally have this as 1

   // double c = 299792458; //speed of light

    //time step:
    //double dt = element_size / (poly_order * poly_order);
    //double dt = element_size * 0.01 / c; //start off with something small so it's certain to be stable; the larger the time step the less likely it is to be stable

    //declaring our functions for inversing matrices:
    void inverting_matrix_2x2(double** arr, double** inv_arr);
    void inverting_matrix_3x3(double** arr, double** inv_arr);
    void inverting_matrix_4x4(double** arr, double** inv_arr); //using adjugate matrix
    void transpose_matrix_2x2(double** arr, double** tran_arr);
    void transpose_matrix_3x3(double** arr, double** tran_arr);
    //void updating_Runge_Kutta_coefficients(double** mass_matrix, double** stiffness_matrix, double** face_matrix_LHS, double** face_matrix_RHS, double** inverse_mass_matrix, double alpha, double element_length);

    //making mass matrix:
    MyArray2D mass_m_arr(2, 2);
    double** mass_m = mass_m_arr.getdata();
    for (int i = 0; i < mass_m_arr.getrowdim(); i++) {
        for (int j = 0; j < mass_m_arr.getcoldim(); j++) {
            if (i == j) {
                mass_m[i][j] = -element_size/3.0;
            }
            else {
                mass_m[i][j] = element_size/6.0;
            }
        }
    }

    //making inverse mass matrix:
    MyArray2D m_inv_arr(2, 2);
    double** m_inv = m_inv_arr.getdata();
    inverting_matrix_2x2(mass_m, m_inv); //inverting e_mass_m to e_m_inv

    //making stiffness matrix:
    MyArray2D stiff_m_arr(2, 2);
    double** stiff_m = stiff_m_arr.getdata();
    for (int i = 0; i < stiff_m_arr.getrowdim(); i++) {
        for (int j = 0; j < stiff_m_arr.getcoldim(); j++) {
            stiff_m[0][0] = -1.0 / 2.0;
            stiff_m[1][0] = -1.0 / 2.0;
            stiff_m[0][1] = 1.0 / 2.0;
            stiff_m[1][1] = 1.0 / 2.0;
        }
    }

    //making face mass matrix:
    MyArray2D face_arr_LHS(2, 2);
    double** face_m_LHS = face_arr_LHS.getdata();
    face_m_LHS[0][0] = 1.0;
    face_m_LHS[1][0] = 0.0;
    face_m_LHS[0][1] = 0.0;
    face_m_LHS[1][1] = 0.0;
    MyArray2D face_arr_RHS(2, 2);
    double** face_m_RHS = face_arr_RHS.getdata();
    face_m_RHS[0][0] = 0.0;
    face_m_RHS[1][0] = 0.0;
    face_m_RHS[0][1] = 0.0;
    face_m_RHS[1][1] = 1.0;

    //conditions of free space (boundary conditions)
    //double epsilon_freespace = 8.85e-12;
    //double mu_freespace = 12.57e-7;
    //double Z_freespace = sqrt(mu_freespace / epsilon_freespace);
    //double Y_freespace = sqrt(epsilon_freespace / mu_freespace);
    //making freespace E-field and H-field 2D arrays
    //MyArray2D E_free_arr(2, time_steps);
    //double** E_free = E_free_arr.getdata();
    //MyArray2D H_free_arr(2, time_steps);
    //double** H_free = H_free_arr.getdata();
    //for (int i = 0; i < 2; i++) {
        //for (int j = 0; j < time_steps; j++) {
            //E_free[i][j] = 0; //setting freespace E-field to 0
            //H_free[i][j] = 1; //setting freespace H-field to 0
        //}
    //}
    //MyArray2D HZ_free_arr(2, time_steps);
    //double** HZ_free = HZ_free_arr.getdata();
    //for (int i = 0; i < 2; i++) {
        //for (int j = 0; j < time_steps; j++) {
            //HZ_free[i][j] = H_free[i][j] * Z_freespace; //setting freespace H-field to 0
        //}
    //}
   // std::cout << HZ_free[0][0];

    //double E_freespace = 0; //assume no electric field in freespace
    //double H_freespace = 0; //assume no magnetic field in freespace
    //double HZ_freespace = H_freespace * Z_freespace; //setting value in freespace


    //Let's make a line of three elements
    //declaring all features of elements to start with
    //Do three elements of different impedance and conductance and on the outside is free space
    //Make element 1
    /*My1DElement ele1(4, 1, 20);
    ele1.mElementNumber = 1;
    ele1.mEpsilon = 2*epsilon_freespace; //pow(1, -10);
    ele1.mMu = 2*mu_freespace; //pow(1, -5);
    ele1.mZlocal = sqrt(ele1.mMu / ele1.mEpsilon);
    ele1.mYlocal = sqrt(ele1.mEpsilon / ele1.mMu);
    //setting initial magnetic and electric fields for both nodes - first order polynomial approximation:
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < time_steps; j++) {
            ele1.mE[i][j] = 1; //setting initial E-field to 1
            ele1.mH[i][j] = 1; //setting initial H-field to 1
        }
    }


    //Make element 2
    My1DElement ele2(3, 1, 20);
    ele2.mElementNumber = 2;
    ele2.mEpsilon = 3*epsilon_freespace;
    ele2.mMu = 3*mu_freespace;
    ele2.mZlocal = sqrt(ele2.mMu / ele2.mEpsilon);
    ele2.mYlocal = sqrt(ele2.mEpsilon / ele2.mMu);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < time_steps; j++) {
            ele2.mE[i][j] = 1; //setting initial E-field to 1
            ele2.mH[i][j] = 1; //setting initial H-field to 1
        }
    }

    //Make element 3
    My1DElement ele3(4, 1, 20);
    ele3.mElementNumber = 3;
    ele3.mEpsilon = 4*epsilon_freespace;
    ele3.mMu = 4*mu_freespace;
    ele3.mZlocal = sqrt(ele3.mMu / ele3.mEpsilon);
    ele3.mYlocal = sqrt(ele3.mEpsilon / ele3.mMu);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < time_steps; j++) {
            ele3.mE[i][j] = 1; //setting initial E-field to 1
            ele3.mH[i][j] = 1; //setting initial H-field to 1
        }
    } 

    */
    //My1DElement elements[100];
    //My1DElement le;
    //double** ele = elements.getdata();

    //for (int i = 0; i < 10; i++) {
      //  elements[i].getdata();
   // };

    //setting neighbouring impedance and conductance conditions:
    //defining our node 1 to be on the left and our node 2 to be on the right
    //ele1.mZ_1 = Z_freespace;
    //ele1.mZ_2 = ele2.mZlocal;
    //ele1.mY_1 = Y_freespace;
    //ele1.mY_2 = ele2.mYlocal;
    //ele2.mZ_1 = ele1.mZlocal;
    //ele2.mZ_2 = ele3.mZlocal;
    //ele2.mY_1 = ele1.mYlocal;
    //ele2.mY_2 = ele3.mYlocal;
    //ele3.mZ_1 = ele2.mZlocal;
    //ele3.mZ_2 = Z_freespace;
    //ele3.mY_1 = ele2.mYlocal;
    //ele3.mY_2 = Y_freespace;

    //if values don't exist for that entry enter zero
    //declaring function
    //void updating_Runge_Kutta_coefficients(double** mass_matrix, double** stiffness_matrix, double** face_matrix, double** inverse_mass_matrix, double** element_E_field, double** element_H_field, double** before_element_E_field, double** before_element_H_field, double** after_element_E_field, double** after_element_H_field, double element_epsilon, double element_mu, double element_impedance, double element_conductance, double before_element_impedance, double before_element_conductance, double after_element_impedance, double after_element_conductance, double alpha); //will calculate all our updated fields!

    //calling function
    //for element 1
    //updating_Runge_Kutta_coefficients(mass_m, stiff_m, face_m, m_inv, ele1.mE, ele1.mH, E_free, H_free, ele2.mE, ele2.mH, ele1.mEpsilon, ele1.mMu, ele1.mZlocal, ele1.mYlocal, Z_freespace, Y_freespace, ele2.mZlocal, ele2.mYlocal, alpha);
    //for element 2
    //updating_Runge_Kutta_coefficients(mass_m, stiff_m, face_m, m_inv, ele2.mE, ele2.mH, ele1.mE, ele1.mH, ele3.mE, ele3.mH, ele2.mEpsilon, ele2.mMu, ele2.mZlocal, ele2.mYlocal, ele1.mZlocal, ele1.mYlocal, ele3.mZlocal, ele3.mYlocal, alpha);
    //for element 3
    //updating_Runge_Kutta_coefficients(mass_m, stiff_m, face_m, m_inv, ele3.mE, ele3.mH, ele2.mE, ele2.mH, E_free, H_free, ele3.mEpsilon, ele3.mMu, ele3.mZlocal, ele3.mYlocal, ele2.mZlocal, ele2.mYlocal, Z_freespace, Y_freespace, alpha);



    //MyArray2D Runge_Kutta_coefs_arr(4, 4);
    //double** Runge_Kutta_coefs = Runge_Kutta_coefs_arr.getdata();
    //Implementing Runge-Kutta using functions:
   // updating_Runge_Kutta_coefficients();


    //example function to see if they accept objects:
    //void example();
    //example() {
    //    ele1.mH[0][0] = 5;
    //};
    //example();
    //std::cout << ele1.mH[0][0];

    //example function:
    //void exam_print(double impedance); 
    //exam_print(ele1.mZlocal);

    //determining the total number of elements
    //double total_number_of_elements = My1DElement::totalObjects();
    //std::cout << "Total objects created: " << My1DElement::totalObjects() << "\n";
    //updating_Runge_Kutta_coefficients(mass_m, stiff_m, face_m_LHS, face_m_RHS, m_inv, alpha, 2);






    //Implement Runge-Kutta to calculate fields:
    //for element 1:
    /*for (int T = 0; T < time_steps; T++) {
        double k01, k02, k11, k12, k21, k22, k31, k32, l01, l02, l11, l12, l21, l22, l31, l32;
        //for first set of equations with time derivative of E
        //for first node which is bordering with freespace
        k01 = dt / (ele1.mEpsilon) * (((m_inv[0][0]*stiff_m[0][0] + m_inv[0][1]*stiff_m[1][0])*ele1.mH[0][T] + (m_inv[0][0]*stiff_m[0][1] + m_inv[0][1]*stiff_m[1][1])*ele1.mH[1][T])*ele1.mH[0][T] + (alpha* (E_freespace - ele1.mE[0][T]) + Z_freespace*(H_freespace - ele1.mH[0][T]))/(Z_freespace + ele1.mZlocal));
        //for second node which is bordering with element 2
        k02 = dt / (ele1.mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * ele1.mH[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * ele1.mH[1][T]) * ele1.mH[1][T] + (alpha * (ele2.mE[0][T] - ele1.mE[1][T]) + ele2.mZlocal * (ele2.mH[0][T] - ele1.mH[1][T])) / (ele2.mZlocal + ele1.mZlocal));
        //for second set of equations with time derivative of H
        l01 = dt / (ele1.mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * ele1.mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * ele1.mE[1][T]) + (alpha*(H_freespace - ele1.mH[0][T]) - Y_freespace*(E_freespace - ele1.mE[0][T]))/(Y_freespace + ele1.mYlocal));
        l02 = dt / (ele1.mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * ele1.mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * ele1.mE[1][T]) + (alpha * (ele2.mH[0][T] - ele1.mH[1][T]) - ele2.mYlocal * (ele2.mE[0][T] - ele1.mE[1][T])) / (ele2.mYlocal + ele1.mYlocal));

        //next coefficients in Runge-Kutta Method
        k11 = (dt/2) / (ele1.mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele1.mH[0][T] + l01/2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele1.mH[1][T] + l02/2)) * (ele1.mH[0][T] + l01/2) + (alpha * (E_freespace - (ele1.mE[0][T] + (k01/2))) + Z_freespace * (H_freespace - (ele1.mH[0][T]) + (l01/2))) / (Z_freespace + ele1.mZlocal));
        k12 = (dt/2) / (ele1.mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele1.mH[0][T] + l01/2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele1.mH[1][T] + l02/2)) * (ele1.mH[1][T] + l02/2)  + (alpha * ((ele2.mE[0][T] + k01/2) - (ele1.mE[1][T] + k02/2)) + ele2.mZlocal * ((ele2.mH[0][T] + l01/2) - (ele1.mH[1][T] + l02/2))) / (ele2.mZlocal + ele1.mZlocal));
        l11 = (dt/2) / (ele1.mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele1.mE[0][T] + k01/2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele1.mE[1][T] + k02/2)) + (alpha * (H_freespace - (ele1.mH[0][T] + l01/2)) - Y_freespace * (E_freespace - (ele1.mE[0][T] + k01/2))) / (Y_freespace + ele1.mYlocal));
        l12 = (dt/2) / (ele1.mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele1.mE[0][T] + k01/2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele1.mE[1][T] + k02/2)) + (alpha * ((ele2.mH[0][T] + l01/2) - (ele1.mH[1][T] + l02/2)) - ele2.mYlocal * ((ele2.mE[0][T] + k01/2) - (ele1.mE[1][T] + k02/2))) / (ele2.mYlocal + ele1.mYlocal));

        //next coefficients in Runge-Kutta Method
        k21 = (dt / 2) / (ele1.mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele1.mH[0][T] + l11 / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele1.mH[1][T] + l12 / 2)) * (ele1.mH[0][T] + l11 / 2) + (alpha * (E_freespace - (ele1.mE[0][T] + (k11 / 2))) + Z_freespace * (H_freespace - (ele1.mH[0][T]) + (l11 / 2))) / (Z_freespace + ele1.mZlocal));
        k22 = (dt / 2) / (ele1.mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele1.mH[0][T] + l11 / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele1.mH[1][T] + l12 / 2)) * (ele1.mH[1][T] + l12 / 2) + (alpha * ((ele2.mE[0][T] + k11 / 2) - (ele1.mE[1][T] + k12 / 2)) + ele2.mZlocal * ((ele2.mH[0][T] + l11 / 2) - (ele1.mH[1][T] + l12 / 2))) / (ele2.mZlocal + ele1.mZlocal));
        l21 = (dt / 2) / (ele1.mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele1.mE[0][T] + k11 / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele1.mE[1][T] + k12 / 2)) + (alpha * (H_freespace - (ele1.mH[0][T] + l11 / 2)) - Y_freespace * (E_freespace - (ele1.mE[0][T] + k11 / 2))) / (Y_freespace + ele1.mYlocal));
        l22 = (dt / 2) / (ele1.mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele1.mE[0][T] + k11 / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele1.mE[1][T] + k12 / 2)) + (alpha * ((ele2.mH[0][T] + l11 / 2) - (ele1.mH[1][T] + l12 / 2)) - ele2.mYlocal * ((ele2.mE[0][T] + k11 / 2) - (ele1.mE[1][T] + k12 / 2))) / (ele2.mYlocal + ele1.mYlocal));

        //last coefficients in Runge-Kutta Method
        k31 = dt / (ele1.mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele1.mH[0][T] + l21) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele1.mH[1][T] + l22 / 2)) * (ele1.mH[0][T] + l21) + (alpha * (E_freespace - (ele1.mE[0][T] + (k21))) + Z_freespace * (H_freespace - (ele1.mH[0][T]) + (l21))) / (Z_freespace + ele1.mZlocal));
        k32 = dt / (ele1.mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele1.mH[0][T] + l21) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele1.mH[1][T] + l22 / 2)) * (ele1.mH[1][T] + l22) + (alpha * ((ele2.mE[0][T] + k21) - (ele1.mE[1][T] + k22)) + ele2.mZlocal * ((ele2.mH[0][T] + l21) - (ele1.mH[1][T] + l22))) / (ele2.mZlocal + ele1.mZlocal));
        l31 = dt / (ele1.mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele1.mE[0][T] + k21) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele1.mE[1][T] + k22 / 2)) + (alpha * (H_freespace - (ele1.mH[0][T] + l21)) - Y_freespace * (E_freespace - (ele1.mE[0][T] + k21))) / (Y_freespace + ele1.mYlocal));
        l32 = dt / (ele1.mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele1.mE[0][T] + k21) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele1.mE[1][T] + k22 / 2)) + (alpha * ((ele2.mH[0][T] + l21) - (ele1.mH[1][T] + l22)) - ele2.mYlocal * ((ele2.mE[0][T] + k21) - (ele1.mE[1][T] + k22))) / (ele2.mYlocal + ele1.mYlocal));

        ele1.mE[0][T + 1] = ele1.mE[0][T] + ((k01 + 2 * k11 + 2 * k21 + k31) / 6);
        ele1.mE[1][T + 1] = ele1.mE[1][T] + ((k02 + 2 * k12 + 2 * k22 + k32) / 6);
        ele1.mH[0][T + 1] = ele1.mH[0][T] + ((l01 + 2 * l11 + 2 * l21 + l31) / 6);
        ele1.mH[1][T + 1] = ele1.mH[1][T] + ((l02 + 2 * l12 + 2 * l22 + l32) / 6);
    }
    //for element 2
    for (int T = 0; T < time_steps; T++) {
        double k01, k02, k11, k12, k21, k22, k31, k32, l01, l02, l11, l12, l21, l22, l31, l32;
        //for first set of equations with time derivative of E
        //for first node which is bordering with freespace
        k01 = dt / (ele2.mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * ele2.mH[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * ele2.mH[1][T]) * ele2.mH[0][T] + (alpha * (ele1.mE[1][T] - ele2.mE[0][T]) + ele1.mZlocal * (ele1.mH[1][T] - ele2.mH[0][T])) / (ele1.mZlocal + ele2.mZlocal));
        //for second node which is bordering with element 2
        k02 = dt / (ele2.mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * ele2.mH[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * ele2.mH[1][T]) * ele2.mH[1][T] + (alpha * (ele3.mE[0][T] - ele2.mE[1][T]) + ele3.mZlocal * (ele3.mH[0][T] - ele2.mH[1][T])) / (ele3.mZlocal + ele2.mZlocal));
        //for second set of equations with time derivative of H
        l01 = dt / (ele2.mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * ele2.mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * ele2.mE[1][T]) + (alpha * (ele1.mH[1][T] - ele2.mH[0][T]) - ele1.mYlocal * (ele1.mE[1][T] - ele2.mE[0][T])) / (ele1.mYlocal + ele2.mYlocal));
        l02 = dt / (ele2.mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * ele2.mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * ele2.mE[1][T]) + (alpha * (ele3.mH[0][T] - ele2.mH[1][T]) - ele3.mYlocal * (ele3.mE[0][T] - ele2.mE[1][T])) / (ele3.mYlocal + ele2.mYlocal));

        //next coefficients in Runge-Kutta Method
        k11 = (dt / 2) / (ele2.mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele2.mH[0][T] + l01 / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele2.mH[1][T] + l02 / 2)) * (ele2.mH[0][T] + l01 / 2) + (alpha * ((ele1.mE[1][T] + k02/2) - (ele2.mE[0][T] + (k01 / 2))) + ele1.mZlocal * ((ele1.mH[1][T] + l02/2) - (ele2.mH[0][T]) + (l01 / 2))) / (ele1.mZlocal + ele2.mZlocal));
        k12 = (dt / 2) / (ele2.mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele2.mH[0][T] + l01 / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele2.mH[1][T] + l02 / 2)) * (ele2.mH[1][T] + l02 / 2) + (alpha * ((ele3.mE[0][T] + k01 / 2) - (ele2.mE[1][T] + k02 / 2)) + ele3.mZlocal * ((ele3.mH[0][T] + l01 / 2) - (ele2.mH[1][T] + l02 / 2))) / (ele3.mZlocal + ele2.mZlocal));
        l11 = (dt / 2) / (ele2.mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele2.mE[0][T] + k01 / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele2.mE[1][T] + k02 / 2)) + (alpha * ((ele1.mH[1][T] + l02 / 2) - (ele2.mH[0][T] + l01 / 2)) - ele1.mYlocal * ((ele1.mE[1][T] + k02/2) - (ele2.mE[0][T] + k01 / 2))) / (ele1.mYlocal + ele2.mYlocal));
        l12 = (dt / 2) / (ele2.mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele2.mE[0][T] + k01 / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele2.mE[1][T] + k02 / 2)) + (alpha * ((ele3.mH[0][T] + l01 / 2) - (ele2.mH[1][T] + l02 / 2)) - ele3.mYlocal * ((ele3.mE[0][T] + k01 / 2) - (ele2.mE[1][T] + k02 / 2))) / (ele3.mYlocal + ele2.mYlocal));

        //next coefficients in Runge-Kutta Method
        k21 = (dt / 2) / (ele2.mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele2.mH[0][T] + l11 / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele2.mH[1][T] + l12 / 2)) * (ele2.mH[0][T] + l11 / 2) + (alpha * ((ele2.mE[1][T] + k12/2) - (ele2.mE[0][T] + (k11 / 2))) + ele1.mZlocal * ((ele1.mH[1][T] + l12/2) - (ele2.mH[0][T]) + (l11 / 2))) / (ele1.mZlocal + ele2.mZlocal));
        k22 = (dt / 2) / (ele2.mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele2.mH[0][T] + l11 / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele2.mH[1][T] + l12 / 2)) * (ele2.mH[1][T] + l12 / 2) + (alpha * ((ele3.mE[0][T] + k11 / 2) - (ele2.mE[1][T] + k12 / 2)) + ele3.mZlocal * ((ele3.mH[0][T] + l11 / 2) - (ele2.mH[1][T] + l12 / 2))) / (ele3.mZlocal + ele2.mZlocal));
        l21 = (dt / 2) / (ele2.mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele2.mE[0][T] + k11 / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele2.mE[1][T] + k12 / 2)) + (alpha * ((ele1.mH[1][T] + l12 /2) - (ele2.mH[0][T] + l11 / 2)) - ele1.mYlocal * ((ele1.mE[1][T] + k12/2) - (ele2.mE[0][T] + k11 / 2))) / (ele1.mYlocal + ele2.mYlocal));
        l22 = (dt / 2) / (ele2.mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele2.mE[0][T] + k11 / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele2.mE[1][T] + k12 / 2)) + (alpha * ((ele3.mH[0][T] + l11 / 2) - (ele2.mH[1][T] + l12 / 2)) - ele3.mYlocal * ((ele3.mE[0][T] + k11 / 2) - (ele2.mE[1][T] + k12 / 2))) / (ele3.mYlocal + ele2.mYlocal));

        //last coefficients in Runge-Kutta Method
        k31 = dt / (ele2.mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele2.mH[0][T] + l21) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele2.mH[1][T] + l22)) * (ele2.mH[0][T] + l21) + (alpha * ((ele2.mE[1][T] + k12) - (ele2.mE[0][T] + (k21))) + ele1.mZlocal * ((ele1.mH[1][T] + l12) - (ele2.mH[0][T]) + (l21))) / (ele1.mZlocal + ele2.mZlocal));
        k32 = dt / (ele2.mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele2.mH[0][T] + l21) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele2.mH[1][T] + l22)) * (ele2.mH[1][T] + l22) + (alpha * ((ele3.mE[0][T] + k21) - (ele2.mE[1][T] + k22)) + ele3.mZlocal * ((ele3.mH[0][T] + l21) - (ele2.mH[1][T] + l22))) / (ele3.mZlocal + ele2.mZlocal));
        l31 = dt / (ele2.mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (ele2.mE[0][T] + k21) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (ele2.mE[1][T] + k22)) + (alpha * ((ele1.mH[1][T] + l22) - (ele2.mH[0][T] + l21)) - ele1.mYlocal * ((ele1.mE[1][T] + k22) - (ele2.mE[0][T] + k21))) / (ele1.mYlocal + ele2.mYlocal));
        l32 = dt / (ele2.mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (ele2.mE[0][T] + k21) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (ele2.mE[1][T] + k22)) + (alpha * ((ele3.mH[0][T] + l21) - (ele2.mH[1][T] + l22)) - ele3.mYlocal * ((ele3.mE[0][T] + k21) - (ele2.mE[1][T] + k22))) / (ele3.mYlocal + ele2.mYlocal));

        ele2.mE[0][T + 1] = ele2.mE[0][T] + ((k01 + 2 * k11 + 2 * k21 + k31) / 6);
        ele2.mE[1][T + 1] = ele2.mE[1][T] + ((k02 + 2 * k12 + 2 * k22 + k32) / 6);
        ele2.mH[0][T + 1] = ele2.mH[0][T] + ((l01 + 2 * l11 + 2 * l21 + l31) / 6);
        ele2.mH[1][T + 1] = ele2.mH[1][T] + ((l02 + 2 * l12 + 2 * l22 + l32) / 6);
    }
    //////add in element 3!!! */

    /*
    //saving output
    std::ofstream MyFile1("E_field_RK_ele1.txt");
    for (int m = 0; m < 2; m++) {
        for (int z = 0; z < time_steps; z++) {
            MyFile1 << ele1.mE[m][z] << " ";
        }
        MyFile1 << "" << "\n";
    }
    MyFile1.close();


    std::ofstream MyFile2("H_field_RK_ele1.txt");
    for (int m = 0; m < 2; m++) {
        for (int z = 0; z < time_steps; z++) {
            MyFile2 << ele1.mH[m][z] << " ";
        }
        MyFile2 << "" << "\n";
    }
    MyFile2.close();

    std::ofstream MyFile3("E_field_RK_ele2.txt");
    for (int m = 0; m < 2; m++) {
        for (int z = 0; z < time_steps; z++) {
            MyFile3 << ele2.mE[m][z] << " ";
        }
        MyFile3 << "" << "\n";
    }
    MyFile3.close();


    std::ofstream MyFile4("H_field_RK_ele2.txt");
    for (int m = 0; m < 2; m++) {
        for (int z = 0; z < time_steps; z++) {
            MyFile4 << ele2.mH[m][z] << " ";
        }
        MyFile4 << "" << "\n";
    }
    MyFile4.close();

    std::ofstream MyFile5("E_field_RK_ele3.txt");
    for (int m = 0; m < 2; m++) {
        for (int z = 0; z < time_steps; z++) {
            MyFile5 << ele3.mE[m][z] << " ";
        }
        MyFile5 << "" << "\n";
    }
    MyFile5.close();


    std::ofstream MyFile6("H_field_RK_ele3.txt");
    for (int m = 0; m < 2; m++) {
        for (int z = 0; z < time_steps; z++) {
            MyFile6 << ele3.mH[m][z] << " ";
        }
        MyFile6 << "" << "\n";
    }
    MyFile6.close(); */


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
   // elements[2].mLength = 1;
    //std::cout << elements[2].mLength;
    /*for (int i = 0; i < 2; i++) {
        elements[i].mLength = 1;
        std::cout << elements[0].mLength;
        std::cout << elements[1].mLength;
        std::cout << elements[i].mLength;
    } */


    for (int i = 0; i < 100; i++) {
        elements[i].mLength = element_length; //setting all their lengths to equal 2
        elements[i].mRefractive_index = 1; //setting all their refractive indices to 1
        elements[i].mEpsilon = 1 * epsilon_freespace;
        elements[i].mMu = 1 * mu_freespace;
        elements[i].mZlocal = sqrt(elements[i].mMu / elements[i].mEpsilon);
        elements[i].mYlocal = sqrt(elements[i].mEpsilon / elements[i].mMu);
    };

    
    
    for (int z = 0; z < 100; z++) {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < time_steps; j++) {
                if (z == 0) {
                    elements[z].mE[i][j] = 1;
                    elements[z].mH[i][j] = 1;
                    elements[z].mZ0H[i][j] = 1;
                }
                else {
                    //setting initial E field and H field //update this so we can propagate a pulse**
                    elements[z].mE[i][j] = 0; //setting initial E-field to 1
                    elements[z].mH[i][j] = 0;
                    //elements[z].mZ0H[i][j] = 1*Z_freespace; //setting initial H-field to 1 (in correct units)
                    elements[z].mZ0H[i][j] = 0; //let's start them off the same size
                }
            }
        }
    }

    //std::cout << dt;
    //std::cout << elements[0].mEpsilon;
    //std::cout << stiff_m[0][0];
    //std::cout << elements[99].mZlocal << "\n";
    //std::cout << elements[99].mZlocal + elements[0].mZlocal << "\n";
    //std::cout << dt / (elements[0].mEpsilon) * ((((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[0].mZ0H[0][0] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[0].mZ0H[1][0])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mE[1][0] - elements[0].mE[0][0]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[1].mE[0][0] - elements[0].mE[1][0])) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mZ0H[1][0] - elements[0].mZ0H[0][0]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[1].mZ0H[0][0] - elements[0].mZ0H[1][0]))) / (elements[99].mZlocal + elements[0].mZlocal));

    //std::cout << dt / (elements[0].mEpsilon) * ((((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[0].mZ0H[0][0])));
    //std::cout << (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[0].mZ0H[1][0];
    //std::cout << dt / (elements[0].mEpsilon) * ((((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[0].mZ0H[0][0] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[0].mZ0H[1][0])));
    //std::cout << (m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0])* elements[0].mZ0H[0][0] << "\n";
    //std::cout << (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1])* elements[0].mZ0H[1][0] << "\n";
    //std::cout << elements[99].mE[1][0] - elements[0].mE[0][0];
   //std::cout << dt / (elements[0].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[0].mZ0H[0][0] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[0].mZ0H[1][0]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mE[1][0] - elements[0].mE[0][0]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[1].mE[0][0] - elements[0].mE[1][0])) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mZ0H[1][0] - elements[0].mZ0H[0][0]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[1].mZ0H[0][0] - elements[0].mZ0H[1][0]))) / (elements[99].mZlocal + elements[0].mZlocal));
    //std::cout << dt / (elements[0].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[0].mZ0H[0][0] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[0].mZ0H[1][0]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[99].mE[1][0] - elements[0].mE[0][0]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[1].mE[0][0] - elements[0].mE[1][0])) + elements[1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[99].mZ0H[1][0] - elements[0].mZ0H[0][0]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[1].mZ0H[0][0] - elements[0].mZ0H[1][0]))) / (elements[1].mZlocal + elements[0].mZlocal));
    //std::cout << dt / (elements[0].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[0].mE[0][0] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[0].mE[1][0]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mZ0H[1][0] - elements[0].mZ0H[0][0]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[1].mZ0H[0][0] - elements[0].mZ0H[1][0])) - elements[99].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mE[1][0] - elements[0].mE[0][0]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[1].mE[0][0] - elements[0].mE[1][0]))) / (elements[99].mYlocal + elements[0].mYlocal));
    //std::cout << dt / (elements[0].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[0].mE[0][0] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[0].mE[1][0]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[99].mZ0H[1][0] - elements[0].mZ0H[0][0]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[1].mZ0H[0][0] - elements[0].mZ0H[1][0])) - elements[1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * elements[99].mE[1][0] - elements[0].mE[0][0] + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[1].mE[0][0] - elements[0].mE[1][0]))) / (elements[1].mYlocal + elements[0].mYlocal));
    //std::cout << (dt / 2) / (elements[0].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[0].mZ0H[0][0] + 0.0199488 / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[0].mZ0H[1][0] + 11.3373 / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][0] + 0.0299762 / 2) - (elements[0].mE[0][0] + 7.55817 / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[1].mE[0][0] + 7.55817 / 2) - (elements[0].mE[1][0] + 0.0299762 / 2))) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][0] + 11.3373 / 2) - (elements[0].mZ0H[0][0] + 0.0199488 / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[1].mZ0H[0][0] + 0.0199488 / 2) - (elements[0].mZ0H[1][0] + 11.3373 / 2)))) / (elements[99].mZlocal + elements[0].mZlocal));


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

//Running the Runge-Kutta method for each element!
    for (int T = 0; T < time_steps; T++) { //running through for each time step
   //for loop for first coefficients:
        for (int i = 0; i < 100; i++) {
            if (i == 0) { //for first element
                k01[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[99].mZlocal + elements[i].mZlocal));
                //std::cout << k01[0];
                k02[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[99].mE[1][T] - elements[i].mE[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                //std::cout << k02[0];
                l01[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[99].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[99].mYlocal + elements[i].mYlocal));
                //std::cout << l01[0];
                l02[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * elements[99].mE[1][T] - elements[i].mE[0][T] + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
                //std::cout << l02[0];
                //std::cout << i << "\n";
            }
            //else if (i != 0 && i != 99) {
            else if (i > 0 && i < 99) {
                k01[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k02[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l01[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l02[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * elements[i - 1].mE[1][T] - elements[i].mE[0][T] + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
                //std::cout << i << "\n";
            }
            //else { //for last element
            else if (i == 99) {
                k01[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[0].mE[0][T] - elements[i].mE[1][T])) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k02[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[0].mE[0][T] - elements[i].mE[1][T])) + elements[0].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[0].mZlocal + elements[i].mZlocal));
                l01[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[0].mE[0][T] - elements[i].mE[1][T]))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l02[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[0].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * elements[i - 1].mE[1][T] - elements[i].mE[0][T] + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[0].mE[0][T] - elements[i].mE[1][T]))) / (elements[0].mYlocal + elements[i].mYlocal));
            }
        }
        //for loop for second coefficients:
        for (int i = 0; i < 100; i++) {
            if (i == 0) { //for first element
                k11[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k02[99] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l02[99] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[99].mZlocal + elements[i].mZlocal));
                k12[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k02[99] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l02[99] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l11[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l02[99] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[99].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k02[99] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[99].mYlocal + elements[i].mYlocal));
                l12[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l02[99] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k02[99] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else if (i != 0 && i != 99) {
            else if (i > 0 && i < 99) {
                k11[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k12[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l11[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l12[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else {
            else if (i == 99) {
                k11[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k12[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[0].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[0].mZlocal + elements[i].mZlocal));
                l11[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l12[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[0].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[0].mYlocal + elements[i].mYlocal));
            }
        }
        //For loop for third coefficients: 
        for (int i = 0; i < 100; i++) {
            if (i == 0) { //for first element
                k21[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[99].mZlocal + elements[i].mZlocal));
                k22[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l21[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[99].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[99].mYlocal + elements[i].mYlocal));
                l22[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else if (i != 0 && i != 99) {
            else if (i > 0 && i < 99) {
                k21[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k22[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l21[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l22[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else {
            else if (i == 99) {
                k21[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k22[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[0].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[0].mZlocal + elements[i].mZlocal));
                l21[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l22[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[0].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[0].mYlocal + elements[i].mYlocal));
            }
        }
        //For loop for fourth coefficients:
        for (int i = 0; i < 100; i++) {
            if (i == 0) { //for first element
                k31[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[99].mZlocal + elements[i].mZlocal));
                k32[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l31[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[99].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[99].mYlocal + elements[i].mYlocal));
                l32[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else if (i != 0 && i != 99) {
            else if (i > 0 && i < 99) {
                k31[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k32[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l31[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l32[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else {
            else if (i == 99) {
                k31[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i]))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k32[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i]))) + elements[0].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[0].mZlocal + elements[i].mZlocal));
                l31[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l32[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[0].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i])))) / (elements[0].mYlocal + elements[i].mYlocal));
                //std::cout << k31[i];
            }
        }
        for (int i = 0; i < 100; i++) {
            //updating E and H fields
            elements[i].mE[0][T + 1] = elements[i].mE[0][T] + ((k01[i] + 2 * k11[i] + 2 * k21[i] + k31[i]) / 6);
            elements[i].mE[1][T + 1] = elements[i].mE[1][T] + ((k02[i] + 2 * k12[i] + 2 * k22[i] + k32[i]) / 6);
            elements[i].mZ0H[0][T + 1] = elements[i].mZ0H[0][T] + ((l01[i] + 2 * l11[i] + 2 * l21[i] + l31[i]) / 6);
            elements[i].mZ0H[1][T + 1] = elements[i].mZ0H[1][T] + ((l02[i] + 2 * l12[i] + 2 * l22[i] + l32[i]) / 6);
        }
    }


                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           
                                           //for (int T = 0; T < 1; T++){
        //for (int i = 0; i < 100; i++) { //running through for each element
        //making arrays of size 100 for the Runge-Kutta coefficients
        //k is for node on LHS, l is for node on RHS
        //set an if condition if i=0 or i=99 as we will use periodic b.c.s!!            /free space (absorbing) b.c.s 
/*
        //calculating first coefficients
            if (i == 0) { //for first element
                k01[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[0][0]*face_m_LHS[0][0] + m_inv[0][1]*face_m_LHS[1][0]) * (elements[99].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0]*face_m_LHS[0][1] + m_inv[0][1]*face_m_LHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[99].mZlocal + elements[i].mZlocal));
                //std::cout << k01[0];
                k02[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[1][0]*face_m_LHS[0][0]+m_inv[1][1]*face_m_RHS[1][0]) * (elements[99].mE[1][T] - elements[i].mE[0][T]) + (m_inv[1][0]*face_m_LHS[0][1]+m_inv[1][1]*face_m_RHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                //std::cout << k02[0];
                l01[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[99].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[99].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[99].mYlocal + elements[i].mYlocal));
                //std::cout << l01[0];
                l02[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[99].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * elements[99].mE[1][T] - elements[i].mE[0][T] + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
                //std::cout << l02[0];
                //std::cout << i << "\n";
            }
            //else if (i != 0 && i != 99) {
            else if (i > 0 && i < 99){
                k01[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k02[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T])) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l01[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l02[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * elements[i - 1].mE[1][T] - elements[i].mE[0][T] + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[i + 1].mE[0][T] - elements[i].mE[1][T]))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
                //std::cout << i << "\n";
            }
            //else { //for last element
            else if (i == 99){
                k01[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[0].mE[0][T] - elements[i].mE[1][T])) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k02[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mZ0H[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mZ0H[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[0].mE[0][T] - elements[i].mE[1][T])) + elements[0].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T]))) / (elements[0].mZlocal + elements[i].mZlocal));
                l01[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * (elements[i - 1].mE[1][T] - elements[i].mE[0][T]) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * (elements[0].mE[0][T] - elements[i].mE[1][T]))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l02[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * elements[i].mE[0][T] + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * elements[i].mE[1][T]) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * (elements[i - 1].mZ0H[1][T] - elements[i].mZ0H[0][T]) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[0].mZ0H[0][T] - elements[i].mZ0H[1][T])) - elements[0].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * elements[i - 1].mE[1][T] - elements[i].mE[0][T] + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * (elements[0].mE[0][T] - elements[i].mE[1][T]))) / (elements[0].mYlocal + elements[i].mYlocal));
            }

            //std::cout << k01[0] << "\n";
            //std::cout << k02[0] << "\n";
            //std::cout << l01[0] << "\n";

        //calculating second coefficients
            if (i == 0) { //for first element
                k11[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k02[99] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l02[99] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[99].mZlocal + elements[i].mZlocal));
                k12[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k02[99] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l02[99] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l11[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l02[99] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[99].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k02[99] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[99].mYlocal + elements[i].mYlocal));
                l12[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l02[99] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k02[99] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else if (i != 0 && i != 99) {
            else if (i > 0 && i < 99) {
                k11[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k12[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l11[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l12[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l01[i + 1] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k01[i + 1] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else {
            else if (i==99){
                k11[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k12[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2))) + elements[0].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2)))) / (elements[0].mZlocal + elements[i].mZlocal));
                l11[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l12[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k01[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k02[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l02[i - 1] / 2) - (elements[i].mZ0H[0][T] + l01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l01[0] / 2) - (elements[i].mZ0H[1][T] + l02[i] / 2))) - elements[0].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k02[i - 1] / 2) - (elements[i].mE[0][T] + k01[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k01[0] / 2) - (elements[i].mE[1][T] + k02[i] / 2)))) / (elements[0].mYlocal + elements[i].mYlocal));
            }

            //std::cout << k01[0] << "\n";
            //std::cout << k02[0] << "\n";
            //std::cout << l01[0] << "\n";

        //calculating third coefficients
            if (i == 0) { //for first element
                k21[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[99].mZlocal + elements[i].mZlocal));
                k22[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l21[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[99].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[99].mYlocal + elements[i].mYlocal));
                l22[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l12[99] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k12[99] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else if (i != 0 && i != 99) {
            else if (i > 0 && i < 99) {
                k21[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k22[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l21[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l22[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l11[i + 1] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k11[i + 1] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else {
            else if (i==99){
                k21[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k22[i] = (dt / 2) / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2))) + elements[0].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2)))) / (elements[0].mZlocal + elements[i].mZlocal));
                l21[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l22[i] = (dt / 2) / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k11[i] / 2) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k12[i] / 2)) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l12[i - 1] / 2) - (elements[i].mZ0H[0][T] + l11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l11[0] / 2) - (elements[i].mZ0H[1][T] + l12[i] / 2))) - elements[0].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k12[i - 1] / 2) - (elements[i].mE[0][T] + k11[i] / 2)) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k11[0] / 2) - (elements[i].mE[1][T] + k12[i] / 2)))) / (elements[0].mYlocal + elements[i].mYlocal));
            }

        //calculating fourth coefficients
            if (i == 0) { //for first element
                k31[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[99].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[99].mZlocal + elements[i].mZlocal));
                k32[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l31[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[99].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[99].mYlocal + elements[i].mYlocal));
                l32[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mZ0H[1][T] + l22[99]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[99].mE[1][T] + k22[99]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else if (i != 0 && i != 99) {
            else if (i > 0 && i < 99) {
                k31[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k32[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i]))) + elements[i + 1].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i + 1].mZlocal + elements[i].mZlocal));
                l31[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l32[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mZ0H[0][T] + l21[i + 1]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i + 1].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[i + 1].mE[0][T] + k21[i + 1]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i + 1].mYlocal + elements[i].mYlocal));
            }
            //else {
            else if (i==99){
                k31[i] = dt / (elements[i].mEpsilon) * (((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i]))) + elements[i - 1].mZlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[i - 1].mZlocal + elements[i].mZlocal));
                k32[i] = dt / (elements[i].mEpsilon) * (((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mZ0H[0][T] + l21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mZ0H[1][T] + l22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i]))) + elements[0].mZlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i])))) / (elements[0].mZlocal + elements[i].mZlocal));
                l31[i] = dt / (elements[i].mMu) * (-((m_inv[0][0] * stiff_m[0][0] + m_inv[0][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[0][0] * stiff_m[0][1] + m_inv[0][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[i - 1].mYlocal * ((m_inv[0][0] * face_m_LHS[0][0] + m_inv[0][1] * face_m_LHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[0][0] * face_m_LHS[0][1] + m_inv[0][1] * face_m_LHS[1][1]) * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i])))) / (elements[i - 1].mYlocal + elements[i].mYlocal));
                l32[i] = dt / (elements[i].mMu) * (-((m_inv[1][0] * stiff_m[0][0] + m_inv[1][1] * stiff_m[1][0]) * (elements[i].mE[0][T] + k21[i]) + (m_inv[1][0] * stiff_m[0][1] + m_inv[1][1] * stiff_m[1][1]) * (elements[i].mE[1][T] + k22[i])) + (alpha * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mZ0H[1][T] + l22[i - 1]) - (elements[i].mZ0H[0][T] + l21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mZ0H[0][T] + l21[0]) - (elements[i].mZ0H[1][T] + l22[i]))) - elements[0].mYlocal * ((m_inv[1][0] * face_m_LHS[0][0] + m_inv[1][1] * face_m_RHS[1][0]) * ((elements[i - 1].mE[1][T] + k22[i - 1]) - (elements[i].mE[0][T] + k21[i])) + (m_inv[1][0] * face_m_LHS[0][1] + m_inv[1][1] * face_m_RHS[1][1]) * ((elements[0].mE[0][T] + k21[0]) - (elements[i].mE[1][T] + k22[i])))) / (elements[0].mYlocal + elements[i].mYlocal));
                //std::cout << k31[i];
            };
            //std::cout << k01[0] << "\n";
            //std::cout << k02[0] << "\n";
            //std::cout << l01[0] << "\n";
            //std::cout << k01[0];
            //std::cout << "k01 " << k01[0] << "\n";
            //std::cout << "k11 " << k11[0] << "\n";
            //std::cout << "k21 " << k21[0] << "\n";
            //std::cout << "k31 " << k31[0] << "\n";
            //std::cout << "Updating coefficient: " << (k01[0] + 2 * k11[0] + 2 * k21[0] + k31[0]) / 6 << "\n";

        //updating E and H fields
            elements[i].mE[0][T + 1] = elements[i].mE[0][T] + ((k01[i] + 2 * k11[i] + 2 * k21[i] + k31[i]) / 6);
            //std::cout << elements[0].mE[0][1];
            elements[i].mE[1][T + 1] = elements[i].mE[1][T] + ((k02[i] + 2 * k12[i] + 2 * k22[i] + k32[i]) / 6);
            elements[i].mZ0H[0][T + 1] = elements[i].mZ0H[0][T] + ((l01[i] + 2 * l11[i] + 2 * l21[i] + l31[i]) / 6);
            elements[i].mZ0H[1][T + 1] = elements[i].mZ0H[1][T] + ((l02[i] + 2 * l12[i] + 2 * l22[i] + l32[i]) / 6);
        }
        
    } */

    //std::cout << "k01 " << k01[0] << "\n";


//save output in text files so we can access it!
//saving E fields
    std::ofstream MyFile1("E_field_RK.txt");
    for (int i = 0; i < 100; i++) {
        for (int m = 0; m < 2; m++) {
            for (int z = 0; z < time_steps; z++) {
                MyFile1 << elements[i].mE[m][z] << " ";
                //std::cout << elements[i].mE[m][z];
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





}
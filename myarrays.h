//making a uniform 1d array
class MyArray1D {
public:
	//declaring constructors
	MyArray1D();
	MyArray1D(int size);
	//declaring destructors
	~MyArray1D();
private: //we will declare all as public for now, is easier
	double* mData;
	int mSize;
public:
	double* getdata() { return mData; } //access function to get mData
	int getsize() { return mSize; } //access function to get mSize

	double operator()(int a) {
		return mData[a];
	};
};

//making an individual 1d element
class My1DElement {
public:
	static int count;
	My1DElement();
	My1DElement(double length, double time_steps); //use if want default order of polynomial(1)/number of nodes (2)
	My1DElement(double length, int order_of_polynomial, double time_steps); //use if want to alter the number of nodes
	~My1DElement();
	//counting the number of elements we produce
	//static int counter = 0;
	//static int getCreated() { return instanceCount; }
	//private:
	double* mData;
	double mRefractive_index;
	double mEpsilon; //electric permitivitty
	double mMu; //electric permeabililty
	double mLength;
	double mElementNumber;
	double mZlocal; //local impedance
	double mYlocal; //local conductance
	double mZ_1; //node 1 boundary condition impedance
	double mZ_2; //node 2 boundary condition impedance
	double mY_1; //node 1 boundary condition conductance
	double mY_2; //node 2 boundary condition conductance
	int mNumber_of_nodes;
	double** mE; //is 2d array of E-field in time, is 2xtime_step as 2 is for coordinate nodes
	double** mH; //is 2d array of H-field in time, is 2xtime_step as 2 is for coordinate nodes
	double** mZ0H; //now it will be in the correct units, we have scaled it to give in same units as E-field
	double mTime_steps;
	double mMid_point; //this will give us the midpoint of the element, which will be super useful for plotting!
	double* getdata() { return mData; } //get data
	double getlength() { return mLength; } //get length of element

	//creating a static member function to get total count:
	static int totalObjects(void) {
		return count;
	}

	//making a function  to get Runge-Kutta coefficients:
};

//making a uniform 2d array
class MyArray2D {
public:
	MyArray2D();
	MyArray2D(int row_dim, int col_dim);
	~MyArray2D();
private:
	double** mData;
	int mrow_dim;
	int mcol_dim;
public:
	double** getdata() { return mData; } //access function to get value mData
	int getrowdim() { return mrow_dim; } //access function to get value mrow_dim
	int getcoldim() { return mcol_dim; } //access function to get value mcol_dim

	double operator()(int a, int b) {
		return mData[a][b];
	};

};





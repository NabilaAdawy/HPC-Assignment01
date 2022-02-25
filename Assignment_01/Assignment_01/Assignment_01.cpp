#include <iostream>
#include <vector>
#include<omp.h>
using std::vector;

// structures of CRS & regular Matrices & a vector.
struct crsMatrix
{
	// Matrix size (N x N)
	int N;
	// Number of non-zero elements
	int NZ;
	// Array of values (size NZ)
	double* Value;
	// Array of column numbers (size NZ)
	int* Col;
	// Array of row indexes (size N + 1)
	int* row_index;
};
struct regMatrix
{
	int N;
	double** P;
};
struct Vec {
	int N;
	double** P;
};

// Allocating Matrices & vector.
void AllocateCRSMatrix(int N, int NZ, crsMatrix& mtx) {
	mtx.N = N;
	mtx.NZ = NZ;
	mtx.Value = new double [NZ];
	mtx.Col = new int [NZ];
	mtx.row_index = new int [N+1];

}
void AllocateRegMatrix(int N, regMatrix& mtx){
	mtx.N = N;
	mtx.P = new double*[N];
	for (int i = 0; i < N; i++) {
		mtx.P[i] = new double[N];
	}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			mtx.P[i][j] = 0;
}
void AllocateVector(int N, Vec& vec) {
	vec.N = N;
	vec.P = new double* [N];
	for (int i = 0; i < N; i++) {
		vec.P[i] = new double[1]();
	}
}

// Deallocating Matrices & vector.
void FreeRegMatrix(regMatrix& mtx) {
	for (int i = 0; i < mtx.N; ++i)
		delete[] mtx.P[i];
	delete[] mtx.P;
}
void FreeCRSMatrix(crsMatrix& mtx) {
	delete[] mtx.Value;
	delete[] mtx.Col;
	delete[] mtx.row_index;
}
void FreeVector(Vec& vec) {
	for (size_t i = 0; i < vec.N; ++i) {
		delete[] vec.P[i];
	}
	delete[] vec.P;
}

// Copying Matrices & vector
void CopyRegMatrix(regMatrix im, regMatrix& om) {
	for (int i = 0; i < im.N; i++)
	{
		for (int j = 0; j < im.N; j++)
		{
			om.P[i][j] = im.P[i][j];
		}
	}
}
void CopyCRSMatrix(crsMatrix im, crsMatrix& om) {
	for (int i = 0; i < im.NZ; i++)
	{
		om.Value[i] = im.Value[i];
		om.Col[i] = im.Col[i];
	}
	for (int i = 0; i < im.N+1; i++)
	{
		om.row_index[i] = im.row_index[i];
	}
}
void CopyVector(int n, Vec& iv, Vec& ov) {
	for (int i = 0; i < n; i++)
	{
			ov.P[i][0] = iv.P[i][0];
	}
}

double next()
{
    return ((double)rand() / (double)RAND_MAX);
}

// Comparing Functions (between 2 vectos, between 2 matrices).
int CompareVectors(Vec& vec1 , Vec& vec2, int n, double& diff)
{
	diff = 0.0;
	for (int i = 0; i < n; i++)
	{
		if (diff < fabs(vec1.P[i] - vec2.P[i]))
		{
			diff = fabs(vec1.P[i] - vec2.P[i]);
		}
	}
	return 0;
}
int CompareMatrices(regMatrix& m1, regMatrix& m2, double& diff) {
	diff = 0.0;
	for (int j = 0; j < m1.N; j++) {
		for (int i = 0; i < m1.N; i++)
		{
			if (diff < fabs(m1.P[j][i] - m2.P[j][i]))
			{
				diff = fabs(m1.P[j][i] - m2.P[j][i]);
			}
		}
	}
	return 0;
}
// Generating Matrices & vector.
void GenerateRegularCRS(int seed, int N, int cntInRow, crsMatrix& mtx)
{
	int i, j, k, f, tmp, notNull, c;

	srand(seed);

	notNull = cntInRow * N;
	AllocateCRSMatrix(N, notNull, mtx);

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < cntInRow; j++)
		{
			do
			{
				mtx.Col[i * cntInRow + j] = rand() % N;
				f = 0;
				for (k = 0; k < j; k++)
					if (mtx.Col[i * cntInRow + j] == mtx.Col[i * cntInRow + k])
						f = 1;
			} while (f == 1);
		}
		for (j = 0; j < cntInRow - 1; j++)
			for (k = 0; k < cntInRow - 1; k++)
				if (mtx.Col[i * cntInRow + k] > mtx.Col[i * cntInRow + k + 1])
				{
					tmp = mtx.Col[i * cntInRow + k];
					mtx.Col[i * cntInRow + k] = mtx.Col[i * cntInRow + k + 1];
					mtx.Col[i * cntInRow + k + 1] = tmp;
				}
	}

	for (i = 0; i < cntInRow * N; i++)
		mtx.Value[i] = next();//* MAX_VAL;

	c = 0;
	for (i = 0; i <= N; i++)
	{
		mtx.row_index[i] = c;
		c += cntInRow;
	}
	/*std::cout << std::endl;
	std::cout << "Gen Regular CRS \n";
	for (int j = 0; j < mtx.NZ; j++) {
		std::cout << mtx.Col[j] << "  ";
	}
	std::cout << std::endl;
	for (int j = 0; j < mtx.NZ; j++) {
		std::cout << mtx.Value[j] << "  ";
	}
	std::cout << std::endl;
	for (int j = 0; j < mtx.N+1; j++) {
		std::cout << mtx.row_index[j] << "  ";
	}
	std::cout << std::endl;*/

}
void GenerateSpecialCRS(int seed, int N, int cntInRow, crsMatrix& mtx)
{
	srand(seed);
	double end = pow((double)cntInRow, 1.0 / 3.0);
	double step = end / N;

	vector<int>* columns = new vector<int>[N];
	int NZ = 0;

	for (int i = 0; i < N; i++)
	{
		int rowNZ = int(pow((double(i + 1) * step), 3) + 1);
		NZ += rowNZ;
		int num1 = (rowNZ - 1) / 2;
		int num2 = rowNZ - 1 - num1;

		if (rowNZ != 0)
		{
			if (i < num1)
			{
				num2 += num1 - i;
				num1 = i;
				for (int j = 0; j < i; j++)
					columns[i].push_back(j);
				columns[i].push_back(i);
				for (int j = 0; j < num2; j++)
					columns[i].push_back(i + 1 + j);
			}
			else
			{
				if (N - i - 1 < num2)
				{
					num1 += num2 - (N - 1 - i);
					num2 = N - i - 1;
				}
				for (int j = 0; j < num1; j++)
					columns[i].push_back(i - num1 + j);
				columns[i].push_back(i);
				for (int j = 0; j < num2; j++)
					columns[i].push_back(i + j + 1);
			}
		}
	}

	AllocateCRSMatrix(N, NZ, mtx);

	int count = 0;
	int sum = 0;
	for (int i = 0; i < N; i++)
	{
		mtx.row_index[i] = sum;
		sum += columns[i].size();
		for (unsigned int j = 0; j < columns[i].size(); j++)
		{
			mtx.Col[count] = columns[i][j];
			mtx.Value[count] = next();
			count++;
		}
	}
	mtx.row_index[N] = sum;

	delete[] columns;
}
void GenerateVector(int seed, int N, Vec& vec) {
	srand(seed);
	for (int i = 0; i < vec.N; i++)
		vec.P[i][0] = next();
}

// Converting CRS matrix into regular one.
double CRStoReg(crsMatrix& iM, regMatrix& oM) {
	oM.N = iM.N;
	for (int i = 0; i < iM.N + 1; i++) {
		for (int j = iM.row_index[i]; j <= iM.row_index[i + 1] - 1; j++) {
			double a = iM.Value[j];
			int b = iM.Col[j];
			oM.P[i][b] = a;
		}
	}
	//std::cout << std::endl;
	//std::cout << "From CRS to Regular\n";
	/*for (int i = 0; i < iM.N; i++) {
		for (int j = 0; j < iM.N; j++) {
			std::cout << oM.P[i][j] << "  ";
		}
		std::cout << "\n";

	}*/
	return **oM.P;
}

// Matrix-Vector Multiplications.
// Sequencial approach.
double SeqRegMult(regMatrix A, Vec& x, Vec& b, double& time) {
	clock_t start = clock();

	for (int i = 0; i < A.N; i++)
	{
		double sum = 0;
		for (int j = 0; j < A.N; j++)
			sum = sum + (A.P[i][j] * x.P[j][0]);
		b.P[i][0] = sum;
	}
	/*std::cout << std::endl;
	std::cout << "Regular Mult \n";
	for (int j = 0; j < A.N; j++) {
		std::cout << b.P[j][0] << "  ";
	}
	std::cout << std::endl;*/
	clock_t finish = clock();
	time = (double)(finish - start) / CLOCKS_PER_SEC;

	return **b.P;
}
double SeqCRSMult(crsMatrix A, Vec& x, Vec& b, double& time) {
	clock_t start = clock();

	for (int i = 0; i < A.N; i++)
	{
		double sum = 0;
		int j1 = A.row_index[i];
		int j2 = A.row_index[i + 1];
		for (int j = j1; j < j2; j++)
			sum = sum + A.Value[j] * x.P[A.Col[j]][0];
		b.P[i][0] = sum;
	}
	/*std::cout << std::endl;
	std::cout << "CRS Mult\n";
	for (int j = 0; j < A.N; j++) {
		std::cout << b.P[j][0] << "  ";
	}
	std::cout << std::endl;*/
	clock_t finish = clock();
	time = (double)(finish - start) / CLOCKS_PER_SEC;

	return **b.P;
}

// Matrix-Vector Multiplications.
// Parallel approach.
double OmpRegMult(regMatrix A, Vec& x, Vec& b, int numThreads, double& time) {
	clock_t start = clock();
	omp_set_num_threads(numThreads);
	#pragma omp parallel for
	for (int i = 0; i < A.N; i++)
	{
		double sum = 0;
		for (int j = 0; j < A.N; j++)
			sum = sum + (A.P[i][j] * x.P[j][0]);
		b.P[i][0] = sum;
	}
	/*std::cout << std::endl;
	std::cout << "Regular Parallel Mult \n";
	for (int j = 0; j < A.N; j++) {
		std::cout << b.P[j][0] << "  ";
	}
	std::cout << std::endl;*/
	clock_t finish = clock();
	time = (double)(finish - start) /(CLOCKS_PER_SEC);

	return **b.P;
}
int OmpCRSMult(crsMatrix A, Vec& x, Vec& b, int numThreads, double& time) {
	clock_t start = clock();
	omp_set_num_threads(numThreads);
	#pragma omp parallel for
	for (int i = 0; i < A.N; i++)
	{
		double sum = 0;
		int j1 = A.row_index[i];
		int j2 = A.row_index[i + 1];
		for (int j = j1; j < j2; j++)
			sum = sum + A.Value[j] * x.P[A.Col[j]][0];
		b.P[i][0] = sum;
	}
	/*std::cout << std::endl;
	std::cout << "CRS parallel Mult\n";
	for (int j = 0; j < A.N; j++) {
		std::cout << b.P[j][0] << "  ";
	}
	std::cout << std::endl;*/
	clock_t finish = clock();
	time = (double)(finish - start) / CLOCKS_PER_SEC;

	return 0;
}

// Regular Matrix-Matrix Multiplications.
// Sequencial approach.
double SeqRegtoRegMult(regMatrix A, regMatrix B, regMatrix& C, double& time) {
	// Multiplying matrix A and B and storing in array C.
	for (int i = 0; i < A.N; ++i)
		for (int j = 0; j < B.N; ++j)
			for (int k = 0; k < A.N; ++k)
			{
				C.P[i][j] += A.P[i][k] * B.P[k][j];
			}

	// Displaying the multiplication of two matrices.
	std::cout << "\n" << "Output Matrix: \n";
	for (int i = 0; i < A.N; ++i)
		for (int j = 0; j < B.N; ++j)
		{
			std:: cout << " " << C.P[i][j];
			if (j == B.N - 1)
				std::cout << "\n";
		}
	return **C.P;
}

// Regular Matrix-Matrix Multiplications.
// Parallel approach.
double OmpRegtoRegMult(regMatrix A, regMatrix B, regMatrix& C, int numThreads, double& time) {
	// Multiplying matrix A and B and storing in array C.
	clock_t start = clock();
	omp_set_num_threads(numThreads);
    #pragma omp parallel for
	for (int i = 0; i < A.N; ++i)
		for (int j = 0; j < B.N; ++j)
			for (int k = 0; k < A.N; ++k)
			{
				C.P[i][j] += A.P[i][k] * B.P[k][j];
			}

	// Displaying the multiplication of two matrices.
	std::cout << "\n" << "Output Matrix: \n";
	for (int i = 0; i < A.N; ++i)
		for (int j = 0; j < B.N; ++j)
		{
			std::cout << " " << C.P[i][j];
			if (j == B.N - 1)
				std::cout << "\n";
		}
	clock_t finish = clock();
	time = (double)(finish - start) / (CLOCKS_PER_SEC);
	return **C.P;
}


double Transpose2(crsMatrix imtx, crsMatrix& omtx)
{
	clock_t start, finish;
	int i, j;

	start = clock();

	AllocateCRSMatrix(imtx.N, imtx.NZ, omtx);

	memset(omtx.row_index, 0, (imtx.N + 1) * sizeof(int));
	for (i = 0; i < imtx.NZ; i++)
		omtx.row_index[imtx.Col[i] + 1]++;

	int S = 0;
	for (i = 1; i <= imtx.N; i++)
	{
		int tmp = omtx.row_index[i];
		omtx.row_index[i] = S;
		S = S + tmp;
	}

	for (i = 0; i < imtx.N; i++)
	{
		int j1 = imtx.row_index[i];
		int j2 = imtx.row_index[i + 1];
		int Col = i; 
		for (j = j1; j < j2; j++)
		{
			double V = imtx.Value[j];  
			int RIndex = imtx.Col[j];  
			int IIndex = omtx.row_index[RIndex + 1];
			omtx.Value[IIndex] = V;
			omtx.Col[IIndex] = Col;
			omtx.row_index[RIndex + 1]++;
		}
	}

	finish = clock();

	return double(finish - start) / CLOCKS_PER_SEC;
}

// CRS Matrix-Matrix Multiplications.
 //Sequencial approach.
int SeqMult(crsMatrix A, crsMatrix B, crsMatrix& C, double& time)
{
	int ZERO_IN_CRS = 0;
	if (A.N != B.N)
		return 1;

	int N = A.N;
	vector<int> columns;
	vector<double> values;
	vector<int> row_index;

	clock_t start = clock();
	int rowNZ;
	row_index.push_back(0);
	for (int i = 0; i < N; i++)
	{
		rowNZ = 0;
		for (int j = 0; j < N; j++)
		{
			double sum = 0;
			for (int k = A.row_index[i]; k < A.row_index[i + 1]; k++)
			{
				for (int l = B.row_index[j]; l < B.row_index[j + 1]; l++)
				{
					if (A.Col[k] == B.Col[l])
					{
						sum += A.Value[k] * B.Value[l];
						break;
					}
				}
			}
			if (fabs(sum) > ZERO_IN_CRS)
			{
				columns.push_back(j);
				values.push_back(sum);
				rowNZ++;
			}
		}
		row_index.push_back(rowNZ + row_index[i]);
	}

	AllocateCRSMatrix(N, columns.size(), C);

	for (unsigned int j = 0; j < columns.size(); j++)
	{
		C.Col[j] = columns[j];
		C.Value[j] = values[j];
	}
	for (int i = 0; i <= N; i++)
		C.row_index[i] = row_index[i];

	clock_t finish = clock();

	time = (double)(finish - start) / CLOCKS_PER_SEC;

	return 0;
} 

// CRS Matrix-Matrix Multiplications.
// OpenMP approach.
int OmpMult(crsMatrix A, crsMatrix B, crsMatrix& C, int numThreads ,double& time)
{	
	int ZERO_IN_CRS = 0;
	if (A.N != B.N)
		return 1;
	omp_set_num_threads(numThreads);
	#pragma omp parallel for
	int N = A.N;
	vector<int> columns;
	vector<double> values;
	vector<int> row_index;

	clock_t start = clock();
	int rowNZ;
	row_index.push_back(0);
	for (int i = 0; i < N; i++)
	{
		rowNZ = 0;
		for (int j = 0; j < N; j++)
		{
			double sum = 0;
			for (int k = A.row_index[i]; k < A.row_index[i + 1]; k++)
			{
				for (int l = B.row_index[j]; l < B.row_index[j + 1]; l++)
				{
					if (A.Col[k] == B.Col[l])
					{
						sum += A.Value[k] * B.Value[l];
						break;
					}
				}
			}
			if (fabs(sum) > ZERO_IN_CRS)
			{
				columns.push_back(j);
				values.push_back(sum);
				rowNZ++;
			}
		}
		row_index.push_back(rowNZ + row_index[i]);
	}

	AllocateCRSMatrix(N, columns.size(), C);

	for (unsigned int j = 0; j < columns.size(); j++)
	{
		C.Col[j] = columns[j];
		C.Value[j] = values[j];
	}
	for (int i = 0; i <= N; i++)
		C.row_index[i] = row_index[i];

	clock_t finish = clock();

	time = (double)(finish - start) / CLOCKS_PER_SEC;

	return 0;
}



int main() {

	int N;
	int NZ;
	int cN;
	std::cout << "Enter matrix Size:\n";
	std::cin >> N;
	std::cout << "Enter number of non zero values:\n";
	std::cin >> NZ;
	std::cout << "Enter number of non zero in row:\n";
	std::cin >> cN;

	crsMatrix M = {};
	regMatrix reg = {};
	AllocateCRSMatrix(N, NZ, M);
	AllocateRegMatrix(N, reg);
	GenerateRegularCRS(5,N, cN, M);
	CRStoReg(M, reg);

	Vec x;
	AllocateVector(N, x);
	GenerateVector(5, N, x);

	Vec b1 = {};
	Vec b2 = {};
	AllocateVector(N, b1);
	AllocateVector(N, b2);
	double t1;
	double t2;
	// sequenctial Matrix-Vector Multiplications.
	SeqRegMult(reg, x, b1, t1);
	SeqCRSMult(M, x, b2, t2);

	double diff;
	CompareVectors(b1,b2,N, diff);
	if (diff > 0.001) {
		std::cout << "For n=100 correct" << std::endl;
	}
	else {
		std::cout << "For n=100 incorrect" << std::endl;
	}

	double t_reg_p;
	double t_crs_p;
	int numThreads = 8;
	// Parallel (using OpenMP) Matrix-Vector Multiplications.
	OmpRegMult(reg, x, b1, numThreads, t_reg_p);
	std::cout << "s_reg time " << t1 << std::endl;
	std::cout << "p_reg time " << t_reg_p << std::endl;
	std::cout << "time diff  " << t1 - t_reg_p << std::endl;
	
	OmpCRSMult(M, x, b2, numThreads ,t_crs_p);
	std::cout << "s_crs time " << t2 << std::endl;
	std::cout << "p_crs time " << t_reg_p << std::endl;
	std::cout << "time diff  " << t1 - t_reg_p << std::endl;

	// Sequencial CRS Matrix-Matrix Multiplications.
	crsMatrix cA, cB, cBT, cC;
	GenerateRegularCRS(1, N, NZ, cA);
	GenerateSpecialCRS(2, N, NZ, cB);
	double timeT = Transpose2(cB, cBT);
	double timeM1;
	SeqMult(cA, cBT, cC, timeM1);
	// Parallel CRS Matrix-Matrix Multiplications.
	crsMatrix cA_P, cB_P, cBT_P, cC_P;
	GenerateRegularCRS(1, N, NZ, cA_P);
	GenerateSpecialCRS(2, N, NZ, cB_P);
	double timeT_P = Transpose2(cB_P, cBT_P);
	double timeM1P;
	OmpMult(cA_P, cBT_P, cC_P ,numThreads,timeM1P);

	// Sequencial Regular Matrix-Matrix Multiplications.
	regMatrix rA = {};
	regMatrix rB = {};
	regMatrix rC;
	regMatrix rcC;
	CRStoReg(cA, rA);
	CRStoReg(cB, rB);
	double timeM2;
	SeqRegtoRegMult(rA, rB, rC, timeM2);
	// Parallel Regular Matrix-Matrix Multiplications.
	regMatrix rA_P = {};
	regMatrix rB_P = {};
	regMatrix rC_P;
	regMatrix rcC_P = {};
	CRStoReg(cA_P, rA_P);
	CRStoReg(cB_P, rB_P);
	double timeM2S;
	OmpRegtoRegMult(rA_P, rB_P, rC_P, numThreads ,timeM2S);
	
	// Comparing Sequencial results.
	double diff_mtxS;
	CRStoReg(cC, rcC);
	CompareMatrices(rC, rcC, diff_mtxS);
	if (diff_mtxS < 0.001)
		printf("OK\n");
	else
		printf("not OK\n");

	// Comparing Parallel results.
	double diff_mtxP;
	CRStoReg(cC_P, rcC_P);
	CompareMatrices(rC_P, rcC_P, diff_mtxP);
	if (diff_mtxP < 0.001)
		printf("OK\n");
	else
		printf("not OK\n");
	// parallel regular plotting
	/*double* T_reg;
	double* N_reg;
	for (int j = 2; j < 10; j = j + 2) {
		for (int i = 500; i < 5000; i = i + 500) {
			reg.N = i;
			x.N = i;
			b1.N = i;
			OmpRegMult(reg, x, b1, j, t_reg_p);
			T_reg[i] = t_reg_p;
		}
	}*/
}

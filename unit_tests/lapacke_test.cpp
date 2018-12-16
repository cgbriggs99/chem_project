#include <cppunit/TestCase.h>
#include <cppunit/TestFixture.h>
#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/XmlOutputter.h>

#include <lapacke.h>
#include <math.h>
using namespace CppUnit;
using namespace std;

class TestBasicMath: public CppUnit::TestFixture {
CPPUNIT_TEST_SUITE(TestBasicMath);
	CPPUNIT_TEST(testDGELS);
	CPPUNIT_TEST(testDGEEV);
	CPPUNIT_TEST_SUITE_END();

//public:
//    void setUp(void);
//    void tearDown(void);

protected:
	void testDGELS(void);
	void testDGEEV(void);

private:
//    CBasicMath *mTestObj;
};

void TestBasicMath::testDGELS(void) {

	double a[5][3] = { 1, 1, 1, 2, 3, 4, 3, 5, 2, 4, 2, 5, 5, 4, 3 };
	double b[5][2] = { -10, -3, 12, 14, 14, 12, 16, 16, 18, 16 };

	lapack_int info, m, n, lda, ldb, nrhs;
	int i, j;

	m = 5;
	n = 3;
	nrhs = 2;
	lda = 3;
	ldb = 2;

	info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, n, nrhs, *a, lda, *b, ldb);

	double expected[3][2] = {  2, 1, 1, 1, 1, 2 };

	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 2; c++) {

			double err = (expected[r][c] - b[r][c]) / expected[r][c];
			err = fabs(err);

			CPPUNIT_ASSERT( err < 1e-3);
		}
	}

}

#define N 5
#define LDA N
#define LDVL N
#define LDVR N

void TestBasicMath::testDGEEV(void) {


	// lapack_int LAPACKE_dgeev( int matrix_layout, char jobvl, char jobvr,
    //lapack_int n, double* a, lapack_int lda, double* wr,
    //double* wi, double* vl, lapack_int ldvl, double* vr,
    //lapack_int ldvr );

	lapack_int n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info;

	/* Local arrays */
	double wr[N], wi[N], vl[LDVL*N], vr[LDVR*N];
	double a[LDA*N] = {
	   -1.01,  0.86, -4.60,  3.31, -4.81,
		3.98,  0.53, -7.04,  5.29,  3.55,
		3.30,  8.26, -3.89,  8.20, -1.51,
		4.43,  4.96, -7.66, -7.33,  6.18,
		7.31, -6.43, -6.16,  2.47,  5.58
	};

	/* Executable statements */
	printf( "LAPACKE_dgeev (row-major, high-level) Example Program Results\n" );

	/* Solve eigenproblem */
	info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'V', 'V', n, a, lda, wr, wi,
					vl, ldvl, vr, ldvr );

	double expected_values[9] = {
			2.86, 10.76,   2.86,-10.76,  -0.69,  4.70,  -0.69, -4.70 -10.46
	 };

	double expected_left[5][9] = {
	   0.04,  0.29,   0.04, -0.29,  -0.13, -0.33,  -0.13,  0.33,   0.04,
	   0.62,  0.00,   0.62,  0.00,   0.69,  0.00,   0.69,  0.00,   0.56,
	  -0.04, -0.58,  -0.04,  0.58,  -0.39, -0.07,  -0.39,  0.07,  -0.13,
	   0.28,  0.01,   0.28, -0.01,  -0.02, -0.19,  -0.02,  0.19,  -0.80,
	  -0.04,  0.34,  -0.04, -0.34,  -0.40,  0.22,  -0.40, -0.22,   0.18
	  };

	double expected_right[5][9] = {
	   0.11,  0.17,   0.11, -0.17,   0.73,  0.00,   0.73,  0.00,   0.46,
	   0.41, -0.26,   0.41,  0.26,  -0.03, -0.02,  -0.03,  0.02,   0.34,
	   0.10, -0.51,   0.10,  0.51,   0.19, -0.29,   0.19,  0.29,   0.31,
	   0.40, -0.09,   0.40,  0.09,  -0.08, -0.08,  -0.08,  0.08,  -0.74,
	   0.54,  0.00,   0.54,  0.00,  -0.29, -0.49,  -0.29,  0.49,   0.16
	};


//	/* Auxiliary routine: printing eigenvectors */
//	void print_eigenvectors( char* desc, MKL_INT n, double* wi, double* v, MKL_INT ldv ) {
//	        MKL_INT i, j;
//	        printf( "\n %s\n", desc );
	   for(int i = 0; i < N; i++ ) {
	      int j = 0;
	      while( j < N ) {
	         if( wi[j] == (double)0.0 ) {
	            printf( " %6.2f", vl[i*ldvl+j] );
	            j++;
	         } else {
	            printf( " (%6.2f,%6.2f)", vl[i*ldvl+j], vl[i*ldvl+(j+1)] );
	            printf( " (%6.2f,%6.2f)", vl[i*ldvl+j], -vl[i*ldvl+(j+1)] );
	            j += 2;
	         }
	      }
	      printf( "\n" );
	   }
}


//void TestBasicMath::setUp(void)
//{
//    mTestObj = new CBasicMath();
//}

//void TestBasicMath::tearDown(void)
//{
//    delete mTestObj;
//}

CPPUNIT_TEST_SUITE_REGISTRATION(TestBasicMath);

int main(int argc, char* argv[]) {
	// informs test-listener about testresults
	CPPUNIT_NS::TestResult testresult;
	// register listener for collecting the test-results
	CPPUNIT_NS::TestResultCollector collectedresults;
	testresult.addListener(&collectedresults);
	// register listener for per-test progress output
	CPPUNIT_NS::BriefTestProgressListener progress;
	testresult.addListener(&progress);

	// insert test-suite at test-runner by registry
	CPPUNIT_NS::TestRunner testrunner;
	testrunner.addTest(
			CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest());
	testrunner.run(testresult);

	// output results in compiler-format
	CPPUNIT_NS::CompilerOutputter compileroutputter(&collectedresults,
			std::cerr);
	compileroutputter.write();
	// Output XML for Jenkins CPPunit plugin
	ofstream xmlFileOut("cppTestBasicMathResults.xml");
	XmlOutputter xmlOut(&collectedresults, xmlFileOut);
	xmlOut.write();

	// return 0 if tests were successful
	return collectedresults.wasSuccessful() ? 0 : 1;
}

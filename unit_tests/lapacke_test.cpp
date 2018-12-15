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
	CPPUNIT_TEST(testDGELS);CPPUNIT_TEST_SUITE_END()
	;

//public:
//    void setUp(void);
//    void tearDown(void);

protected:
	void testDGELS(void);

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

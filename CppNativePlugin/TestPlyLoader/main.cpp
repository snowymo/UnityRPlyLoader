#include "../PlyLoader/PlyLoader.h"

int main() {

	PlyFileObject* test = LoadPly("C:\\Scans\\currentScan.ply");

	unsigned int vCnt = 0;
	float* verts = new float[vCnt * 3];
	verts = GetPlyVerts(test, vCnt);

	unsigned int cCnt = 0;
	unsigned char* colors = new unsigned char[cCnt * 3];
	colors = GetPlyColors(test, cCnt);

	unsigned int iCnt = 0;
	unsigned int* faces = new unsigned int[iCnt * 3];
	faces = GetPlyIndexs(test, iCnt);

	delete[] verts;
	delete[] colors;
	delete[] faces;

	return 0;
}
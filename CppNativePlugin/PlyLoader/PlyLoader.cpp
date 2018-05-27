#include "PlyLoader.h"

#define   STB_IMAGE_IMPLEMENTATION

#include <Eigen/Core>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <map>
#include <memory>
#include "stb_image.h"
extern "C" {
#include "rply.h"
}



#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <vcg/math/base.h>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/complex/algorithms/attribute_seam.h>

using namespace vcg;

// global var and funcs

const double delta = 6.0 / 29.0;
const double delta3 = delta * delta * delta;

//reference white
double Xn = 62289;
double Yn = 65535;
double Zn = 71356;

double f(double t)
{
	if (t > delta3)
		return pow(t, 1.0 / 3.0);
	else
		return t / (3 * delta * delta) + 4.0 / 29.0;
}


Eigen::Matrix<unsigned short, 3, 1> RGBToLab(const Eigen::Matrix<unsigned short, 3, 1>& rgb)
{
	Eigen::Matrix3d mRGBToXYZ;
	mRGBToXYZ <<
		0.4124564, 0.3575761, 0.1804375,
		0.2126729, 0.7151522, 0.0721750,
		0.0193339, 0.1191920, 0.9503041;

	auto rgbDbl = rgb.cast<double>();
	auto XYZ = mRGBToXYZ * rgbDbl;

	double fx = f(XYZ.x() / Xn);
	double fy = f(XYZ.y() / Yn);
	double fz = f(XYZ.z() / Zn);

	double L = 116 * fy - 16;
	double a = 500 * (fx - fy);
	double b = 200 * (fy - fz);

	L = std::max(0.0, std::min(65535.0, 655.35 * L));
	a = std::max(0.0, std::min(65535.0, 327.67 * a + 32768));
	b = std::max(0.0, std::min(65535.0, 327.67 * b + 32768));

	return Eigen::Matrix<unsigned short, 3, 1>((unsigned short)L, (unsigned short)a, (unsigned short)b);
}


inline std::vector<std::string> &str_tokenize(const std::string &s, char delim, std::vector<std::string> &elems, bool include_empty = false) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
		if (!item.empty() || include_empty)
			elems.push_back(item);
	return elems;
}

inline std::vector<std::string> str_tokenize(const std::string &s, char delim, bool include_empty) {
	std::vector<std::string> elems;
	str_tokenize(s, delim, elems, include_empty);
	return elems;
}
inline uint32_t str_to_uint32_t(const std::string &str) {
	char *end_ptr = nullptr;
	uint32_t result = (uint32_t)strtoul(str.c_str(), &end_ptr, 10);
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse unsigned integer \"" + str + "\"");
	return result;
}

inline bool file_exists(const std::string& path)
{
	//implementation from https://stackoverflow.com/a/12774387/1210053
	if (FILE *file = fopen(path.c_str(), "r"))
	{
		fclose(file);
		return true;
	}
	else
	{
		return false;
	}
}

// end of global var and funcs

/////// FUNCTIONS NEEDED BY "UV WEDGE TO VERTEX" FILTER
inline void ExtractVertex(const CMeshO & srcMesh, const CMeshO::FaceType & f, int whichWedge, const CMeshO & dstMesh, CMeshO::VertexType & v)
{
	(void)srcMesh;
	(void)dstMesh;
	// This is done to preserve every single perVertex property
	// perVextex Texture Coordinate is instead obtained from perWedge one.
	v.ImportData(*f.cV(whichWedge));
	v.T() = f.cWT(whichWedge);
}

inline bool CompareVertex(const CMeshO & m, const CMeshO::VertexType & vA, const CMeshO::VertexType & vB)
{
	(void)m;
	return (vA.cT() == vB.cT());
}

PlyFileObject::PlyFileObject(const char* fileName, int downsample)
{
	load_ply(fileName, F, V, N, C, true);
}

void PlyFileObject::Enable(int openingFileMask)
{
	if (openingFileMask & tri::io::Mask::IOM_VERTTEXCOORD)
		updateDataMask(MM_VERTTEXCOORD);
	if (openingFileMask & tri::io::Mask::IOM_WEDGTEXCOORD)
		updateDataMask(MM_WEDGTEXCOORD);
	if (openingFileMask & tri::io::Mask::IOM_VERTCOLOR)
		updateDataMask(MM_VERTCOLOR);
	if (openingFileMask & tri::io::Mask::IOM_FACECOLOR)
		updateDataMask(MM_FACECOLOR);
	if (openingFileMask & tri::io::Mask::IOM_VERTRADIUS) updateDataMask(MM_VERTRADIUS);
	if (openingFileMask & tri::io::Mask::IOM_CAMERA) updateDataMask(MM_CAMERA);
	if (openingFileMask & tri::io::Mask::IOM_VERTQUALITY) updateDataMask(MM_VERTQUALITY);
	if (openingFileMask & tri::io::Mask::IOM_FACEQUALITY) updateDataMask(MM_FACEQUALITY);
	if (openingFileMask & tri::io::Mask::IOM_BITPOLYGONAL) updateDataMask(MM_POLYGONAL);
}

void PlyFileObject::updateDataMask(int neededDataMask)
{
	if ((neededDataMask & MM_FACEFACETOPO) != 0)
	{
		cm.face.EnableFFAdjacency();
		tri::UpdateTopology<CMeshO>::FaceFace(cm);
	}
	if ((neededDataMask & MM_VERTFACETOPO) != 0)
	{
		cm.vert.EnableVFAdjacency();
		cm.face.EnableVFAdjacency();
		tri::UpdateTopology<CMeshO>::VertexFace(cm);
	}

	if ((neededDataMask & MM_WEDGTEXCOORD) != 0)
		cm.face.EnableWedgeTexCoord();
	if ((neededDataMask & MM_FACECOLOR) != 0)
		cm.face.EnableColor();
	if ((neededDataMask & MM_FACEQUALITY) != 0)
		cm.face.EnableQuality();
	if ((neededDataMask & MM_FACECURVDIR) != 0)
		cm.face.EnableCurvatureDir();
	if ((neededDataMask & MM_FACEMARK) != 0)
		cm.face.EnableMark();
	if ((neededDataMask & MM_VERTMARK) != 0)
		cm.vert.EnableMark();
	if ((neededDataMask & MM_VERTCURV) != 0)
		cm.vert.EnableCurvature();
	if ((neededDataMask & MM_VERTCURVDIR) != 0)
		cm.vert.EnableCurvatureDir();
	if ((neededDataMask & MM_VERTRADIUS) != 0)
		cm.vert.EnableRadius();
	if ((neededDataMask & MM_VERTTEXCOORD) != 0)
		cm.vert.EnableTexCoord();
}

void PlyFileObject::load_ply(const std::string &filename, MatrixXu &F, Matrix3Xf &V, Matrix3Xf &N, Matrix3Xuc &C, bool pointcloud)
{
	auto message_cb = [](p_ply ply, const char *msg) { std::cerr << "rply: " << msg << std::endl; };

	std::cout << "Loading \"" << filename << "\" .. ";
	std::cout.flush();

	p_ply ply = ply_open(filename.c_str(), message_cb, 0, nullptr);
	if (!ply)
		throw std::runtime_error("Unable to open PLY file \"" + filename + "\"!");

	if (!ply_read_header(ply)) {
		ply_close(ply);
		throw std::runtime_error("Unable to open PLY header of \"" + filename + "\"!");
	}

	const float gamma = 2.2f;

	p_ply_element element = nullptr;
	uint32_t vertexCount = 0, faceCount = 0;

	p_ply_property prop = nullptr;
	const char* pname;
	e_ply_type type, length, value;

	bool hasUv = false;
	bool hasNormals = false;
	bool hasColor = false;

	/* Inspect the structure of the PLY file, load number of faces if avaliable */
	while ((element = ply_get_next_element(ply, element)) != nullptr) {
		const char *name;
		long nInstances;

		ply_get_element_info(element, &name, &nInstances);
		if (!strcmp(name, "vertex"))
		{
			vertexCount = (uint32_t)nInstances;
			while ((prop = ply_get_next_property(element, prop)) != nullptr)
			{
				ply_get_property_info(prop, &pname, &type, &length, &value);
				if (!strcmp(pname, "texture_u"))
					hasUv = true;
				else if (!strcmp(pname, "nx"))
					hasNormals = true;
				else if (!strcmp(pname, "red"))
					hasColor = true;
			}
		}
		else if (!strcmp(name, "face"))
			faceCount = (uint32_t)nInstances;
	}

	if (vertexCount == 0 && faceCount == 0)
		throw std::runtime_error("PLY file \"" + filename + "\" is invalid! No face/vertex/elements found!");
	else if (faceCount == 0)
		pointcloud = true;

	F.resize(3, faceCount);
	V.resize(3, vertexCount);
	N.resize(3, hasNormals ? vertexCount : 0);
	C.resize(3, hasColor ? vertexCount : 0);
	Eigen::Matrix<float, 2, Eigen::Dynamic> uv(2, hasUv ? vertexCount : 0);


	struct VertexCallbackData {
		Matrix3Xf &V;
		VertexCallbackData(Matrix3Xf &V)
			: V(V) { }
	};

	struct FaceCallbackData {
		MatrixXu &F;
		FaceCallbackData(MatrixXu &F)
			: F(F) { }
	};

	struct VertexNormalCallbackData {
		Matrix3Xf &N;
		VertexNormalCallbackData(Matrix3Xf &_N)
			: N(_N) { }
	};

	struct VertexUVCallbackData {
		Eigen::Matrix<float, 2, Eigen::Dynamic> &uv;
		VertexUVCallbackData(Eigen::Matrix<float, 2, Eigen::Dynamic> &_uv)
			: uv(_uv) { }
	};

	struct VertexColorCallbackData {
		Eigen::Matrix<unsigned char, 3, Eigen::Dynamic> &c;
		VertexColorCallbackData(Eigen::Matrix<unsigned char, 3, Eigen::Dynamic> &c)
			: c(c) { }
	};

	auto rply_vertex_cb = [](p_ply_argument argument) -> int {
		VertexCallbackData *data; long index, coord;
		ply_get_argument_user_data(argument, (void **)&data, &coord);
		ply_get_argument_element(argument, nullptr, &index);
		data->V(coord, index) = (Float)ply_get_argument_value(argument);
		return 1;
	};

	auto rply_vertex_normal_cb = [](p_ply_argument argument) -> int {
		VertexNormalCallbackData *data; long index, coord;
		ply_get_argument_user_data(argument, (void **)&data, &coord);
		ply_get_argument_element(argument, nullptr, &index);
		data->N(coord, index) = (Float)ply_get_argument_value(argument);
		return 1;
	};

	auto rply_vertex_uv_cb = [](p_ply_argument argument) -> int
	{
		VertexUVCallbackData *data; long index, coord;
		ply_get_argument_user_data(argument, (void **)&data, &coord);
		ply_get_argument_element(argument, nullptr, &index);
		data->uv(coord, index) = (Float)ply_get_argument_value(argument);
		return 1;
	};

	auto rply_vertex_color_cb = [](p_ply_argument argument) -> int
	{
		VertexColorCallbackData *data; long index, coord;
		ply_get_argument_user_data(argument, (void **)&data, &coord);
		ply_get_argument_element(argument, nullptr, &index);
		auto colorInPly = ply_get_argument_value(argument);
		data->c(coord, index) = (unsigned short)(std::pow(colorInPly / 255.0, 2.2) * 65535);
		return 1;
	};

	auto rply_index_cb = [](p_ply_argument argument) -> int {
		FaceCallbackData *data;
		long length, value_index, index;
		ply_get_argument_property(argument, nullptr, &length, &value_index);

		if (length != 3)
			throw std::runtime_error("Only triangle faces are supported!");

		ply_get_argument_user_data(argument, (void **)&data, nullptr);
		ply_get_argument_element(argument, nullptr, &index);

		if (value_index >= 0)
			data->F(value_index, index) = (uint32_t)ply_get_argument_value(argument);

		return 1;
	};

	VertexCallbackData vcbData(V);
	FaceCallbackData fcbData(F);
	VertexNormalCallbackData vncbData(N);
	VertexUVCallbackData vuvcbData(uv);
	VertexColorCallbackData vcData(C);

	ply_set_read_cb(ply, "vertex", "x", rply_vertex_cb, &vcbData, 0);
	ply_set_read_cb(ply, "vertex", "y", rply_vertex_cb, &vcbData, 1);
	ply_set_read_cb(ply, "vertex", "z", rply_vertex_cb, &vcbData, 2);

	ply_set_read_cb(ply, "vertex", "nx", rply_vertex_normal_cb, &vncbData, 0);
	ply_set_read_cb(ply, "vertex", "ny", rply_vertex_normal_cb, &vncbData, 1);
	ply_set_read_cb(ply, "vertex", "nz", rply_vertex_normal_cb, &vncbData, 2);

	if (faceCount > 0)
	{
		if (!ply_set_read_cb(ply, "face", "vertex_index", rply_index_cb, &fcbData, 0) && !ply_set_read_cb(ply, "face", "vertex_indices", rply_index_cb, &fcbData, 0))
		{
			ply_close(ply);
			throw std::runtime_error("PLY file \"" + filename + "\" does not contain vertex indices!");
		}
	}

	ply_set_read_cb(ply, "vertex", "texture_u", rply_vertex_uv_cb, &vuvcbData, 0);
	ply_set_read_cb(ply, "vertex", "texture_v", rply_vertex_uv_cb, &vuvcbData, 1);

	ply_set_read_cb(ply, "vertex", "red", rply_vertex_color_cb, &vcData, 0);
	ply_set_read_cb(ply, "vertex", "green", rply_vertex_color_cb, &vcData, 1);
	ply_set_read_cb(ply, "vertex", "blue", rply_vertex_color_cb, &vcData, 2);

	if (!ply_read(ply)) {
		ply_close(ply);
		throw std::runtime_error("Error while loading PLY data from \"" + filename + "\"!");
	}

	ply_close(ply);

	if (hasUv && !hasColor)
	{
		std::string texturePath = filename;
		//change extension from ply to png
		texturePath[texturePath.length() - 2] = 'n';
		texturePath[texturePath.length() - 1] = 'g';

		if (file_exists(texturePath))
		{
			std::cout << "loading texture .. ";
			C.resizeLike(V);
			using handleType = std::unique_ptr<uint8_t[], void(*)(void*)>;
			int w, h, n;
			handleType textureData(stbi_load(texturePath.c_str(), &w, &h, &n, 3), stbi_image_free);
			if (textureData)
			{
#pragma omp parallel for
				for (int i = 0; i < V.cols(); ++i)
				{
					auto vuv = uv.col(i);
					int x = (int)std::floor(w * vuv.x());
					int y = (int)std::floor(h * (1 - vuv.y()));

					if (x >= w)
						x = w - 1;
					if (y >= h)
						y = h - 1;

					int offset = n * (x + w * y);

					//Gamma de-correction

					C.col(i) <<
						(unsigned char)(std::pow((float)textureData.get()[offset + 0] / 255.0f, gamma) * 255.0f),
						(unsigned char)(std::pow((float)textureData.get()[offset + 1] / 255.0f, gamma) * 255.0f),
						(unsigned char)(std::pow((float)textureData.get()[offset + 2] / 255.0f, gamma) * 255.0f);
// 					if (i % 100 == 0) {
// 						std::cout << "i:" << i << "\tcolor\t" << C.col(i) << "\t" << (int)C(0,i) << "\n";
// 							}
					 			
				}
			}
		}

	}

// #pragma omp parallel for
// 	for (int i = 0; i < C.cols(); ++i)
// 	{
// 		// i prefer rgb than lab
// 		/*C.col(i) = RGBToLab(C.col(i));*/
// 		C.col(i) /= 256;
// 		if (i % 100 == 0)
// 			std::cout << "i:" << i << "\tcolor\t" << C.col(i) << "\n";
// 	}

	std::cout << "done. (V=" << vertexCount;
	if (faceCount > 0)
		std::cout << ", F=" << faceCount;
}

__declspec(dllexport) PlyFileObject* LoadPly(const char* fileName)
{
	return new PlyFileObject(fileName);
}

PlyFileObject* LoadPlyDownSample(const char* fileName, int downsample)
{
	return new PlyFileObject(fileName, downsample);
}

__declspec(dllexport) void UnLoadPly(PlyFileObject* plyFile)
{
	delete plyFile;
}

__declspec(dllexport) float* GetPlyVerts(PlyFileObject* plyFile, unsigned int& count)
{
	if (plyFile == NULL)
		return NULL;

	count = (unsigned int)plyFile->V.cols();
	return (plyFile->V).data();
}

__declspec(dllexport) float* GetPlyNormals(PlyFileObject* plyFile, unsigned int& count)
{
	if (plyFile == NULL)
		return NULL;

	count = (unsigned int)plyFile->N.cols();
	return plyFile->N.data();
}

__declspec(dllexport) unsigned char* GetPlyColors(PlyFileObject* plyFile, unsigned int& count)
{
	if (plyFile == NULL)
		return NULL;

	count = (unsigned int)(plyFile->C).cols();
	return (plyFile->C).data();
}

__declspec(dllexport) unsigned int* GetPlyIndexs(PlyFileObject* plyFile, unsigned int& count)
{
	if (plyFile == NULL)
		return NULL;

	count = (unsigned int)plyFile->F.cols();
	return (plyFile->F).data();
}

__declspec(dllexport) float* GetPlyUvs(PlyFileObject* plyFile, unsigned int& count)
{
	if (plyFile == NULL)
		return NULL;

	count = (unsigned int)plyFile->uvCoords.size();
	return plyFile->uvCoords.data();
}

__declspec(dllexport) const char* GetPlyTextureName(PlyFileObject* plyFile)
{
	if (plyFile == NULL || plyFile->cm.textures.size() == 0)
		return NULL;

	return plyFile->textureName.c_str();
}
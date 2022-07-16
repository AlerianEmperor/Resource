#ifndef _SCENE_H_
#define _SCENE_H_
#define _CRT_SECURE_NO_WARNINGS
//#include "GGX_Material.h"
#include "AreaLight.h"
//#include "EnvLight.h"
//#include "Triangle.h"
//#include "Triangle_Save_Space.h"
#include "Triangle_16.h"
#include "BBox.h"
#include "material.h"
#include "BSDF.h"
//#include "RoughPlastic.h"
#include "Enviroment.h"
//#include "All_Material.h"

#include "onb.h"
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <fstream>
#include <time.h>
//#include <omp.h>

#define max_depth 5
#define num_bin 12
#define max_bin 60
#define imax_bin 1 / max_bin
#define ct 1
#define ci 4

#define Light_Path_Length 2
#define Eye_Path_Length 3

using namespace std;

#define pi 3.1415926f
#define tau 6.2831853f
#define ipi 1.0f / pi
#define itau 1.0f / tau

#define i65535 1.0f / 65535

void fixIndex(int& v, const int& n)
{
	v = v < 0 ? v + n : v > 0 ? v - 1 : -1;
}

void fixIndex(int v[3], const int& n)
{
	v[0] = v[0] < 0 ? v[0] + n : v[0] > 0 ? v[0] - 1 : -1;
	v[1] = v[1] < 0 ? v[1] + n : v[1] > 0 ? v[1] - 1 : -1;
	v[2] = v[2] < 0 ? v[2] + n : v[2] > 0 ? v[2] - 1 : -1;

}

void SkipSpace(char *&t)
{
	t += strspn(t, " \t");
}

void getdirection(const string& filepath, string& direction)
{
	size_t found = filepath.find_last_of('/\\');
	direction = filepath.substr(0, found + 1);
}

void getfilename(const string& filepath, string& filename)
{

	size_t found = filepath.find_last_of('/\\');
	filename = filepath.substr(found + 1);
}

void getdDirectionAndName(const string& filepath, string& direction, string& filename)
{
	size_t found = filepath.find_last_of('/\\');
	direction = filepath.substr(0, found + 1);
	filename = filepath.substr(found + 1);
}


struct SimpleTriangle
{
	vec3 c;
	int i;
	BBox b;
	SimpleTriangle() {}
	SimpleTriangle(vec3 c_, int i_) : c(c_), i(i_) {}
};

struct node
{
	BBox box;

	node *nodes[2];// = { NULL };

	unsigned start : 25;
	unsigned range : 4;
	unsigned axis : 2;
	unsigned leaf : 1;


	node() { }
	node(int s, int r) : start(s), range(r) {}
};

struct Normal_Struct
{
	vec3 normal = vec3(0.0f, 0.0f, 0.0f);
	int count = 0;
};


void delete_tree(node* n)
{
	if (n->leaf)
	{
		free(n);
		return;
	}
	delete_tree(n->nodes[0]);
	delete_tree(n->nodes[1]);

	free(n);
}

void delete_bvh(node** n)
{
	delete_tree(*n);
	*n = NULL;
}

enum face
{
	Single_Line, Double_Line
};

void get_face_index(char*& t, const int& vs, const int& vts, const int& vns, vector<Face>& faces)
{
	string s = t;
	int length = s.find_last_of("0123456789");
	s = s.substr(0, length + 1);

	int sign = 1;
	int count = 0;
	vector<int> index;
	//vector<int> indices;
	int face_type = Single_Line;

	int num_data_per_vertex = 0;
	bool found_num_data_per_vertex = false;

	for (int i = 0; i <= length + 1; ++i)
	{
		if (s[i] == '-')
		{
			sign = -1;
		}
		else if (isdigit(s[i]))
		{
			count = 10 * count + s[i] - '0';
		}
		else if (s[i] == '/')
		{
			if (!found_num_data_per_vertex)
			{
				++num_data_per_vertex;
			}
			face_type = Single_Line;
			index.emplace_back(sign * count);
			sign = 1;
			count = 0;

			if (s[i + 1] == '/')
			{
				face_type = Double_Line;
				++i;
			}
		}
		else if (s[i] == ' ')
		{

			index.emplace_back(sign * count);
			sign = 1;
			count = 0;
			if (!found_num_data_per_vertex)
			{
				++num_data_per_vertex;
				found_num_data_per_vertex = true;
			}
		}
		else if (i == length + 1)
		{
			index.emplace_back(sign * count);
			sign = 1;
			break;
		}
	}

	//cout << index.size() << "\n";
	//cout << "\n";
	int size = index.size();

	//cout << num_data_per_vertex << "\n";
	//if (size % 2 == 0 && size % 3 == 0) //6 12 18 24 30
	//{
	if (num_data_per_vertex == 3) // v/vt/vn case   12 18 24 30 
	{
		//cout << "line 3\n";
		for (int i = 0; i < size; i += 3)
		{
			fixIndex(index[i], vs);
			fixIndex(index[i + 1], vts);
			fixIndex(index[i + 2], vns);
		}

		int start_v = 0;
		int start_vt = 1;
		int start_vn = 2;

		int num_Triangle = size / 3 - 2;

		for (int i = 0; i < num_Triangle; ++i)
		{
			Face f(index[start_v], index[start_vt], index[start_vn], index[3 * i + 3], index[3 * i + 4], index[3 * i + 5], index[3 * i + 6], index[3 * i + 7], index[3 * i + 8]);
			faces.emplace_back(f);
		}
	}
	else if (num_data_per_vertex == 2)  //  v/vt or v//vn
	{
		//cout << "line 2\n";
		if (face_type == Single_Line) // v / vt
		{
			//cout << "single\n";
			for (int i = 0; i < size; i += 2)
			{
				fixIndex(index[i], vs);
				fixIndex(index[i + 1], vts);
			}

			int num_Triangle = size / 2 - 2;

			int start_v = 0;
			int start_vt = 1;

			for (int i = 0; i < num_Triangle; ++i)
			{
				//01 23 45    01 45 67
				Face f(index[start_v], index[start_vt], -1, index[2 * i + 2], index[2 * i + 3], -1, index[2 * i + 4], index[2 * i + 5], -1);
				faces.emplace_back(f);
			}
		}
		else if (face_type == Double_Line)// v // vn
		{
			//cout << "double\n";
			for (int i = 0; i < size; i += 2)
			{
				fixIndex(index[i], vs);
				fixIndex(index[i + 1], vns);
			}

			int num_Triangle = size / 2 - 2;

			int start_v = 0;
			int start_vn = 1;

			for (int i = 0; i < num_Triangle; ++i)
			{
				Face f(index[start_v], -1, index[start_vn], index[2 * i + 2], -1, index[2 * i + 3], index[2 * i + 4], -1, index[2 * i + 5]);
				faces.emplace_back(f);
			}
		}
	}
}



static float RGB_To_GrayScale(const float& R, const float& G, const float& B)
{
	return 0.299f * R + 0.587f * G + 0.114f * B;
}

/*static int clamp(const int& i, const int& low, const int& high)
{
return minf(maxf(low, i), high);
}*/

static int get_pixel_index(const int& i, const int& j, const int& w, const int& h)
{
	int u = clamp(i, 0, w - 1);
	int v = clamp(j, 0, h - 1);

	return v * w + u;
}

static vector<vec3> Bump_Map_To_Normal_Map(const int& w, const int& h, const int& n, const float& bump_strength, const vector<float>& grayscale)
{
	vector<vec3> normal_map;

	normal_map.resize(w * h);

	for (int i = 0; i < w; ++i)
	{
		for (int j = 0; j < h; ++j)
		{
			int top_left = get_pixel_index(i - 1, j - 1, w, h);
			int top = get_pixel_index(i - 1, j, w, h);
			int top_right = get_pixel_index(i - 1, j + 1, w, h);

			int left = get_pixel_index(i, j - 1, w, h);
			int center = get_pixel_index(i, j, w, h);
			int right = get_pixel_index(i, j + 1, w, h);


			int bottom_left = get_pixel_index(i + 1, j - 1, w, h);
			int bottom = get_pixel_index(i + 1, j, w, h);
			int bottom_right = get_pixel_index(i + 1, j + 1, w, h);

			const float tl = grayscale[top_left];
			const float t = grayscale[top];
			const float tr = grayscale[top_right];

			const float l = grayscale[left];
			const float c = grayscale[center];
			const float r = grayscale[right];

			const float bl = grayscale[bottom_left];
			const float b = grayscale[bottom];
			const float br = grayscale[bottom_right];

			// sobel filter
			const float dx = (tr + 2.0 * r + br) - (tl + 2.0 * l + bl);
			const float dy = (bl + 2.0 * b + br) - (tl + 2.0 * t + tr);
			const float dz = 1.0f / (bump_strength);

			vec3 n(dx, dy, dz);

			/*
			n.normalize();

			n = vec3((n.x + 1.0f) * 0.5f, (n.y + 1.0f) * 0.5f, (n.z + 1.0f) * 0.5f);

			normal_map[j * w + i] = n.norm();

			*/

			n.normalize();

			n = vec3((n.x + 1.0f) * 0.5f, (n.y + 1.0f) * 0.5f, (n.z + 1.0f) * 0.5f);

			normal_map[j * w + i] = n;

		}
	}
	return normal_map;
}

/*long long int hash_function( int& a,  int& b,  int& c)
{
	return ((a << 16) | (b << 8) | (c));
}
*/

/*long long int hash_function(const char* s, unsigned salt)
{
	unsigned h = salt;
	while (*s)
		h = h * 101 + (unsigned)*s++;
	return h;
}*/

//https://stackoverflow.com/questions/1358468/how-to-create-unique-integer-number-from-3-different-integers-numbers1-oracle-l/1362712
//Hash Function Formula
//h = (a * P1 + b) * P2 + c 
//P1 and P2 is very large Prime

//maximum posible for SanMiguel is 6721532 vertex normal size
//so I will choose P1 and P2 larger than this maxmimum
//consider chossing bigger prime if number of vertex grow larger than one of these two
//Prime taken from here
//http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php
//last two primes in the bottom right conner
//6727583
//6727741

/*
const int P1 = 6727583;
const int P2 = 6727741;

long long int prime_hash_function(const int& a, const int& b, const int& c)
{
	//return (a * P1 + b) * P2 + c;
	return a * P1 + b;
}
*/

//murmur hash
/*
static inline uint32_t murmur_32_scramble(uint32_t k) {
	k *= 0xcc9e2d51;
	k = (k << 15) | (k >> 17);
	k *= 0x1b873593;
	return k;
}
uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed)
{
	uint32_t h = seed;
	uint32_t k;
	//Read in groups of 4. 
	for (size_t i = len >> 2; i; i--) {
		// Here is a source of differing results across endiannesses.
		// A swap here has no effects on hash properties though.
		memcpy(&k, key, sizeof(uint32_t));
		key += sizeof(uint32_t);
		h ^= murmur_32_scramble(k);
		h = (h << 13) | (h >> 19);
		h = h * 5 + 0xe6546b64;
	}
	/Read the rest.
	k = 0;
	for (size_t i = len & 3; i; i--) {
		k <<= 8;
		k |= key[i - 1];
	}
	// A swap is *not* necessary here because the preceding loop already
	// places the low bytes in the low places according to whatever endianness
	// we use. Swaps only apply when the memory is copied in a chunk.
	h ^= murmur_32_scramble(k);
	//Finalize
	h ^= len;
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;
	return h;
}
*/
/*
static long long int compute_hash(const string& s) {
	const int p = 6727583;
	const int m = 1e9 + 9;
	long long hash_value = 0;
	long long p_pow = 1;
	for (char c : s) {
		hash_value = (hash_value + (c - '0' + 1) * p_pow) % m;
		p_pow = (p_pow * p) % m;
	}
	return hash_value;
}
*/

//http://www.partow.net/programming/hashfunctions/
//failed
unsigned int DJBHash(const char* str, unsigned int length)
{
	unsigned int hash = 6727583;//5381;
	unsigned int i = 0;

	for (i = 0; i < length; ++str, ++i)
	{
		hash = ((hash << 5) + hash) + (*str);
	}

	return hash;
}
//failed
unsigned int RSHash(const char* str, unsigned int length)
{
	//6727583
	//6727741
	unsigned int b = 6727583;//378551;
	unsigned int a = 6727741;//63689;
	unsigned int hash = 0;
	unsigned int i = 0;

	for (i = 0; i < length; ++str, ++i)
	{
		hash = hash * a + (*str);
		a = a * b;
	}

	return hash;
}
//failed
unsigned int ELFHash(const char* str, unsigned int length)
{
	unsigned int hash = 0;
	unsigned int x = 0;
	unsigned int i = 0;

	for (i = 0; i < length; ++str, ++i)
	{
		hash = (hash << 4) + (*str);

		if ((x = hash & 0xF0000000L) != 0)
		{
			hash ^= (x >> 24);
		}

		hash &= ~x;
	}

	return hash;
}
//failed
long long int sfold(string s, int M) {
	int intLength = s.length() / 4;
	long sum = 0;
	for (int j = 0; j < intLength; j++) {
		const char *c = s.substr(j * 4, (j * 4) + 4).c_str();
		long mult = 1;

		int length = sizeof(c) / sizeof(char);
		for (int k = 0; k < length; k++) {
			sum += c[k] * mult;
			mult *= 256;
		}
	}

	const char *c = s.substr(intLength * 4).c_str();
	long mult = 1;
	int length = sizeof(c) / sizeof(char);
	for (int k = 0; k < length; k++) {
		sum += c[k] * mult;
		mult *= 256;
	}

	return(abs(sum) % M);
}

struct vertex_node
{
	vec3 position;
	vec3 normal;
	vec3 color;
	vec3 radiance;
	int triangle_index;
};

void add_unused_face(unordered_set<string>& unused_mtl)
{
	//phan tren nay ko giup gi duoc
	unused_mtl.insert("wire_125000000");
	unused_mtl.insert("wire_176026026");
	unused_mtl.insert("wire_154215229");
	unused_mtl.insert("comVases_decors_00");
	unused_mtl.insert("Line009");
	unused_mtl.insert("Object705440804");
	unused_mtl.insert("wire_225198087");
	unused_mtl.insert("Obj3d66_877474_90_12");
	unused_mtl.insert("Obj3d66_877474_61_715");
	unused_mtl.insert("Obj3d66_877474_62_106");
	unused_mtl.insert("Obj3d66_877474_63_424");
	unused_mtl.insert("Obj3d66_877474_65_596");
	unused_mtl.insert("Obj3d66_877474_66_57");
	unused_mtl.insert("Obj3d66_877474_67_290");
	unused_mtl.insert("Obj3d66_877474_68_506");
	unused_mtl.insert("Obj3d66_877474_80_344");
	unused_mtl.insert("Obj3d66_877474_81_711");
	unused_mtl.insert("Obj3d66_877474_83_437");
	//unused_mtl.insert("Shape004");
	unused_mtl.insert("xiadele_127_027_obj_11");
	unused_mtl.insert("xiadele_127_027_obj_12");
	unused_mtl.insert("xiadele_127_027_obj_13");
	//minoti lawrance nay co the nhin thay nhung khong ro
	unused_mtl.insert("Minotti_Minotti_Lawrence_Clan_Seating_4_107");
	unused_mtl.insert("Obj3d66_877474_82_327");
	unused_mtl.insert("Obj3d66_877474_86_667");
	unused_mtl.insert("Obj3d66_877474_87_684");
	unused_mtl.insert("Obj3d66_877474_91_297");
	unused_mtl.insert("Obj3d66_877474_92_420");
	unused_mtl.insert("Obj3d66_877474_84_274");
	unused_mtl.insert("Obj3d66_877474_94_591");
	unused_mtl.insert("Obj3d66_877474_95_891");
	unused_mtl.insert("Obj3d66_877474_96_757");
	unused_mtl.insert("Obj3d66_877474_97_400");
	unused_mtl.insert("Obj3d66_877474_98_282");
	unused_mtl.insert("Obj3d66_877474_99_928");
	unused_mtl.insert("Obj3d66_877474_100_241");
	unused_mtl.insert("Obj3d66_877474_101_903");
	unused_mtl.insert("Obj3d66_877474_102_920");
	unused_mtl.insert("Obj3d66_877474_85_046");
	unused_mtl.insert("Obj3d66_877474_93_211");
	unused_mtl.insert("Obj3d66_877474_86_666");
	unused_mtl.insert("Obj3d66_877474_87_683");
	unused_mtl.insert("Obj3d66_877474_80_343");
	unused_mtl.insert("Obj3d66_877474_81_710");
	unused_mtl.insert("Obj3d66_877474_82_326");
	unused_mtl.insert("Obj3d66_877474_83_436");
	unused_mtl.insert("Obj3d66_877474_84_273");
	unused_mtl.insert("Obj3d66_877474_85_45");
	unused_mtl.insert("Obj3d66_877474_173_633");
	unused_mtl.insert("Obj3d66_877474_174_511");
	unused_mtl.insert("Obj3d66_877474_175_717");
	unused_mtl.insert("Obj3d66_877474_176_807");
	unused_mtl.insert("Obj3d66_877474_177_934");
	unused_mtl.insert("wire_115115115");
	unused_mtl.insert("Arc006");
	unused_mtl.insert("wire_224198087");
	unused_mtl.insert("Md_House_168");
	unused_mtl.insert("Md_House_218");
	unused_mtl.insert("Md_House_219");
	unused_mtl.insert("Md_House_220");
	unused_mtl.insert("Vitra_S20_004");
	unused_mtl.insert("Vitra_S20_005");
	unused_mtl.insert("Vitra_S20_006");
	unused_mtl.insert("wire_153228184");
	//unused_mtl.insert("Md_House_083");
	unused_mtl.insert("Obj3d66_676235_47_503");
	unused_mtl.insert("wire_199199199");
	unused_mtl.insert("Cylinder887368");
	unused_mtl.insert("Cylinder887374");
	unused_mtl.insert("book_02");
	unused_mtl.insert("Obj3d66_981776_154_168");
	unused_mtl.insert("wire_027177148");
	unused_mtl.insert("ekko_r6_body074");
	unused_mtl.insert("ekko_r6_body075");
	unused_mtl.insert("wire_141171233");
	unused_mtl.insert("tiet110163");
	unused_mtl.insert("Rectangle2106816070");
	//unused_mtl.insert("Rectangle2106816073");

	//phan nay ko giup duoc gi
	unused_mtl.insert("wire_255255255");
	unused_mtl.insert("Obj3d66_848522_151_735");
	unused_mtl.insert("Obj3d66_848522_148_415");
	unused_mtl.insert("Plane2115543301");
	unused_mtl.insert("Plane2115543302");

	unused_mtl.insert("wire_213154229");
	unused_mtl.insert("wire_008008136");
	unused_mtl.insert("wire_228184153");
	unused_mtl.insert("wire_006135113");
	unused_mtl.insert("wire_177028149");
	
	unused_mtl.insert("wire_177088027");
	unused_mtl.insert("wire_138008110");
	unused_mtl.insert("wire_135110008");
	unused_mtl.insert("wire_113135006");
	unused_mtl.insert("wire_085028177");
	unused_mtl.insert("wire_225088199");

	unused_mtl.insert("Rectangle2106816069");
	unused_mtl.insert("PictureGlass");
	//unused_mtl.insert("Curtain");
	unused_mtl.insert("GlassDoor");

	//phan nam tren den am tran
	//unused_mtl.insert("Rectangle2106816065");

	//mai lon xon chua 2 cai den
	unused_mtl.insert("Rectangle2133390648");


	
	//Md_House_217

	//unused_mtl.insert("RectangleLight");

	//unused_mtl.insert("Rectangle006");
	//unused_mtl.insert("Rectangle005");

	//unused_mtl.insert("wire_184228153");


	//rem
	//unused_mtl.insert("wire_143224087");


	//phong ngoai
	//unused_mtl.insert("wire_006135058");

	unused_mtl.insert("ekko_r6_body161");

	unused_mtl.insert("SlideDoor006");

	unused_mtl.insert("clock_glass");

	//unused_mtl.insert("watch_glass");
}

void random_point_on_circle(float& r, float& x, float& z)
{
	float R = r * sqrt14(randf());

	float phi = tau * randf();

	float cos_phi, sin_phi;

	FTA::sincos(phi, &sin_phi, &cos_phi);

	x = R * cos_phi;
	z = R * sin_phi;
}

void random_point_on_sphere(float& r, vec3& point)
{
	float cos_theta = randf();
	float sin_theta = sqrt14(1.0f - cos_theta * cos_theta);

	float phi = tau * randf();

	float cos_phi, sin_phi;

	FTA::sincos(phi, &sin_phi, &cos_phi);

	point = vec3(r * sin_theta * cos_phi, r * sin_theta * sin_phi, r * cos_theta);
}

vec3 SphericalCoordinate(float& sin_theta, float& cos_theta, float& phi, vec3& x, vec3& y, vec3& z)
{
	float sin_phi, cos_phi;

	FTA::sincos(phi, &sin_phi, &cos_phi);

	return sin_theta * cos_phi * x + sin_theta * sin_phi * y + cos_theta * z;
}

vec3 sigma_s = vec3(0.2, 0.2, 0.2);
vec3 sigma_t = vec3(0.6f, 0.6f, 0.6f);

float g = -0.5f;

float abs_g = g > 0.0f ? g : -g;

bool isScatter(float& intersect_distance)
{
	int channel = minf(int(randf() * 3), 2);

	float sample_distance = -log(randf()) / sigma_t[channel];

	if (sample_distance < intersect_distance)
	{
		intersect_distance = sample_distance;
		return true;
	}
	return false;
}

//float phase_HG(float& cos_theta)
//float phase_HG(float& cos_theta)
float HG_pdf(vec3& wi, vec3& wo)
//float HG_pdf(float& cos_theta)
{
	//float denom = 1.0 + g * (g + 2.0 * cos_theta);

	float denom = 1.0 + g * (g + 2.0 * wi.dot(wo));

	return i4pi * (1.0 - g * g) / (denom * sqrt14(denom));
}

//HenyeyGreenStein Compute the sample angle cos_theta
//Phase_HG compute the probability cos_theta will be sample
//similar to bsdf and pdf
float HG_Sample(vec3& wi, vec3& wo)
{
	float cos_theta;

	float r1 = randf();

	if (abs_g < 1e-3)
		cos_theta = 1.0 - 2.0 * r1;
	else
	{
		float sqrt_term = (1.0 - g * g) / (1.0 - g + 2.0 * g * r1);

		cos_theta = (1.0 + g * g - sqrt_term * sqrt_term) / (2.0 * g);
	}

	float sin_theta = sqrt14(1.0 - cos_theta * cos_theta);

	float phi = tau * randf();

	onb onb(wi);

	wo = SphericalCoordinate(sin_theta, cos_theta, phi, onb.u, onb.v, onb.w);

	//return phase_HG(cos_theta);
	//return phase_HG(-wi, wo);

	return HG_pdf(-wi, wo);
	//return HG_pdf(cos_theta);
}

struct Scene
{
	Scene(const string& filename)
	{
		local_sphere = onb(vec3(0, -1, 0));

		round_light.resize(4);
		round_light_area.resize(4);

		add_unused_face(unused_mtl);
		ReadObj2(filename, "");

		compute_area();
	}
	Scene(const string& filename, string enviroment_path, float rotation)
	{
		local_sphere = onb(vec3(0, -1, 0));

		round_light.resize(4);
		round_light_area.resize(4);

		add_unused_face(unused_mtl);

		cout << "a\n";
		rotation_u = rotation;

		ReadObj2(filename, enviroment_path);

		compute_area();
	}

	Scene(const string& filename, vec3& look_from, string enviroment_path, float rotation)
	{
		local_sphere = onb(vec3(0, -1, 0));

		round_light.resize(4);
		round_light_area.resize(4);

		add_unused_face(unused_mtl);
		ReadObj_Remove_Back_Face(filename, look_from, enviroment_path);

		compute_area();
	}

	Scene(const string& filename, bool use_area_light_, bool use_directional_light_, vec3 light_direction_, float directional_blur_radius_ = 0.0f, bool use_sun_color_ = true, bool compute_sky_color_ = true, string enviroment_path = "", float rotation = 0, vec3& look_from = vec3(0), bool remove_back_face = false)
	{
		local_sphere = onb(vec3(0, -1, 0));

		round_light.resize(4);
		round_light_area.resize(4);

		add_unused_face(unused_mtl);

		cout << "b\n";
		rotation_u = rotation;
		use_area_light = use_area_light_;
		use_sun_color = use_sun_color_;
		compute_sky_color = compute_sky_color_;

		directional_blur_radius = directional_blur_radius_ - floorf(directional_blur_radius_);//only get the fractional part 

		blur_directional_light = directional_blur_radius > 0.0f;

		inifinite_light_direction = light_direction_.norm();
		use_directional_light = use_directional_light_;


		if (!remove_back_face)
			ReadObj2(filename, enviroment_path);
		else
			ReadObj_Remove_Back_Face(filename, look_from, enviroment_path);

		compute_area();
		/*if (enviroment_path != "")
		{
			compute_sky_color = false;
		}*/
	}

	Scene(const string& filename, vec3& delta_light_position_, vec3& delta_light_Ke_)
	{
		//vec3 c1 = vec3(618.993, 332.914, -5640.13);
		//vec3 c2 = vec3(658.968, 338.402, -5286.09);
		//vec3 c3 = vec3(929.277, 332.525, -5374.87);
		//vec3 c4 = vec3(928.717, 333.161, -5324.65);

		local_sphere = onb(vec3(0, -1, 0));

		round_light.resize(4);
		round_light_area.resize(4);

		add_unused_face(unused_mtl);

		delta_light_position = delta_light_position_;
		delta_light_Ke = delta_light_Ke_;

		ReadObj2(filename, "");

		compute_area();
	}
	/*
		531.382 2798.47 - 10924.8          x : 86.054
		617.436 2800.15 - 10852.4          z : 72.4

		2037.63 2798.47 - 10924.8          x : 72.5
		2110.13 2798.48 - 10852.4          z : 72.4

		3411.47 2798.47 - 10924.8          x : 72.5
		3483.97 2798.48 - 10852.4          z : 72.4

		3411.47 2798.47 - 8728.62          x : 86.06
		3497.53 2798.48 - 8656.18          z : 72.44
	*/

	/*vec3 c1 = vec3(574.409f, 2798.44f, -10888.6f);
	vec3 c2 = vec3(2073.88f, 2798.44f, -10888.6f);
	vec3 c3 = vec3(3447.72f, 2798.44f, -10888.6f);
	vec3 c4 = vec3(3454.5f, 2798.44f, -8692.4f);
	*/

	vec3 light_center[4] =
	{
		{574.409f, 2798.2f, -10888.6f},
		{2073.88f, 2798.2f, -10888.6f},
		{3447.72f, 2798.2f, -10888.6f},
		{3454.5f, 2798.2f, -8692.4f}
	};

	float light_radius = 36.2f;

	float pdf_light_radius = 1.0f / (4 * pi * light_radius * light_radius);

	//float limit_y = 2421.4804f;

	float cos_cone_light = 0.6f;

	//vec3 picture_light_position = vec3(2083.6, 2237.5, -11611.8);

	//up 600
	float limit_y = 2890.875732f;
	//up 700
	//float limit_y = 2990.875732f;

	//up 800
	//float limit_y = 3090.875732f;
	//up 1000
	//float limit_y = 3290.875732f;
	//up 1400
	//float limit_y = 3490.875732f;

	vec3 picture_light_position = //vec3(2063.422363, 2206.701904, -11416.070313);//vec3(2063.422363, 2406.701904, -11616.070313);

	//pos cone
	//vec3(2049.7150f, 2403.1689f, -11607.6386f);

	//forward 200
	//vec3(2049.7150f, 2403.1689f, -11807.6386f);
	//forward 200 up 200
	//vec3(2049.7150f, 2603.1689f, -11807.6386f);

	//pos 1
	//vec3(2062.427246f, 2290.875732f, -11534.795898f);


	//pos 2
	//move
	//vec3(2062.427246f, 2290.875732f, -11420.795898f);
	//up 200
	//vec3(2062.427246f, 2490.875732f, -11420.795898f);
	//move back 300
	//vec3(2062.427246f, 2290.875732f, -11234.795898f);
	//move back 400 up 400
	//vec3(2062.427246f, 2690.875732f, -11134.795898f);
	//move back 400 up 350
	//vec3(2062.427246f, 2640.875732f, -11134.795898f);
	//move back 400 up 400
	//vec3(2062.427246f, 2690.875732f, -11134.795898f);
	//move back 500 up 400
	//vec3(2062.427246f, 2690.875732f, -11034.795898f);
	//back 700 up 600
	//vec3(2062.427246f, limit_y, -10834.795898f);
	//back 800 up 600
	vec3(2062.427246f, limit_y, -10834.795898f);
	//back 800 up 700
	//vec3(2062.427246f, limit_y, -10734.795898f);
	//back 1000 up 800
	//vec3(2062.427246f, limit_y, -10534.795898f);
	//back 1000 up 1000
	//vec3(2062.427246f, limit_y, -10534.795898f);
	//back 1200 up 1000
	//vec3(2062.427246f, limit_y, -10334.795898f);
	//back 1400 up 1400
	//vec3(2062.427246f, limit_y, -10134.795898f);
	//back 1600 up 1400
	//vec3(2062.427246f, limit_y, -9934.795898f);

	//vec3(2075.0148f, 2421.4804f, -11634.6035f);

	float pic_light_radius = 6.0f;

	float pdf_light_cylinder = 1.0f / (4.0f * pi * pic_light_radius * pic_light_radius * pic_light_radius);

	//vector<vec3> v;
	//vector<vec2> vt;
	//vector<vec3> vn;

	vector<uint16_t> v;
	vector<int16_t> vt;
	vector<int16_t> vn;

	//vector<Triangle_save_space> trs;
	//vector<Triangle_24> trs;

	vector<Index_List> indices;
	vector<Triangle_16> trs;

	vector<Material> mats;
	vector<BSDF_Sampling*> bsdfs;

	vector<Texture> texs;

	//vector<vec3> sheen;

	//Texture sheen;

	//circle light
	vector<vector<int>> round_light;

	vector<vector<float>> round_light_area;

	//pipe light
	vec3 cylinder_right = vec3(2276, 2416, -11606);//vec3(2276.425293, 2419.020508, -11615.982422);
	vec3 cylinder_left = vec3(1843, 2416, -11606);//vec3(1843.658569, 2416.413086, -11606.794922)
	
	

	/*vec3 light_center[4] =
	{
		{ 618.993f, 332.914f, -5640.13f },
		{ 658.968f, 338.402f, -5286.09f },
		{ 929.277f, 332.525f, -5374.87f },
		{ 928.717f, 333.161f, -5324.65f }
	};*/

	//vec3 light_position = vec3(1618.993f, 260.914f, -4640.13f);
	
	//vec3 light_position = vec3(582.733, 2790.058, -10873);
	//vec3 light_position = vec3(1438, 1529, -8741);

	float round_light_pdf[4] = { 0.00065274151f, 0.00065316786f, 0.00065316786f, 0.00064432989f };

	vec3 round_light_Ke = vec3(1.0f, 0.4f, 0.0f);

	float sphere_radius = 60;

	onb local_sphere;
	/*
	light round 0 : 1532
	light round 1 : 1531
	light round 2 : 1531
	light round 3 : 1552
	*/


	float cylinder_length_x = 433.0f;

	float cylinder_radius = 6.0f;

	unsigned int first_shadow_node = 0;
	unsigned int second_shadow_node = 1;

	//AREA LIGHT
	AreaLight Ls;

	unordered_map<int, int> light_map;

	unordered_set<string> unused_mtl;

	//DIRECTIONAL LIGHT

	//default
	//sun color = (500, 400, 100)
	//sky color = (50, 80, 100)

	//these part are for directional light only
	//start of directional light variable
	bool use_area_light = true;
	bool use_directional_light = true;//false;
	bool use_sun_color = true;
	bool compute_sky_color = true;
	bool blur_directional_light = true;//false;


	float directional_blur_radius = 0.27f;//0.5f;//prev 0.3 //0.0f;
	float small_blur_radius = 0.027f;

	//253 243 198 x 4
	
	//vec3 sun_color = vec3(15.6, 15.2, 12.4);

	//vec3 sun_color = vec3(3.9686, 3.8117, 3.1058);
	
	//vec3 sun_color = vec3(3.9686, 3.1117, 2.8058);


	//244, 233, 155
	//0.9568 0.9137 0.6078

	//vec3 sun_color = vec3(5.7408, 5.4822, 3.6468);


	//252, 229, 112

	//vec3 sun_color = vec3(9.9215, 8.9803, 4.3921);
	
	//
	//vec3 sun_color = vec3(20 * 0.9882, 20 * 0.8980, 20 * 0.4392);
	vec3 sun_color = vec3(6 * 0.9882, 6 * 0.8980, 6 * 0.4392);

	//vec3 sun_color = vec3(6, 4.5, 4);
	
	//vec3 sun_color = vec3(0, 0, 0);

	//vec3 sun_color = vec3(30, 20, 10);

	//vec3 sun_color = vec3(3, 2.2, 1.2);

	//used in 2_6 directional
	//vec3 sun_color = vec3(1.5, 1.1, 0.6);

	//vec3 sun_color = vec3(3, 2.2, 1.2);

	
	//vec3 sun_color = vec3(0.65, 0.55, 0.35);

	//vec3 sun_color = vec3(2.6, 2.2, 1.4);

	//vec3 sun_color = vec3(50, 40, 10);//vec3(5, 5, 5);//(5, 4, 1);

	//Dinning Room use this color
	//vec3 sun_color = vec3(30, 30, 30);

	//vec3 inifinite_light_direction = vec3(1, 0.0, -0.2);
	
	vec3 inifinite_light_direction = vec3(1, 0.6, -0.5);//vec3(1, 0.0, -0.5);

	//end of directional light variable;

	//DELTA LIGHT OR POINT LIGHT

	vec3 delta_light_position;
	vec3 delta_light_Ke;

	int leaf = 4;
	bool use_eviroment = false;

	Enviroment Env;

	//float Env_Scale = 10.0f;

	float Env_Scale = 6.0f;

	//120 = best so far
	float rotation_u = 0.0f;// = -160 / 360.0f;

							//normal map is disable by default to save space
							//normal map is not always desirable since not having normal map result in better image
							//turn on to use
	bool use_bump_map = false;

	//normal mapping strength
	//default 1.0, increase this for greater deep

	float bum_strength = 2.0f;

	vec3 v_min = vec3(1e20);
	vec3 v_max = vec3(-1e20);
	vec3 v_extend;
	vec3 v_extend_i65535;

	string direction;

	void delete_before_render()
	{
		//vector<vec3>().swap(v);
		//vector<vec3>().swap(vn);
	}

	void delete_struct()
	{
		vector<uint16_t>().swap(v);
		vector<int16_t>().swap(vt);
		vector<int16_t>().swap(vn);
		vector<Triangle_16>().swap(trs);
		//vector<Triangle_24>().swap(trs);
		//vector<Triangle_save_space>().swap(trs);
		//vector<Triangle>().swap(trs);

		trs.shrink_to_fit();

		vector<Material>().swap(mats);
		for (auto& v : texs)
		{
			vector<float>().swap(v.a);
			vector<float>().swap(v.b);
			vector<uint8_t>().swap(v.c);
			vector<float>().swap(v.r);
		}

		vector<Texture>().swap(texs);
		vector<vec3>().swap(Env.c);
	}

	void Read_Sheen(const string& s)
	{
		int w, h, n;


		unsigned char* c = stbi_load((direction + s).c_str(), &w, &h, &n, 0);

		cout << w << " " << h << "\n";

		float a = 1.4f;
		if (n == 3)
		{
			Sheen_w = w;
			Sheen_h = h;

			Sheen_Albedo.resize(3 * w * h);

			for (int e = 0; e < 3 * w * h; ++e)
			{
				Sheen_Albedo[e] = c[e] / 255.0f;
			}
		}
		else if (n == 4)
		{
			Sheen_w = w;
			Sheen_h = h;
			

			Sheen_Albedo.resize(3 * w * h);

			int i = 0;
		
			for (int e = 0; e < 4 * w * h; e += 4)
			{
				Sheen_Albedo[i] = c[e];
				Sheen_Albedo[i + 1] = c[e + 1] / 255.0f;
				Sheen_Albedo[i + 2] = c[e + 2] / 255.0f;
				i += 3;
			}
		}
	}

	void ReadMtl(const string& filename, unordered_map<string, int>& mtlMap, const string& EnviromentName = "")//mtlMap chi dung cho giai doan doc file, ender thi khong can)
	{
		//only HDR file extension support

		if (EnviromentName != "")
		{
			cout << "use env\n";

			use_eviroment = true;
			ifstream f(direction + EnviromentName);

			int extension_position = EnviromentName.find_last_of(".");

			string extension = EnviromentName.substr(extension_position + 1);

			//cout << extension << "\n";

			/*if (extension != "hdr")
			{
			cout << "Unsupport Enviroment Extension\n";
			use_eviroment = false;
			}
			else
			{*/
			int w, h, n;

			unsigned char* c = stbi_load((direction + EnviromentName).c_str(), &w, &h, &n, 3);

			cout << direction + EnviromentName << "\n";

			vector<vec3> enviroment_value;

			enviroment_value.resize(w * h);

			for (int i = 0; i < w * h; ++i)
			{
				float x = c[3 * i] / 255.0f;
				float y = c[3 * i + 1] / 255.0f;
				float z = c[3 * i + 2] / 255.0f;

				enviroment_value[i] = vec3(Env_Scale * x, Env_Scale * y, Env_Scale * z);
			}

			Env.width = w;
			Env.height = h;
			Env.c = enviroment_value;

			Env.CDF(w, h, enviroment_value);

			free(c);

			vector<vec3>().swap(enviroment_value);
			//}
		}
		//cout << "texture path:" << direction + filename << "\n";
		ifstream f(direction + filename);
		if (!f)
			cout << "Mtl file not exist\n";



		char line[256];
		int numberOfMaterial = 0;
		while (f.getline(line, 256))
		{
			if (line[0] == 'n')
				numberOfMaterial++;
		}

		mats.resize(numberOfMaterial);
		bsdfs.resize(numberOfMaterial);
		f.clear();
		f.seekg(0, ios::beg);

		int countTexture = -1;
		int countMaterial = -1;

		char mtlname[64];
		char texturename[64];
		bool getdirection = true;
		unordered_map<string, int> temp_tex;

		Texture tex;

		while (f.getline(line, 256))
		{
			char* t = line;

			SkipSpace(t);

			if (strncmp(t, "newmtl", 6) == 0)
			{
				countMaterial++;
				sscanf_s(t += 7, "%s", mtlname);

				mtlMap[mtlname] = countMaterial;
				mats[countMaterial].Matname = mtlname;
			}
			else if (strncmp(t, "type", 4) == 0)
			{
				char type[64];
				sscanf_s(t += 5, "%s", type);
				if (strncmp(type, "Mirror", 6) == 0)
					mats[countMaterial].type = Mirror_type;
				else if (strncmp(type, "GlassColor", 10) == 0)
					mats[countMaterial].type = Glass_Color_type;
				else if (strncmp(type, "Glass", 5) == 0)
					mats[countMaterial].type = Glass_type;
				else if (strncmp(type, "GlassNoRefract", 14) == 0)
					mats[countMaterial].type = Glass_No_Refract_type;
				else if (strncmp(type, "RoughDielectric", 15) == 0)
					mats[countMaterial].type = Rough_Dielectric_type;
				else if (strncmp(type, "Conductor", 9) == 0)
				{
					mats[countMaterial].type = Conductor_type;
					t += 10;
					char complex_ior_name[64];
					sscanf_s(t, "%s", complex_ior_name);

					string ior_name = complex_ior_name;

					//cout << "Conductor Test\n";
					//cout << ior_name << "\n----------\n";

					for (int i = 0; i < 40; ++i)
					{
						if (ior_name == Ior_List[i].name)
						{
							mats[countMaterial].complex_ior_index = i;
							break;
						}
					}
				}

				else if (strncmp(type, "RoughConductor", 14) == 0)
				{
					//cout << "a\n";
					mats[countMaterial].type = RoughConductor_type;
					t += 15;
					char complex_ior_name[64];
					sscanf_s(t, "%s", complex_ior_name);//vi sao ko co t thi no ra GGX

					string ior_name = complex_ior_name;

					//cout << "test here!\n";
					//cout << complex_ior_name << "\n";
					//cout << ior_name << "\n";

					for (int i = 0; i < 40; ++i)
					{
						if (ior_name == Ior_List[i].name)
						{
							mats[countMaterial].complex_ior_index = i;
							break;
						}
					}
				}
				/*else if (strncmp(type, "RoughConductorSimple", 20) == 0)
				{
				//Al by deafult
				//not use any way
				mats[countMaterial].type = RoughConductorSimple_type;
				mats[countMaterial].complex_ior = Ior_List[2];
				}*/
				else if (strncmp(type, "PlasticSpecular", 15) == 0)
					mats[countMaterial].type = PlasticSpecular_type;
				else if (strncmp(type, "Plastic", 7) == 0)
					mats[countMaterial].type = Plastic_type;
				else if (strncmp(type, "RoughPlastic", 12) == 0)
					mats[countMaterial].type = RoughPlastic_type;
				else if (strncmp(type, "ThinSheet", 9) == 0)
					mats[countMaterial].type = ThinSheet_type;
				else if (strncmp(type, "ThinDielectric", 14) == 0)
				{
					//cout << "Im so thin!\n";
					mats[countMaterial].type = ThinDielectric_type;
				}
				else if (strncmp(type, "SmoothCoat", 10) == 0)
				{
					mats[countMaterial].type = SmoothCoat_type;
				}
				//else if (strncmp(type, "Transparent", 11) == 0)
				//{
				//	mats[countMaterial].type = Transparent_type;
				//}
				else if (strncmp(type, "LightDiffuse", 12) == 0)
					mats[countMaterial].type = Light_Diffuse_type;
				//else if (strncmp(type, "LightTransparent", 16) == 0)
					//{
					//	cout << "Trans";
				//	mats[countMaterial].type = Light_Transparent_type;
				//}
				else if (strncmp(type, "BackGround", 10) == 0)
					mats[countMaterial].type = Back_Ground_type;
				else if (strncmp(type, "Velvet", 6) == 0)
					mats[countMaterial].type = Velvet_type;
				else if (strncmp(type, "ClothDiffuse", 11) == 0)
					mats[countMaterial].type = Cloth_Diffuse_type;
				else if (strncmp(type, "Cloth", 5) == 0)
					mats[countMaterial].type = Cloth_type;
				//else if (strncmp(type, "Blinn", 5) == 0)
				//	mats[countMaterial].type = Blinn_type;
			}
			else if (strncmp(t, "ior", 3) == 0)
			{
				//float ior;
				//sscanf_s(t += 4, "%f", &ior);
				//mats[countMaterial].ior = ior;

				sscanf_s(t += 4, "%f", &mats[countMaterial].ior);
			}
			else if (strncmp(t, "gammar_texture", 14) == 0)
				sscanf_s(t += 15, "%f", &mats[countMaterial].gamma_texture);
			else if (strncmp(t, "roughness", 9) == 0)
			{
				sscanf_s(t += 10, "%f", &mats[countMaterial].roughness);
			}
			else if (strncmp(t, "thickness", 9) == 0)
			{
				sscanf_s(t += 10, "%f", &mats[countMaterial].thickness);
			}
			else if (strncmp(t, "sigma", 5) == 0)
				sscanf_s(t += 6, "%f", &mats[countMaterial].sigma);
			else if (strncmp(t, "spec_coeff", 10) == 0)
				sscanf_s(t += 11, "%f", &mats[countMaterial].specular_coefficient);
			else if (strncmp(t, "distribution", 12) == 0)
			{
				char distribution[32];
				sscanf_s(t += 13, "%s", distribution);

				if (strncmp(distribution, "GGX", 3) == 0)
					mats[countMaterial].distribution_type = GGX;
				if (strncmp(distribution, "Phong", 5) == 0)
					mats[countMaterial].distribution_type = Phong;
				if (strncmp(distribution, "Beckmann", 8) == 0)
					mats[countMaterial].distribution_type = Beckmann;
			}
			else if (t[0] == 'N')
			{
				if (t[1] == 'i')
				{
					sscanf_s(t += 3, "%f", &mats[countMaterial].Ni);
					float Ni = mats[countMaterial].Ni;
					float R0 = (Ni - 1) / (Ni + 1);
					//mats[countMaterial].R0 = R0 * R0;
				}
				else if (t[1] == 's')
					sscanf_s(t += 3, "%f", &mats[countMaterial].Ns);
			}
			//else if (t[0] == 'd')
			//	sscanf_s(t += 2, "%d", &mats[countMaterial].d);
			/*else if (t[0] == 'T')
			{
				if (t[1] == 'r')
				{
					mats[countMaterial].alpha_material = true;
					float Tr;
					sscanf_s(t += 3, "%f", &Tr);
					//cout << Tr << "\n";
					mats[countMaterial].Tr = Tr;
					//mats[countMaterial].opaque = 1.0f - Tr;
				}
				else if (t[1] == 'f')
				{
					//this can be used as well for color filter
					//however, all color channel are mostly filter with the same value and mostly 1 
					//so I only use single coefficient for maximum perfomance
					
					mats[countMaterial].alpha_material = true;
					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Tf.x, &mats[countMaterial].Tf.y, &mats[countMaterial].Tf.z);
					mats[countMaterial].Tr = 1.0f;
				}
			}*/
			else if (strncmp(t, "TransGirl", 9) == 0)
			{
				mats[countMaterial].alpha_material = true;
				sscanf_s(t += 10, "%f", &mats[countMaterial].Tr);
			}
			else if (t[0] == 'i')
			{
				int illum;
				sscanf_s(t += 6, "%d", &illum);
				mats[countMaterial].illum = illum;
				//mats[countMaterial].type = Diffuse;//|= illum == 7 ? Fres : illum == 5 ? PRefl : NonS;
			}
			else if (t[0] == 'K')
			{
				//if (t[1] == 'a')
				//	sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Ka.x, &mats[countMaterial].Ka.y, &mats[countMaterial].Ka.z);
				//else if (t[1] == 'd')
				if(t[1] == 'd')
				{

					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Kd.x, &mats[countMaterial].Kd.y, &mats[countMaterial].Kd.z);
				}
				else if (t[1] == 's')
				{
					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Ks.x, &mats[countMaterial].Ks.y, &mats[countMaterial].Ks.z);
				}
				else if (t[1] == 'e')
				{
					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Ke.x, &mats[countMaterial].Ke.y, &mats[countMaterial].Ke.z);
					if (mats[countMaterial].Ke.maxc() > 0)
					{
						//cout << "Light Detected! \n";

						//cout << mats[countMaterial].Matname << "\n";
						mats[countMaterial].isLight = true;
					}
					else
					{
						mats[countMaterial].isLight = false;
					}
				}

			}

			else if (strncmp(t, "map_Kd", 6) == 0)
			{
				mats[countMaterial].use_texture = true;
				sscanf_s(t += 7, "%s", texturename);
				string realname;

				getfilename(string(texturename), realname);

				if (temp_tex.find(realname) == temp_tex.end())
				{
					//cout << realname << "\n";
					//cout << direction << "\n";
					int w, h, n;


					unsigned char* c = stbi_load((direction + realname).c_str(), &w, &h, &n, 0);

					//cout << direction + realname << "\n";

					/*
					if (n == 3 || n == 4)
					{
						tex.w = w;
						tex.h = h;
						tex.n = 3;

						tex.c.resize(3 * w * h);

						for (int e = 0; e < 3 * w * h; ++e)
						{
							tex.c[e] = 255 * powf(c[e] / 255.0f, 2.6f);
						}

						if (n == 4)
						{
							//check alpha material
							float first_alpha_value = c[3] / 255.0f;
							
							bool use_alpha_texture = false;

							for (int i = 7; i < 4 * w * h; i += 4)
							{
								float second_alpha_value = c[i] / 255.0f;

								if (first_alpha_value != second_alpha_value)
								{
									use_alpha_texture = true;
									break;
								}
							}

							if (use_alpha_texture)
							{
								//cout << "Alpha Texture!\n";
								mats[countMaterial].use_alpha_texture = true;

								tex.a.resize(w * h);

								tex.w2 = w;
								tex.h2 = h;
								tex.n2 = 1;

								for (int i = 0; i < w * h; ++i)
									tex.a[i] = c[4 * i + 3] / 255.0f;
							}
						}
					}
					else if(n == 1)
					{
						tex.w = w;
						tex.h = h;
						tex.n = 3;

						tex.c.resize(3 * w * h);

						for (int e = 0; e < n * w * h; ++e)
						{
							tex.c[3 * e] = tex.c[3 * e + 1] = tex.c[3 * e + 2] = 255 * powf(c[e] / 255.0f, 2.6f);
						}
					}
					*/

					float a = 1.4f;
					if (n == 3)
					{
						tex.w = w;
						tex.h = h;
						tex.n = 3;

						tex.c.resize(3 * w * h);

						float b = a;
						if (realname == "Md_House_Cloth5.jpg")
							b = 2.0f;
							//b = 0.7f;
						if (realname == "background_cloth.jpg")
							b = 0.7f;

						if (realname == "dark_blue_velvet.jpg" || realname == "velvet.jpg")
							b = 1.6f;
							//b = 2.6f;
						
						if (realname == "SC_Enzo_Bed_difuse2_bk2.jpg")
							b = 2.4;

						//if (realname == "clock.jpg")
						//	b = 6.0f;

						for (int e = 0; e < 3 * w * h; ++e)
						{
							tex.c[e] = 255 * powf(c[e] / 255.0f, b);
						}
					}
					else if (n == 4)
					{
						tex.w = w;
						tex.h = h;
						tex.n = 3;

						tex.c.resize(3 * w * h);

						int s = 0;
						for (int e = 0; e < 4 * w * h; e += 4)
						{
							tex.c[s] = 255 * powf(c[e] / 255.0f, a);
							tex.c[s + 1] = 255 * powf(c[e + 1] / 255.0f, a);
							tex.c[s + 2] = 255 * powf(c[e + 2] / 255.0f, a);
							s += 3;
						}

						//check alpha material
						float first_alpha_value = c[3] / 255.0f;

						bool use_alpha_texture = false;

						for (int i = 7; i < 4 * w * h; i += 4)
						{
							float second_alpha_value = c[i] / 255.0f;

							if (first_alpha_value != second_alpha_value)
							{
								use_alpha_texture = true;
								break;
							}
						}
						
						

						if (use_alpha_texture)
						{
							//cout << "Alpha Texture!\n";
							mats[countMaterial].use_alpha_texture = true;

							tex.a.resize(w * h);

							//tex.w2 = w;
							//tex.h2 = h;
							//tex.n2 = 1;
							//cout << realname << "\n";
							//bay gio chi chua transparent
							for (int i = 0; i < w * h; ++i)
							{
								float alpha = c[4 * i + 3] / 255.0f;
								//if (alpha > 0.0f && alpha != 1.0f)
								//	alpha = 0.0f;

								/*if (alpha >= 0.9f)
									alpha = 1.0f;
								else
									alpha = 0.0f;*/

										//cout << alpha << "\n";
								//tex.a[i] = c[4 * i + 3] / 255.0f;
								
								//if (alpha > 0.0f && alpha != 1.0f)
								//{
								//tex.a[i] = 1.0f - (c[4 * i + 3] / 255.0f);// > 0.0f ? 1.0f : 0.0f;
								tex.a[i] = 1.0f - alpha;
								//cout << tex.a[i] << "\n";
								//}
							}
						}
					}
					
					else if (n == 1)
					{
						tex.w = w;
						tex.h = h;
						tex.n = 3;

						tex.c.resize(3 * w * h);

						for (int e = 0; e < n * w * h; ++e)
						{
							tex.c[3 * e] = tex.c[3 * e + 1] = tex.c[3 * e + 2] = 255 * powf(c[e] / 255.0f, a);
						}
					}


					//vector<float>().swap(tex.c);

					//tex.c = value;

					free(c);
					//free(value);


					//t.Texname = realname;
					++countTexture;

					texs.emplace_back(tex);
					temp_tex[realname] = countTexture;
					mats[countMaterial].TexInd = countTexture;
					//mats[countMaterial].Texname = realname;
				}
				else
				{
					int ind = temp_tex[realname];
					mats[countMaterial].TexInd = ind;
				}
			}

			else if (strncmp(t, "map_d", 5) == 0)
			{
				//mats[countMaterial].alpha_material = true;
				mats[countMaterial].use_alpha_texture = true;
				sscanf_s(t += 6, "%s", texturename);
				string realname;

				getfilename(string(texturename), realname);

				int w2, h2, n2;
				//stbi_load((direction + realname).c_str(), &w2, &h2, &n2, 3);

				unsigned char* a = stbi_load((direction + realname).c_str(), &w2, &h2, &n2, 0);

				//tex.w2 = w2;
				//tex.h2 = h2;
				//tex.n2 = n2;

				tex.a.resize(n2 * w2 * h2);

				//int i = 0;
				for (int e = 0; e < n2 * w2 * h2; ++e)
				{
					tex.a[e] = a[e] / 255.0f > 0.1f ? 1.0f : 0.0f;
				}

				free(a);
				texs[countTexture] = tex;
			}

			else if (use_bump_map && strncmp(t, "map_bump", 8) == 0)
			{
				mats[countMaterial].use_Bump = true;
				sscanf(t += 9, "%s", texturename);
				string realname;

				getfilename(string(texturename), realname);

				//cout << realname << "\n";

				int w3, h3, n3;
				//stbi_load((direction + realname).c_str(), &w3, &h3, &n3, 3);

				unsigned char* bump = stbi_load((direction + realname).c_str(), &w3, &h3, &n3, 0);

				//cout << n3 << "\n";

				//tex.w3 = w3;
				//tex.h3 = h3;
				tex.n3 = 3;

				tex.b.resize(3 * w3 * h3);

				vector<float> gray_scale;

				gray_scale.resize(w3 * h3);

				bool is_normal_map;

				if (n3 == 1)
				{
					is_normal_map = false;
					for (int e = 0; e < w3 * h3; ++e)
					{
						gray_scale[e] = bump[e] / 255.0f;
					}
				}
				else
				{
					//n3 = 3 or n3 = 4 can be a bump map or normal map
					//a texture is normal map if 90% of its pixel have dominance blue channel
					//bump map is always greyscale and will never pass this test
					//check if this texture is bump map or normal map
					//if it is bump map, we have to convert to norml map

					//check Bump Map or Normal Map
					int count = 0;
					for (int i = 0; i < w3; ++i)
					{
						for (int j = 0; j < h3; ++j)
						{
							int index = (j * w3 + i);

							float r = bump[n3 * index] / 255.0f;
							float g = bump[n3 * index + 1] / 255.0f;
							float b = bump[n3 * index + 2] / 255.0f;

							gray_scale[index] = 0.299f * r + 0.587f * g + 0.114f * b;

							if (b > g && b > r)
								++count;
						}
					}
					//truly normal map
					if (count > 0.9 * w3 * h3)
					{
						is_normal_map = true;
						for (int i = 0; i < w3 * h3; ++i)
						{
							float r = bump[n3 * i] / 255.0f;
							float g = bump[n3 * i + 1] / 255.0f;
							float b = bump[n3 * i + 2] / 255.0f;

							tex.b[3 * i] = r;
							tex.b[3 * i + 1] = g;
							tex.b[3 * i + 2] = b;
						}
					}
					else
						is_normal_map = false;
				}
				//here this is a bump map convert it to a normal map before storing
				if (!is_normal_map)
				{
					vector<vec3> normal_map = Bump_Map_To_Normal_Map(w3, h3, n3, bum_strength, gray_scale);

					for (int i = 0; i < w3 * h3; ++i)
					{
						tex.b[3 * i] = normal_map[i].x;
						tex.b[3 * i + 1] = normal_map[i].y;
						tex.b[3 * i + 2] = normal_map[i].z;
					}
				}

				vector<float>().swap(gray_scale);
				free(bump);
				texs[countTexture] = tex;
			}
			else if (strncmp(t, "map_Kr", 6) == 0)
			{
				//cout << "Kr here\n";
				mats[countMaterial].use_Roughness = true;
				sscanf_s(t += 7, "%s", texturename);
				string realname;

				getfilename(string(texturename), realname);

				if (temp_tex.find(realname) == temp_tex.end())
				{
					//cout << realname << "\n";

					int w4, h4, n4;


					unsigned char* r = stbi_load((direction + realname).c_str(), &w4, &h4, &n4, 0);
					int count = 0;
					for (int i = 0; i < w4 * h4; ++i)
					{
						float x = r[n4 * i] / 255.0f;
						float y = r[n4 * i + 1] / 255.0f;
						float z = r[n4 * i + 2] / 255.0f;
						//if (x != y || y != z || z != x)
						//	++count;
						//cout << x << " " << y << " " << z << "\n";
					}

					//tex.w4 = w4;
					//tex.h4 = h4;
					tex.n4 = 1;

					tex.r.resize(w4 * h4);

					for (int e = 0; e < n4 * w4 * h4; ++e)
					{
						tex.r[e] = r[e] / 255.0f;
					}


					//vector<float>().swap(tex.c);

					//tex.c = value;

					free(r);
					//free(value);


					//t.Texname = realname;
					++countTexture;

					texs.emplace_back(tex);
					temp_tex[realname] = countTexture;
					mats[countMaterial].TexInd = countTexture;
					//mats[countMaterial].Texname = realname;
				}
				else
				{
					int ind = temp_tex[realname];
					mats[countMaterial].TexInd = ind;
				}
			}

			else if (strncmp(t, "enable_interference", 19) == 0)
			{
				char b[8];
				sscanf_s(t += 20, "s", b);

				if (strncmp(b, "false", 5) == 0)
					mats[countMaterial].use_InterFerence = false;
				if (strncmp(b, "true", 4) == 0)
					mats[countMaterial].use_InterFerence = true;
			}
		}


		for (int i = 0; i < numberOfMaterial; ++i)
		{
			if (!mats[i].isLight)
			{
				if (mats[i].type == Diffuse_type)
					bsdfs[i] = new Diffuse();
				else if (mats[i].type == Mirror_type)
				{
					bsdfs[i] = new Mirror();
				}
				else if (mats[i].type == Glass_type)
					//bsdfs[i] = new Glass(mats[i].ior);
					//bsdfs[i] = new Glass_Fresnel_Rough(mats[i].ior, 0.02);

					//bsdfs[i] = new Rough_Dielectric_Original2_Final_V(mats[i].ior, mats[i].roughness);

					//bsdfs[i] = new ThinDielectric(mats[i].ior);

					//failed
					//bsdfs[i] = new ThinSheet2();
					//bsdfs[i] = new ThinSheet(1.5f, 1.0f, 0.0f, true);

					//bsdfs[i] = new Glass_Green_Tint(mats[i].ior);
					//bsdfs[i] = new Glass_Green(mats[i].ior);
					//bsdfs[i] = new Glass(mats[i].ior);
					//bsdfs[i] = new Glass_Schlick(mats[i].ior);


					bsdfs[i] = new Glass_Fresnel(mats[i].ior);
				else if (mats[i].type == Glass_Color_type)
					bsdfs[i] = new Glass_Color(mats[i].Kd, mats[i].ior);
				else if (mats[i].type == Glass_No_Refract_type)
					bsdfs[i] = new Glass_No_Refract(mats[i].ior);
				//bsdfs[i] = new Glass_Beer_Law(mats[i].ior);
				else if (mats[i].type == Rough_Dielectric_type)
					bsdfs[i] = new Rough_Glass_Fresnel(mats[i].ior, mats[i].roughness);

				//bsdfs[i] = new Glass_Fresnel(mats[i].ior);

				//bsdfs[i] = new Glass(mats[i].ior);

				//bsdfs[i] = new Glass_Fresnel_Rough(mats[i].ior, mats[i].roughness);



				//bsdfs[i] = new ThinDielectric(mats[i].ior);

				//bsdfs[i] = new Rough_Dielectric_Benedik(mats[i].ior, mats[i].roughness);

				//bsdfs[i] = new Rough_Dielectric_Paper(mats[i].ior, mats[i].roughness);

				//bsdfs[i] = new Rough_Dielectric_Sang(mats[i].roughness, mats[i].ior);

				//bsdfs[i] = new Mirror();
				//bsdfs[i] = new Rough_Dielectric_Original2_Final_V(mats[i].ior, mats[i].roughness);

				//bsdfs[i] = new Rough_Dielectric_Original_Good_Upgrade(1.5f, 0.12f);

				//bsdfs[i] = new Rough_Dielectric(1.5, 0.12);

				//bsdfs[i] = new Rough_Dielectric(mats[i].ior, mats[i].roughness);

				//bsdfs[i] = new Rough_Dielectric_Y2(mats[i].ior, mats[i].roughness);

				//bsdfs[i] = new Rough_Dielectric_Y2(mats[i].roughness, mats[i].ior);
				else if (mats[i].type == Conductor_type)
					bsdfs[i] = new Conductor();
				else if (mats[i].type == RoughConductor_type || mats[i].type == RoughConductorComplex_type)
				{
					//bsdfs[i] = new Rough_Conductor_TungSten_Absolute(mats[i].roughness, mats[i].complex_ior);


					//bsdfs[i] = new Final_Rough_Conductor_Absolute(mats[i].roughness, mats[i].complex_ior);
					//bsdfs[i] = new Rough_Conductor_TungSten_fix_black(mats[i].roughness, mats[i].complex_ior);
					//bsdfs[i] = new Final_Rough_Conductor(mats[i].roughness, mats[i].complex_ior);

					//if (!mats[i].use_Roughness)
					//{
					//bsdfs[i] = new Diffuse();
					bsdfs[i] = new Rough_Conductor_TungSten(mats[i].roughness, Ior_List[mats[i].complex_ior_index]);
					/*}
					else
					{
					//cout << "here\n";
					//bsdfs[i] = new Final_Rough_Conductor_Roughness_Mapping(mats[i].complex_ior);
					bsdfs[i] = new Final_Rough_Conductor(mats[i].roughness, mats[i].complex_ior);

					}*/

					//bsdfs[i] = new Rough_Conductor_Yoan_Blien(mats[i].roughness)
					//bsdfs[i] = new Final_Rough_Conductor_modify(mats[i].roughness, mats[i].complex_ior);



					//bsdfs[i] = new Conductor(mats[i].complex_ior);





					//bsdfs[i] = new Rough_Dielectric_Y(mats[i].roughness, mats[i].ior);

					//bsdfs[i] = new Rough_Dielectric_Sang(mats[i].roughness, mats[i].ior);

					//bsdfs[i] = new Rough_Dielectric(mats[i].roughness, mats[i].ior);



					//bsdfs[i] = new Rough_Conductor_Visible_Normal(mats[i].roughness, mats[i].complex_ior);



					//bsdfs[i] = new Glossy(mats[i].roughness, mats[i].complex_ior);

					//bsdfs[i] = new GGX_Material(mats[i].ior, mats[i].roughness);


					//bsdfs[i] = new Rough_Conductor(mats[i].roughness, mats[i].complex_ior);

					//bsdfs[i] = new Rough_Conductor_Paper(mats[i].roughness, mats[i].complex_ior);


				}
				else if (mats[i].type == RoughConductorSimple_type)
				{
					//final_rough_conductor chi tinh theo color
					//khong tinh theo complex_ior
					bsdfs[i] = new Final_Rough_Conductor(mats[i].roughness, Ior_List[mats[i].complex_ior_index]);
				}
				else if (mats[i].type == Plastic_type)
				{
					//bsdfs[i] = new Plastic_Complex(mats[i].ior, mats[i].thickness);
					
					//bsdfs[i] = new Plastic_Diffuse_Thickness(mats[i].ior);

					bsdfs[i] = new Plastic_Specular(mats[i].ior);


					//bsdfs[i] = new Plastic_Diffuse(mats[i].ior);

					//bsdfs[i] = new New_Plastic(mats[i].roughness, mats[i].distribution_type, mats[i].ior, mats[i].Ks);

					//bsdfs[i] = new Plastic_Complex(mats[i].ior, 0.0f);
				}
				else if (mats[i].type == RoughPlastic_type)
				{
					//bsdfs[i] = new Rough_Plastic_Tungsten(mats[i].ior, mats[i].roughness, mats[i].thickness);

					//Good
					bsdfs[i] = new Rough_Plastic_Diffuse(mats[i].ior, mats[i].roughness);

					if(mats[i].Matname == "Curtain")
						bsdfs[i] = new Rough_Plastic_Diffuse(mats[i].ior, mats[i].roughness, vec3(0.9f, 0.9f, 0.9f));

					//bsdfs[i] = new Rough_Plastic_Diffuse_Thickness(mats[i].ior, mats[i].roughness);

					//bsdfs[i] = new RoughPlastic_Original_Diffuse_Version(mats[i].ior, mats[i].roughness);

					//bsdfs[i] = new Diffuse();

					//bsdfs[i] = new RoughPlastic_Original_Glossy_Version(mats[i].ior, mats[i].roughness);

					//bsdfs[i] = new Rough_Plastic_Specular(mats[i].ior, mats[i].roughness);


					//bsdfs[i] = new Rough_Plastic_Specular(mats[i].ior, mats[i].roughness);


					//float d = mats[i].Kd.maxc();
					//float s = mats[i].Ks.maxc();

					//float wd = d / (d + s);
					//float ws = 1.0f - wd;


					//bsdfs[i] = new Rough_Plastic(mats[i].ior, mats[i].roughness);		
				}

				else if (mats[i].type == ThinDielectric_type)
					bsdfs[i] = new ThinDielectric(mats[i].ior);
				else if (mats[i].type == ThinSheet_type)
					bsdfs[i] = new ThinSheet(mats[i].ior, 1.0f, 0.0f, false);
				//bsdfs[i] = new ThinDielectric(mats[i].ior);
				//bsdfs[i] = new Glass_Fresnel(mats[i].ior);
				else if (mats[i].type == SmoothCoat_type)
				{
					//BSDF_Sampling* sub_material = new Diffuse();//Plastic_Diffuse(1.9f);
					//bsdfs[i] = new SmoothCoat(mats[i].ior, mats[i].sigma, sub_material);

					bsdfs[i] = new SmoothCoatTungsten();
				}
				else if (mats[i].type == Transparent_type)
				{
					mats[i].alpha_material = true;
					bsdfs[i] = new Transparent(mats[i].Tr);
				}
				else if (mats[i].type == Light_Diffuse_type)
				{
					bsdfs[i] = new Diffuse();//new Light_Diffuse(mats[i].Kd);
				}
				else if(mats[i].type == Velvet_type)
					bsdfs[i] = new Velvet();
				else if(mats[i].type == Cloth_Diffuse_type)
					bsdfs[i] = new Cloth_Diffuse();
				else if (mats[i].type == Cloth_type)
				{
					//if (mats[i].Matname == "Curtain")
					//	bsdfs[i] = new DieletricLambert();
					
					//chay duoc
					//else
						bsdfs[i] = new Ashkhmin2();
					//bsdfs[i] = new Velvet();
				}
					//hoi dom
					//bsdfs[i] = new Cloth2();


					//dom
					//bsdfs[i] = new Cloth_Albedo_Scalling();

					//dom
					//bsdfs[i] = new Corona_Cloth();

					//dom
					//bsdfs[i] = new Ashkmin_Shirley2();

					//chay duoc nhung do
					//bsdfs[i] = new Anisotropic_Cloth();
					
					//bsdfs[i] = new Diffuse();
					//
					//bsds[i] = new Ashkhmin();
					//bsdfs[i] = new Cloth();
					
			}
		}
	}
	/*
	void ReadObj(const string& p, const string& enviroment_path)
	{
		string filename;
		getdDirectionAndName(p, direction, filename);
		getdirection(p, direction);


		ifstream f(p);
		if (!f)
			cout << "obj file not exist";
		char line[1024], name[256];
		char mtlib[64];
		vec3 vtn;

		int mtl = -1;
		bool isLight = false;
		int count = 0;
		unordered_map<string, int> mtlMap;
		bool read_mtl = false;

		int count_light = 0;
		int num_v = 0;
		int num_vt = 0;
		int num_vn = 0;

		int fi = -1;
		int li = 0;

		vec3 v_min(1e20);
		vec3 v_max(-1e20);

		//vec3 v1(1e20);
		//vec3 v2(-1e20);


		int count_triangle = 0;
		int count_exceed_127 = 0;

		int v_max_diff = -1e9;
		int vt_max_diff = -1e9;
		int vn_max_diff = -1e9;

		unordered_map<string, int> index_map;

		int ind = 0;
		//first pass only read f
		//become map become so big it will eat up all the memory
		//second pass read v, vt, vn
		
		while (f.getline(line, 1024))
		{
			char* t = line;
			int prev_space = strspn(t, " \t");
			t += prev_space;

			if (strncmp(t, "v", 1) == 0)
			{
				t += 1;

				if (strncmp(t, " ", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					float x, y, z;
					sscanf_s(t, "%f %f %f", &x, &y, &z);

					//cout << x << " " << y << " " << z << "\n";
					vec3 vertex(x, y, z);

					v.emplace_back(vec3(x, y, z));
					//cout << v.size() << "\n";

					++num_v;

					v_min = min_vec(v_min, vertex);
					v_max = max_vec(v_max, vertex);
				}
				else if (strncmp(t, "t", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					float x, y;
					sscanf_s(t, "%f %f", &x, &y);
					//vt.emplace_back(vec2(x, y));
					
					int16_t vtx = (int16_t)(x * 10000);
					int16_t vty = (int16_t)(y * 10000);

					vt.emplace_back(vtx);
					vt.emplace_back(vty);

					++num_vt;
				}
				else if (strncmp(t, "n", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					float x, y, z;
					sscanf_s(t, "%f %f %f", &x, &y, &z);
					
					int16_t vnx = (int16_t)(x * 10000);
					int16_t vny = (int16_t)(y * 10000);
					int16_t vnz = (int16_t)(z * 10000);
				

					vn.emplace_back(vnx);
					vn.emplace_back(vny);
					vn.emplace_back(vnx);
				
					
					++num_vn;
				}
			}
			else if (strncmp(t, "f", 1) == 0)
			{
				if (count_triangle > 0 && count_triangle % 100000 == 0)
					cout << count_triangle << "\n";
				//t += 1;
				int post_space = strspn(t + 1, " \t");

				t += post_space + 1;


				int face_type;
				vector<Face> faces;
				get_face_index(t, num_v, num_vt, num_vn, faces);

				count_triangle += faces.size();

				for (int i = 0; i < faces.size(); ++i)
				{

					++fi;
					int v0 = faces[i].v[0];
					int v1 = faces[i].v[1];
					int v2 = faces[i].v[2];

					int abs_v0_v1 = abs(v0 - v1);
					int abs_v1_v2 = abs(v1 - v2);
					int abs_v2_v0 = abs(v2 - v0);

					if (abs_v0_v1 > v_max_diff)
						v_max_diff = abs_v0_v1;

					if (abs_v1_v2 > v_max_diff)
						v_max_diff = abs_v1_v2;

					if (abs_v2_v0 > v_max_diff)
						v_max_diff = abs_v2_v0;

					int vt0 = faces[i].vt[0];
					int vt1 = faces[i].vt[1];
					int vt2 = faces[i].vt[2];

					int abs_vt0_vt1 = abs(vt0 - vt1);
					int abs_vt1_vt2 = abs(vt1 - vt2);
					int abs_vt2_vt0 = abs(vt2 - vt0);

					if (abs_vt0_vt1 > vt_max_diff)
						vt_max_diff = abs_vt0_vt1;

					if (abs_vt1_vt2 > vt_max_diff)
						vt_max_diff = abs_vt1_vt2;

					if (abs_vt2_vt0 > vt_max_diff)
						vt_max_diff = abs_vt2_vt0;

					int vn0 = faces[i].vn[0];
					int vn1 = faces[i].vn[1];
					int vn2 = faces[i].vn[2];

					int abs_vn0_vn1 = abs(vn0 - vn1);
					int abs_vn1_vn2 = abs(vn1 - vn2);
					int abs_vn2_vn0 = abs(vn2 - vn0);

					if (abs_vn0_vn1 > vn_max_diff)
						vn_max_diff = abs_vn0_vn1;

					if (abs_vn1_vn2 > vn_max_diff)
						vn_max_diff = abs_vn1_vn2;

					if (abs_vn2_vn0 > vn_max_diff)
						vn_max_diff = abs_vn2_vn0;

					if (abs_v0_v1 > 127 || abs_v1_v2 > 127 || abs_v2_v0 > 127 ||
						abs_vt0_vt1 > 127 || abs_vt1_vt2 > 127 || abs_vt2_vt0 > 127 ||
						abs_vn0_vn1 > 127 || abs_vn1_vn2 > 127 || abs_vn2_vn0 > 127)

						++count_exceed_127;

					string s0 = to_string(v0) + "|" + to_string(vt0) + "|" + to_string(vn0);
					string s1 = to_string(v1) + "|" + to_string(vt1) + "|" + to_string(vn1);
					string s2 = to_string(v2) + "|" + to_string(vt2) + "|" + to_string(vn2);

					

					int ind0, ind1, ind2;

					if (index_map.find(s0) != index_map.end())
						ind0 = index_map[s0];
					else
					{
						indices.emplace_back(Index_List(v0, vt0, vn0));
						ind0 = ind;
						index_map[s0] = ind++;
					}
					if (index_map.find(s1) != index_map.end())
						ind1 = index_map[s1];
					else
					{
						indices.emplace_back(Index_List(v1, vt1, vn1));
						ind1 = ind;
						index_map[s1] = ind++;
					}
					if (index_map.find(s2) != index_map.end())
						ind2 = index_map[s2];
					else
					{
						indices.emplace_back(Index_List(v2, vt2, vn2));
						ind2 = ind;
						index_map[s2] = ind++;
					}

					Triangle_16 tri(ind0, ind1, ind2, mtl);
					trs.emplace_back(tri);

					//Triangle_24 tri(v0, v1, v2, vt0, vt1, vt2, vn0, vn1, vn2, mtl);
					//trs.emplace_back(tri);

					//Triangle_save_space tri(v0, v1, v2, vt0, vt1, vt2, vn0, vn1, vn2, mtl);
					//trs.emplace_back(tri);

				}
			}
			


			else if (strncmp(t, "usemtl", 6) == 0)//(t[0] == 'u')
			{
				sscanf_s(t += 7, "%s", name);
				string realname = (string)name;

				mtl = mtlMap[realname];

				if (mtlMap.find(realname) == mtlMap.end())
					cout << realname << " not exist in obj file\n";

			}
			else if (t[0] == 'm' && !read_mtl)//mtllib
			{
				sscanf_s(t += 7, "%s", mtlib);

				string realname = (string)mtlib;

				ReadMtl(realname, mtlMap, enviroment_path);
				//ReadMtl2(realname, mtlMap);// , enviroment_path);
				read_mtl = true;
			}
		}

		index_map.clear();

		/*cout << "High Poly SanMiguel: \n";

		cout << "Num Triangle: " << count_triangle << "\n";

		cout << "v size : " << num_v << "\n";
		cout << "vt size: " << num_vt << "\n";
		cout << "vn size: " << num_vn << "\n";

		cout << "v  max different: " << v_max_diff << "\n";
		cout << "vt max different: " << vt_max_diff << "\n";
		cout << "vn max different: " << vn_max_diff << "\n";

		cout << "Num Triangles Exceed: " << count_exceed_127 << "\n";
		cout << "Memory Require If All is 24 bytes                            : " << 24 * count_triangle << "\n";
		cout << "Memory Require If 16 bytes for most plus 16 bytes for exceed : " << 16 * count_triangle + 16 * count_exceed_127 << "\n";
		*/

		/*
		Low Poly SanMiguel:
		Num Triangle: 5617451
		v size : 3738829
		vt size: 844670
		vn size: 4517249
		v  max different: 82798
		vt max different: 89030
		vn max different: 200026
		Memory Require If All is 24 bytes                            : 134818824
		Memory Require If 16 bytes for most plus 16 bytes for exceed : 108393504
		Reading time: 53.502s
		36
		0
		v  size: 3738829
		vt size: 844670
		vn size: 4517249
		Building time: 0.004s
		Delete used data
		*/

		/*
		High Poly SanMiguel:
		Num Triangle: 9980699
		v size : 5933233
		vt size: 2221669
		vn size: 6721532
		v  max different: 82798
		vt max different: 89030
		vn max different: 200026
		Num Triangles Exceed: 2011819
		Memory Require If All is 24 bytes                            : 239536776
		Memory Require If 16 bytes for most plus 16 bytes for exceed : 191880288
		Reading time: 70.264s
		36
		0
		v  size: 5933233
		vt size: 2221669
		vn size: 6721532
		Building time: 0.003s
		Delete used data
		*/


	//}
	
	void ReadObj2(const string& p, const string& enviroment_path)
	{
		string filename;
		getdDirectionAndName(p, direction, filename);
		getdirection(p, direction);


		ifstream f(p);
		if (!f)
			cout << "obj file not exist";
		char line[1024], name[256];
		char mtlib[64];
		vec3 vtn;

		int mtl = -1;
		bool isLight = false;
		int count = 0;
		unordered_map<string, int> mtlMap;
		bool read_mtl = false;

		int count_light = 0;
		int num_v = 0;
		int num_vt = 0;
		int num_vn = 0;

		int fi = -1;
		int li = 0;

		//vec3 v_min(1e20);
		//vec3 v_max(-1e20);

		//vec3 v1(1e20);
		//vec3 v2(-1e20);


		int count_triangle = 0;
	
		//int v_max_diff = -1e9;
		//int vt_max_diff = -1e9;
		//int vn_max_diff = -1e9;

		unordered_map<string, int> index_map;

		int ind = 0;
		//first pass only read f
		//become map become so big it will eat up all the memory
		//second pass read v, vt, vn

		//int max_v_vn_diff = -1e20f;

		bool add_face = true;

		vector<vec3> current_v;

		vec3 max1 = vec3(INT_MIN, INT_MIN, INT_MIN), min1 = vec3(INT_MAX, INT_MAX, INT_MAX);
		vec3 max2 = vec3(INT_MIN, INT_MIN, INT_MIN), min2 = vec3(INT_MAX, INT_MAX, INT_MAX);
		vec3 max3 = vec3(INT_MIN, INT_MIN, INT_MIN), min3 = vec3(INT_MAX, INT_MAX, INT_MAX);
		vec3 max4 = vec3(INT_MIN, INT_MIN, INT_MIN), min4 = vec3(INT_MAX, INT_MAX, INT_MAX);

		string material_name = "";

		while (f.getline(line, 1024))
		{
			char* t = line;
			int prev_space = strspn(t, " \t");
			t += prev_space;
			if (strncmp(t, "v ", 2) == 0)
			{

				t += 2;
				int post_space = strspn(t, " \t");
				t += post_space;

				float x, y, z;
				sscanf_s(t, "%f %f %f", &x, &y, &z);
			
				vec3 vertex(x, y, z);
		
				current_v.push_back(vertex);

				v_min = min_vec(v_min, vertex);
				v_max = max_vec(v_max, vertex);

				num_v++;
			}
			else if (strncmp(t, "vt", 2) == 0)
			{
				t += 2;
				int post_space = strspn(t, " \t");
				t += post_space;

				float x, y;
				sscanf_s(t, "%f %f", &x, &y);

				int16_t vtx = (int16_t)(x * 10000);
				int16_t vty = (int16_t)(y * 10000);

				vt.emplace_back(vtx);
				vt.emplace_back(vty);

				//vt.emplace_back(vec2(x, y));
				++num_vt;
			}
			else if (strncmp(t, "vn", 2) == 0)
			{
				t += 2;
				int post_space = strspn(t, " \t");
				t += post_space;

				float x, y, z;
				sscanf_s(t, "%f %f %f", &x, &y, &z);

				int16_t vnx = (int16_t)(x * 10000);
				int16_t vny = (int16_t)(y * 10000);
				int16_t vnz = (int16_t)(z * 10000);


				vn.emplace_back(vnx);
				vn.emplace_back(vny);
				vn.emplace_back(vnz);

				//vn.emplace_back(vec3(x, y, z));
				++num_vn;
			}
			else if (strncmp(t, "f", 1) == 0 && add_face)
			{
				//if (count_triangle > 0 && count_triangle % 100000 == 0)
				//	cout << count_triangle << "\n";
				//t += 1;
				int post_space = strspn(t + 1, " \t");

				t += post_space + 1;

				int face_type;
				vector<Face> faces;
				get_face_index(t, num_v, num_vt, num_vn, faces);

				//count_triangle += faces.size();

				for (int i = 0; i < faces.size(); ++i)
				{

					++fi;
					int v0 = faces[i].v[0];
					int v1 = faces[i].v[1];
					int v2 = faces[i].v[2];

					vec3 vertex1 = current_v[v0];
					vec3 vertex2 = current_v[v1];
					vec3 vertex3 = current_v[v2];
					if (material_name == "ekko_r6_body146")
					{
						min1 = min_vec(min1, vertex1);
						max1 = max_vec(max1, vertex1);

						min1 = min_vec(min1, vertex2);
						max1 = max_vec(max1, vertex2);

						min1 = min_vec(min1, vertex3);
						max1 = max_vec(max1, vertex3);
					}
					else if (material_name == "ekko_r6_body169")
					{
						min2 = min_vec(min2, vertex1);
						max2 = max_vec(max2, vertex1);

						min2 = min_vec(min2, vertex2);
						max2 = max_vec(max2, vertex2);

						min2 = min_vec(min2, vertex3);
						max2 = max_vec(max2, vertex3);
					}
					else if (material_name == "ekko_r6_body171")
					{
						min3 = min_vec(min3, vertex1);
						max3 = max_vec(max3, vertex1);

						min3 = min_vec(min3, vertex2);
						max3 = max_vec(max3, vertex2);

						min3 = min_vec(min3, vertex3);
						max3 = max_vec(max3, vertex3);
					}
					else if (material_name == "ekko_r6_body172")
					{
						min4 = min_vec(min4, vertex1);
						max4 = max_vec(max4, vertex1);

						min4 = min_vec(min4, vertex2);
						max4 = max_vec(max4, vertex2);

						min4 = min_vec(min4, vertex3);
						max4 = max_vec(max4, vertex3);
					}
					
					int vt0 = faces[i].vt[0];
					int vt1 = faces[i].vt[1];
					int vt2 = faces[i].vt[2];

					

					int vn0 = faces[i].vn[0];
					int vn1 = faces[i].vn[1];
					int vn2 = faces[i].vn[2];

					//int abs_v0_vn0 = abs(v0 - vn0);
					//int abs_v1_vn1 = abs(v1 - vn1);
					//int abs_v2_vn2 = abs(v2 - vn2);

					//max_v_vn_diff = maxf(max_v_vn_diff, abs_v0_vn0);
					//max_v_vn_diff = maxf(max_v_vn_diff, abs_v1_vn1);
					//max_v_vn_diff = maxf(max_v_vn_diff, abs_v1_vn1);

					string str0 = to_string(v0) + "|" + to_string(vt0);
					string str1 = to_string(v1) + "|" + to_string(vt1);
					string str2 = to_string(v2) + "|" + to_string(vt2);


					int ind0, ind1, ind2;

					if (index_map.find(str0) != index_map.end())
					{
						ind0 = index_map[str0];

						int vn_in_map = indices[ind0].vn;

						int vn_current = vn0;
						
						//vn_in_map == vn_current ind0 = index_map[str0]
						//so I only need to check vn_in_map != vn_current case
						
						if (vn_in_map != vn_current)
						{
							indices.emplace_back(Index_List(v0, vt0, vn0));
							ind0 = ind++;
							//index_map[str0] = ind++;
						}

					}
					else
					{
						indices.emplace_back(Index_List(v0, vt0, vn0));
						ind0 = ind;
						index_map[str0] = ind++;
					}
					
					if (index_map.find(str1) != index_map.end())
					{
						ind1 = index_map[str1];

						int vn_in_map = indices[ind1].vn;

						int vn_current = vn1;

						//vn_in_map == vn_current ind0 = index_map[str0]
						//so I only need to check vn_in_map != vn_current case

						if (vn_in_map != vn_current)
						{
							indices.emplace_back(Index_List(v1, vt1, vn1));
							ind1 = ind++;
							//index_map[str0] = ind++;
						}
					}
					else
					{
						indices.emplace_back(Index_List(v1, vt1, vn1));
						ind1 = ind;
						index_map[str1] = ind++;
					}
					if (index_map.find(str2) != index_map.end())
					{
						ind2 = index_map[str2];

						int vn_in_map = indices[ind2].vn;

						int vn_current = vn2;

						//vn_in_map == vn_current ind0 = index_map[str0]
						//so I only need to check vn_in_map != vn_current case

						if (vn_in_map != vn_current)
						{
							indices.emplace_back(Index_List(v2, vt2, vn2));
							ind2 = ind++;
							//index_map[str0] = ind++;
						}
					}
					else
					{
						indices.emplace_back(Index_List(v2, vt2, vn2));
						ind2 = ind;
						index_map[str2] = ind++;
					}

					Triangle_16 tri(ind0, ind1, ind2, mtl);

					

					trs.emplace_back(tri);

					//Triangle_24 tri(v0, v1, v2, vt0, vt1, vt2, vn0, vn1, vn2, mtl);
					//trs.emplace_back(tri);

					//Triangle_save_space tri(v0, v1, v2, vt0, vt1, vt2, vn0, vn1, vn2, mtl);
					//trs.emplace_back(tri);

				}
			}



			else if (strncmp(t, "usemtl", 6) == 0)//(t[0] == 'u')
			{
				sscanf_s(t += 7, "%s", name);
				string realname = (string)name;

				add_face = !unused_mtl.count(realname);

				mtl = mtlMap[realname];

				//material_name = mats[mtl].Matname;

				/*if (material_name == "ekko_r6_body146")
				{
					cout << "ahha";
				}*/

				if (mtlMap.find(realname) == mtlMap.end())
					cout << realname << " not exist in obj file\n";

			}
			else if (t[0] == 'm' && !read_mtl)//mtllib
			{
				sscanf_s(t += 7, "%s", mtlib);

				string realname = (string)mtlib;

				ReadMtl(realname, mtlMap, enviroment_path);
				//ReadMtl2(realname, mtlMap);// , enviroment_path);
				read_mtl = true;
			}
		}

		/*cout << min1.x << " " << min1.y << " " << min1.z << "\n";
		cout << max1.x << " " << max1.y << " " << max1.z << "\n";

		cout << min2.x << " " << min2.y << " " << min2.z << "\n";
		cout << max2.x << " " << max2.y << " " << max2.z << "\n";

		cout << min3.x << " " << min3.y << " " << min3.z << "\n";
		cout << max3.x << " " << max3.y << " " << max3.z << "\n";

		cout << min4.x << " " << min3.y << " " << min4.z << "\n";
		cout << max4.x << " " << max3.y << " " << max4.z << "\n";*/
		/*c1 = (min1 + max1) * 0.5f;
		c2 = (min2 + max2) * 0.5f;
		c3 = (min3 + max3) * 0.5f;
		c4 = (min4 + max4) * 0.5f;

		cout << "c1: " << c1.x << " " << c1.y << " " << c1.z << "\n";
		cout << "c2: " << c2.x << " " << c2.y << " " << c2.z << "\n";
		cout << "c3: " << c3.x << " " << c3.y << " " << c3.z << "\n";
		cout << "c4: " << c4.x << " " << c4.y << " " << c4.z << "\n";

		cout << "light round 0 : " << round_light[0].size() << "\n";
		cout << "light round 1 : " << round_light[1].size() << "\n";
		cout << "light round 2 : " << round_light[2].size() << "\n";
		cout << "light round 3 : " << round_light[3].size() << "\n";*/

		/*for (int i = 0; i < trs.size(); ++i)
		{
			int mtl = trs[i].mtl;

			string s = mats[mtl].Matname;

			if (s == "ekko_r6_body146")
				round_light[0].emplace_back(i);
			
			else if (mats[mtl].Matname == "ekko_r6_body169")		
				round_light[1].emplace_back(i);		
			else if (mats[mtl].Matname == "ekko_r6_body171")		
				round_light[2].emplace_back(i);			
			else if (mats[mtl].Matname == "ekko_r6_body172")		
				round_light[3].emplace_back(i);	
		}*/

		index_map.clear();

		unordered_map<string, int>().swap(index_map);

		index_map.rehash(0);

		f.clear();
		f.seekg(0, ios::beg);

		//cout << "max_diff v vn: " << max_v_vn_diff << "\n";

		cout << "vmin: " << v_min.x << " " << v_min.y << " " << v_min.z << "\n";
		cout << "vmax: " << v_max.x << " " << v_max.y << " " << v_max.z << "\n";

		v_extend = v_max - v_min;
		v_extend_i65535 = v_extend / 65535.0f;
		
		vec3 scale = vec3(65535.0f) / v_extend;

		cout << "Start Transforming vertex data: \n";

		for (auto& c : current_v)
		{
			vec3 vertex_scale = (c - v_min) * scale;

			uint16_t vx = (uint16_t)(vertex_scale.x);
			uint16_t vy = (uint16_t)(vertex_scale.y);
			uint16_t vz = (uint16_t)(vertex_scale.z);

			v.emplace_back(vx);
			v.emplace_back(vy);
			v.emplace_back(vz);
		}

		vector<vec3>().swap(current_v);

		//vec3 box_extend = v_max - v_min;

		//vec3 inv_dimension(1.0f / v_extend.x, 1.0f / v_extend.y, 1.0f / v_extend.z);
		
		//vec3 scale = vec3(65535.0f) * inv_dimension;

		

		
		/*while (f.getline(line, 1024))
		{
			char* t = line;
			int prev_space = strspn(t, " \t");
			t += prev_space;

			if (strncmp(t, "v", 1) == 0)
			{
				t += 1;

				if (strncmp(t, " ", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					float x, y, z;
					sscanf_s(t, "%f %f %f", &x, &y, &z);
				
					vec3 vertex(x, y, z);
				

					vec3 vertex_scale = (vertex - v_min) * scale;
					
					uint16_t vx = (uint16_t)(vertex_scale.x);
					uint16_t vy = (uint16_t)(vertex_scale.y);
					uint16_t vz = (uint16_t)(vertex_scale.z);

					v.emplace_back(vx);
					v.emplace_back(vy);
					v.emplace_back(vz);
				}
				else if (strncmp(t, "t", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					float x, y;
					sscanf_s(t, "%f %f", &x, &y);
					
					int16_t vtx = (int16_t)(x * 10000);
					int16_t vty = (int16_t)(y * 10000);

					vt.emplace_back(vtx);
					vt.emplace_back(vty);

					//vt.emplace_back(vec2(x, y));
					++num_vt;
				}
				else if (strncmp(t, "n", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					float x, y, z;
					sscanf_s(t, "%f %f %f", &x, &y, &z);
					
					int16_t vnx = (int16_t)(x * 10000);
					int16_t vny = (int16_t)(y * 10000);
					int16_t vnz = (int16_t)(z * 10000);


					vn.emplace_back(vnx);
					vn.emplace_back(vny);
					vn.emplace_back(vnz);

					//vn.emplace_back(vec3(x, y, z));
					++num_vn;
				}
			}
		}*/

		/*cout << "High Poly SanMiguel: \n";

		cout << "Num Triangle: " << count_triangle << "\n";

		cout << "v size : " << num_v << "\n";
		cout << "vt size: " << num_vt << "\n";
		cout << "vn size: " << num_vn << "\n";

		cout << "v  max different: " << v_max_diff << "\n";
		cout << "vt max different: " << vt_max_diff << "\n";
		cout << "vn max different: " << vn_max_diff << "\n";

		cout << "Num Triangles Exceed: " << count_exceed_127 << "\n";
		cout << "Memory Require If All is 24 bytes                            : " << 24 * count_triangle << "\n";
		cout << "Memory Require If 16 bytes for most plus 16 bytes for exceed : " << 16 * count_triangle + 16 * count_exceed_127 << "\n";
		*/

		/*
		Low Poly SanMiguel:
		Num Triangle: 5617451
		v size : 3738829
		vt size: 844670
		vn size: 4517249
		v  max different: 82798
		vt max different: 89030
		vn max different: 200026
		Memory Require If All is 24 bytes                            : 134818824
		Memory Require If 16 bytes for most plus 16 bytes for exceed : 108393504
		Reading time: 53.502s
		36
		0
		v  size: 3738829
		vt size: 844670
		vn size: 4517249
		Building time: 0.004s
		Delete used data
		*/

		/*
		High Poly SanMiguel:
		Num Triangle: 9980699
		v size : 5933233
		vt size: 2221669
		vn size: 6721532
		v  max different: 82798
		vt max different: 89030
		vn max different: 200026
		Num Triangles Exceed: 2011819
		Memory Require If All is 24 bytes                            : 239536776
		Memory Require If 16 bytes for most plus 16 bytes for exceed : 191880288
		Reading time: 70.264s
		36
		0
		v  size: 5933233
		vt size: 2221669
		vn size: 6721532
		Building time: 0.003s
		Delete used data
		*/


	}

	void ReadObj_Remove_Back_Face(const string& p, vec3& look_from, const string& enviroment_path)
	{
		string filename;
		getdDirectionAndName(p, direction, filename);
		getdirection(p, direction);


		ifstream f(p);
		if (!f)
			cout << "obj file not exist";
		char line[1024], name[256];
		char mtlib[64];
		vec3 vtn;

		int mtl = -1;
		bool isLight = false;
		
		unordered_map<string, int> mtlMap;
		bool read_mtl = false;

		int count_light = 0;
		int num_v = 0;
		int num_vt = 0;
		int num_vn = 0;

		int fi = -1;
		int li = 0;

		//vec3 v_min(1e20);
		//vec3 v_max(-1e20);

		//vec3 v1(1e20);
		//vec3 v2(-1e20);


		int count_triangle = 0;

		//int v_max_diff = -1e9;
		//int vt_max_diff = -1e9;
		//int vn_max_diff = -1e9;

		unordered_map<string, int> index_map;

		int ind = 0;
		//first pass only read f
		//become map become so big it will eat up all the memory
		//second pass read v, vt, vn

		int max_v_vn_diff = -1e20f;

		while (f.getline(line, 1024))
		{
			char* t = line;
			int prev_space = strspn(t, " \t");
			t += prev_space;
			if (strncmp(t, "v ", 2) == 0)
			{

				t += 1;
				int post_space = strspn(t, " \t");
				t += post_space;

				float x, y, z;
				sscanf_s(t, "%f %f %f", &x, &y, &z);

				vec3 vertex(x, y, z);

				v_min = min_vec(v_min, vertex);
				v_max = max_vec(v_max, vertex);

				num_v++;
			}
			else if (strncmp(t, "f", 1) == 0)
			{
				if (count_triangle > 0 && count_triangle % 100000 == 0)
					cout << count_triangle << "\n";
				//t += 1;
				int post_space = strspn(t + 1, " \t");

				t += post_space + 1;


				int face_type;
				vector<Face> faces;
				get_face_index(t, num_v, num_vt, num_vn, faces);

				count_triangle += faces.size();

				for (int i = 0; i < faces.size(); ++i)
				{

					++fi;
					int v0 = faces[i].v[0];
					int v1 = faces[i].v[1];
					int v2 = faces[i].v[2];


					int vt0 = faces[i].vt[0];
					int vt1 = faces[i].vt[1];
					int vt2 = faces[i].vt[2];



					int vn0 = faces[i].vn[0];
					int vn1 = faces[i].vn[1];
					int vn2 = faces[i].vn[2];

					int abs_v0_vn0 = abs(v0 - vn0);
					int abs_v1_vn1 = abs(v1 - vn1);
					int abs_v2_vn2 = abs(v2 - vn2);

					max_v_vn_diff = maxf(max_v_vn_diff, abs_v0_vn0);
					max_v_vn_diff = maxf(max_v_vn_diff, abs_v1_vn1);
					max_v_vn_diff = maxf(max_v_vn_diff, abs_v1_vn1);

					string str0 = to_string(v0) + "|" + to_string(vt0);
					string str1 = to_string(v1) + "|" + to_string(vt1);
					string str2 = to_string(v2) + "|" + to_string(vt2);


					int ind0, ind1, ind2;

					if (index_map.find(str0) != index_map.end())
					{
						ind0 = index_map[str0];

						int vn_in_map = indices[ind0].vn;

						int vn_current = vn0;

						//vn_in_map == vn_current ind0 = index_map[str0]
						//so I only need to check vn_in_map != vn_current case

						if (vn_in_map != vn_current)
						{
							indices.emplace_back(Index_List(v0, vt0, vn0));
							ind0 = ind++;
							//index_map[str0] = ind++;
						}

					}
					else
					{
						indices.emplace_back(Index_List(v0, vt0, vn0));
						ind0 = ind;
						index_map[str0] = ind++;
					}

					if (index_map.find(str1) != index_map.end())
					{
						ind1 = index_map[str1];

						int vn_in_map = indices[ind1].vn;

						int vn_current = vn1;

						//vn_in_map == vn_current ind0 = index_map[str0]
						//so I only need to check vn_in_map != vn_current case

						if (vn_in_map != vn_current)
						{
							indices.emplace_back(Index_List(v1, vt1, vn1));
							ind1 = ind++;
							//index_map[str0] = ind++;
						}
					}
					else
					{
						indices.emplace_back(Index_List(v1, vt1, vn1));
						ind1 = ind;
						index_map[str1] = ind++;
					}
					if (index_map.find(str2) != index_map.end())
					{
						ind2 = index_map[str2];

						int vn_in_map = indices[ind2].vn;

						int vn_current = vn2;

						//vn_in_map == vn_current ind0 = index_map[str0]
						//so I only need to check vn_in_map != vn_current case

						if (vn_in_map != vn_current)
						{
							indices.emplace_back(Index_List(v2, vt2, vn2));
							ind2 = ind++;
							//index_map[str0] = ind++;
						}
					}
					else
					{
						indices.emplace_back(Index_List(v2, vt2, vn2));
						ind2 = ind;
						index_map[str2] = ind++;
					}

					Triangle_16 tri(ind0, ind1, ind2, mtl);
					trs.emplace_back(tri);

					//Triangle_24 tri(v0, v1, v2, vt0, vt1, vt2, vn0, vn1, vn2, mtl);
					//trs.emplace_back(tri);

					//Triangle_save_space tri(v0, v1, v2, vt0, vt1, vt2, vn0, vn1, vn2, mtl);
					//trs.emplace_back(tri);

				}
			}



			else if (strncmp(t, "usemtl", 6) == 0)//(t[0] == 'u')
			{
				sscanf_s(t += 7, "%s", name);
				string realname = (string)name;

				mtl = mtlMap[realname];

				if (mtlMap.find(realname) == mtlMap.end())
					cout << realname << " not exist in obj file\n";

			}
			else if (t[0] == 'm' && !read_mtl)//mtllib
			{
				sscanf_s(t += 7, "%s", mtlib);

				string realname = (string)mtlib;

				ReadMtl(realname, mtlMap, enviroment_path);
				//ReadMtl2(realname, mtlMap);// , enviroment_path);
				read_mtl = true;
			}
		}

		index_map.clear();

		unordered_map<string, int>().swap(index_map);

		index_map.rehash(0);

		f.clear();
		f.seekg(0, ios::beg);

		cout << "max_diff v vn: " << max_v_vn_diff << "\n";

		cout << "vmin: " << v_min.x << " " << v_min.y << " " << v_min.z << "\n";
		cout << "vmax: " << v_max.x << " " << v_max.y << " " << v_max.z << "\n";

		v_extend = v_max - v_min;
		v_extend_i65535 = v_extend / 65535.0f;

		//vec3 box_extend = v_max - v_min;

		//vec3 inv_dimension(1.0f / v_extend.x, 1.0f / v_extend.y, 1.0f / v_extend.z);

		//vec3 scale = vec3(65535.0f) * inv_dimension;

		vec3 scale = vec3(65535.0f) / v_extend;

		cout << "Start Reading vertex data: \n";
		while (f.getline(line, 1024))
		{
			char* t = line;
			int prev_space = strspn(t, " \t");
			t += prev_space;

			if (strncmp(t, "v", 1) == 0)
			{
				t += 1;

				if (strncmp(t, " ", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					float x, y, z;
					sscanf_s(t, "%f %f %f", &x, &y, &z);

					vec3 vertex(x, y, z);


					vec3 vertex_scale = (vertex - v_min) * scale;

					uint16_t vx = (uint16_t)(vertex_scale.x);
					uint16_t vy = (uint16_t)(vertex_scale.y);
					uint16_t vz = (uint16_t)(vertex_scale.z);

					v.emplace_back(vx);
					v.emplace_back(vy);
					v.emplace_back(vz);
				}
				else if (strncmp(t, "t", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					float x, y;
					sscanf_s(t, "%f %f", &x, &y);

					int16_t vtx = (int16_t)(x * 10000);
					int16_t vty = (int16_t)(y * 10000);

					vt.emplace_back(vtx);
					vt.emplace_back(vty);

					//vt.emplace_back(vec2(x, y));
					++num_vt;
				}
				else if (strncmp(t, "n", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					float x, y, z;
					sscanf_s(t, "%f %f %f", &x, &y, &z);

					int16_t vnx = (int16_t)(x * 10000);
					int16_t vny = (int16_t)(y * 10000);
					int16_t vnz = (int16_t)(z * 10000);


					vn.emplace_back(vnx);
					vn.emplace_back(vny);
					vn.emplace_back(vnz);

					//vn.emplace_back(vec3(x, y, z));
					++num_vn;
				}
			}
		}

		cout << "Remove Back Face:\n ";
	
		int s = trs.size();

		int count_front_face = 0;

		for (int i = 0; i < s; ++i)
		{
			int ind0 = trs[i].ind_v0;
			int ind1 = trs[i].ind_v1;
			int ind2 = trs[i].ind_v2;

			int v0 = indices[ind0].v;
			int v1 = indices[ind1].v;
			int v2 = indices[ind2].v;

			vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
			vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
			vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);

			
			p0 *= v_extend_i65535;
			p1 *= v_extend_i65535;
			p2 *= v_extend_i65535;

			p0 += v_min;
			p1 += v_min;
			p2 += v_min;

			vec3 normal = (p1 - p0).cross(p2 - p0);

			vec3 center = (p0 + p1 + p1) * 0.3333f;

			if (normal.dot(center - look_from) <= 0.0f)
				trs[count_front_face++] = trs[i];
		}

		auto first = trs.cbegin() + count_front_face;
		auto last = trs.cbegin() + trs.size();

		int count_back_face = trs.size() - count_front_face;

		trs.erase(first, last);
		cout << "Number of Back Face Triangles Removed: " << count_back_face << "\n";
	}

	void sah(node*& n, const int& start, const int& range, vector<BBox>& boxes, vector<SimpleTriangle>& simp)
	{
		n = new node(start, range);//range
		for (auto i = start; i < start + range; ++i)//<start + range
		{
			int ind = simp[i].i;
			simp[i].b = boxes[ind];
			n->box.expand(boxes[ind]);
		}
		if (range < leaf)
		{
			n->leaf = 1;
			n->axis = n->box.maxDim();
			int axis = n->axis;
			sort(simp.begin() + start, simp.begin() + start + range, [axis](const SimpleTriangle& s1, const SimpleTriangle& s2)
			{
				return s1.c[axis] < s2.c[axis];
			});
			return;
		}
		else
		{
			n->leaf = 0;
			n->range = 0;
			int best_split = 0, best_axis = -1;

			float best_cost = ci * range;
			float area = n->box.area();
			vec3 vmin = n->box.bbox[0], vmax = n->box.bbox[1];
			vec3 extend(vmax - vmin);

			for (int a = 0; a < 3; ++a)
			{
				sort(simp.begin() + start, simp.begin() + start + range, [a](const SimpleTriangle& s1, const SimpleTriangle& s2)
				{
					return s1.c[a] <= s2.c[a];
				});

				//float min_box = n->box.bbox[0][a], length = n->box.bbox[1][a] - min_box;
				float length = n->box.bbox[1][a] - n->box.bbox[0][a];
				//if (length < 0.000001f)
				//	continue;

				vector<BBox> right_boxes;
				right_boxes.resize(range);

				BBox left = simp[start + 0].b;

				right_boxes[range - 1] = simp[start + range - 1].b;

				for (int j = range - 2; j >= 0; --j)
				{
					right_boxes[j] = right_boxes[j + 1].expand_box(simp[start + j].b);
				}

				float extend = length / range;

				int count_left = 1;
				int count_right = range - 1;
				float inv_a = 1.0f / area;
				for (int i = 0; i < range - 1; ++i)
				{
					float left_area = left.area();
					float right_area = right_boxes[i + 1].area();

					BBox right = right_boxes[i + 1];

					float cost = ct + ci * (left_area * count_left + right_area * count_right) * inv_a;

					if (cost < best_cost)
					{
						best_cost = cost;
						best_axis = a;
						best_split = count_left;
					}
					++count_left;
					--count_right;
					left.expand(simp[start + i + 1].b);
				}
			}
			if (best_cost == ci * range)
			{
				n->leaf = 1;
				n->range = range;
				n->axis = n->box.maxDim();
				int axis = n->axis;
				sort(simp.begin() + start, simp.begin() + start + range, [axis](const SimpleTriangle& s1, const SimpleTriangle& s2)
				{
					return s1.c[axis] < s2.c[axis];
				});
				return;
			}

			n->axis = best_axis;

			sort(simp.begin() + start, simp.begin() + start + range, [best_axis](const SimpleTriangle& s1, const SimpleTriangle& s2)
			{
				return s1.c[best_axis] < s2.c[best_axis];
			});

			sah(n->nodes[0], start, best_split, boxes, simp);

			sah(n->nodes[1], start + best_split, range - best_split, boxes, simp);
		}
	}

	
	void sah_bin(node*& n, const int& start, const int& range, vector<BBox>& boxes, vector<SimpleTriangle>& simp)
	{
		n = new node(start, range);//range
		for (auto i = start; i < start + range; ++i)//<start + range
		{
			int ind = simp[i].i;
			simp[i].b = boxes[ind];
			n->box.expand(boxes[ind]);
		}
		if (range < leaf)
		{
			n->leaf = 1;
			n->axis = n->box.maxDim();
			int axis = n->axis;
			sort(simp.begin() + start, simp.begin() + start + range, [axis](const SimpleTriangle& s1, const SimpleTriangle& s2)
			{
				return s1.c[axis] < s2.c[axis];
			});
			return;
		}
		else
		{
			n->leaf = 0;
			n->range = 0;
			int best_split = 0, best_axis = -1;

			float best_cost = ci * range;
			float area = n->box.area();

			float inv_area = 1.0f / area;

			vec3 vmin = n->box.bbox[0], vmax = n->box.bbox[1];
			vec3 extend(vmax - vmin);

			for (int a = 0; a < 3; ++a)
			{
				sort(simp.begin() + start, simp.begin() + start + range, [a](const SimpleTriangle& s1, const SimpleTriangle& s2)
				{
					return s1.c[a] <= s2.c[a];
				});

			
				float length = n->box.bbox[1][a] - n->box.bbox[0][a];
				
				float bin_length = length / num_bin;

				vector<vec3> split_position;

				split_position.resize(num_bin - 1);

				for (int i = 1; i < num_bin; ++i)
					split_position[i - 1] = v_min + i * bin_length;

				vector<BBox> left_box;
				vector<BBox> right_box;

				vector<int> count_triangle_left(num_bin - 1, 0);
				vector<int> count_triangle_right(num_bin - 1, 0);


				left_box.resize(num_bin - 1);
				right_box.resize(num_bin - 1);
				
				//triangle loop outside split position loop
				for (int i = start; i < start + range; ++i)
				{
					for (int b = 0; b < num_bin - 1; ++b)
					{
						//if it on the left of split point, this box go to left box
						//we put the triangle loop outside to achive memory coherent
						//result it faster build
						if (simp[i].b.c()[a] < split_position[b][a])
						{
							left_box[b].expand(simp[i].b);
							++count_triangle_left[b];
						}
						else
						{
							right_box[b].expand(simp[i].b);
							++count_triangle_right[b];
						}
					}
				}			

				for (int b = 0; b < num_bin - 1; ++b)
				{
					float cost = ct + ci * (left_box[b].area() * count_triangle_left[b] + right_box[b].area() * count_triangle_right[b]) * inv_area;
				
					if (cost < best_cost)
					{
						best_axis = a;
						best_cost = cost;
						best_split = count_triangle_left[b];
					}
				}
			}
			if (best_cost == ci * range)
			{
				n->leaf = 1;
				n->range = range;
				n->axis = n->box.maxDim();
				int axis = n->axis;
				sort(simp.begin() + start, simp.begin() + start + range, [axis](const SimpleTriangle& s1, const SimpleTriangle& s2)
				{
					return s1.c[axis] < s2.c[axis];
				});
				return;
			}

			n->axis = best_axis;

			sort(simp.begin() + start, simp.begin() + start + range, [best_axis](const SimpleTriangle& s1, const SimpleTriangle& s2)
			{
				return s1.c[best_axis] < s2.c[best_axis];
			});

			sah(n->nodes[0], start, best_split, boxes, simp);

			sah(n->nodes[1], start + best_split, range - best_split, boxes, simp);
		}
	}


	void build_bvh(node*& n, const int& l)
	{
		leaf = l;
		int s = trs.size();

		vector<BBox> boxes;
		vector<SimpleTriangle> simp;
		boxes.resize(s);
		simp.resize(s);

		cout << "Num Triangle: " << s << "\n";
		
		//cout << "vertex size: " << v.size() << "\n";

		for (int i = 0; i < s; ++i)
		{
			int ind0 = trs[i].ind_v0;
			int ind1 = trs[i].ind_v1;
			int ind2 = trs[i].ind_v2;

			int v0 = indices[ind0].v;
			int v1 = indices[ind1].v;
			int v2 = indices[ind2].v;

			//int v0 = trs[i].v0;
			//int v1 = v0 + trs[i].v0_v1_offset;
			//int v2 = v0 + trs[i].v0_v2_offset;

			//cout << "vertex index: " << v0 << "\n";

			vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
			vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
			vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);

			/*p0 *= v_extend_i65535;
			p1 *= v_extend_i65535;
			p2 *= v_extend_i65535;

			p0 += v_min;
			p1 += v_min;
			p2 += v_min;
			*/

			//vec3 extend(v_max - v_min);

			//p0 *= vec3(65535.0f / extend.x, 65535.0f / extend.y, 65535.0f / extend.z)  ;
			//p1 *= vec3(65535.0f / extend.x, 65535.0f / extend.y, 65535.0f / extend.z);
			//p2 *= vec3(65535.0f / extend.x, 65535.0f / extend.y, 65535.0f / extend.z);

			p0 *= v_extend_i65535;//v_extend * i65535;
			p1 *= v_extend_i65535;//v_extend * i65535;
			p2 *= v_extend_i65535;//v_extend * i65535;

			//cout << v_min.x << " " << v_min.y << " " << v_min.z << "\n";

			p0 += v_min;
			p1 += v_min;
			p2 += v_min;

			//cout << p0.x << " " << p0.y << " " << p0.z << "\n";

			boxes[i].expand(p0);
			boxes[i].expand(p1);
			boxes[i].expand(p2);
			simp[i] = { boxes[i].c(), i };
		}

		sah(n, 0, s, boxes, simp);

		//sah_bin(n, 0, s, boxes, simp);

		vector<Triangle_16> new_trs;

		new_trs.resize(s);

		for (int i = 0; i < s; ++i)
		{
			new_trs[i] = trs[simp[i].i];
		}
		trs = new_trs;


		vector<SimpleTriangle>().swap(simp);


		vector<Triangle_16>().swap(new_trs);

		//Add Light Source
		if (use_area_light)
		{
			int li = 0;
			for (int i = 0; i < trs.size(); ++i)
			{
				int mtl = trs[i].mtl;

				if (mats[mtl].isLight)
				{
					int ind0 = trs[i].ind_v0;
					int ind1 = trs[i].ind_v1;
					int ind2 = trs[i].ind_v2;

					int v0 = indices[ind0].v;
					int v1 = indices[ind1].v;
					int v2 = indices[ind2].v;

					//int v0 = trs[i].v0;
					//int v1 = v0 + trs[i].v0_v1_offset;
					//int v2 = v0 + trs[i].v0_v2_offset;

					//vec3 p0 = v[v0];
					//vec3 p1 = v[v1];
					//vec3 p2 = v[v2];

					vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
					vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
					vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);

					p0 *= v_extend_i65535;
					p1 *= v_extend_i65535;
					p2 *= v_extend_i65535;

					p0 += v_min;
					p1 += v_min;
					p2 += v_min;

					Ls.add_light(p0, p1, p2);
					Ls.Ke.emplace_back(mats[mtl].Ke);

					light_map[i] = li++;

				}
			}

			if (Ls.Ke.size() == 0)
			{

				//cout << "Im here\n";
				//cout << "Add Light On Top\n";
				//vertex1
				vec3 v1(v_min.x, v_max.y, v_min.z);
				//vertex2
				vec3 v2(v_min.x, v_max.y, v_max.z);
				//vertex3
				vec3 v3(v_max.x, v_max.y, v_max.z);
				//vertex4
				vec3 v4(v_max.x, v_max.y, v_min.z);

				//convert to scale coodinate

				vec3 box_extend = v_max - v_min;

				vec3 inv_dimension(1.0f / box_extend.x, 1.0f / box_extend.y, 1.0f / box_extend.z);

				vec3 scale = vec3(65535.0f) * inv_dimension;

				//converte_vertex1
				
				vec3 v1_scale = (v1 - v_min) * scale;
				vec3 v2_scale = (v2 - v_min) * scale;
				vec3 v3_scale = (v3 - v_min) * scale;
				vec3 v4_scale = (v4 - v_min) * scale;



				Ls.add_light(v1, v3, v2);
				Ls.add_light(v1, v4, v3);

				//failed
				//Ls.add_light(v1, v2, v4);
				//Ls.add_light(v2, v3, v4);


				vec3 ke1(4.0f);
				vec3 ke2(4.0f);

				//Ls.Ke.resize(2);

				Ls.Ke.emplace_back(ke1);
				Ls.Ke.emplace_back(ke2);

				//Ls.Ke[0] = ke1;
				//Ls.Ke[1] = ke2;

				//still have directional light
				//but when you dont want to use enviroment
				//this function will use area light instead
				if (!use_eviroment)
				{
					//int v_size = v.size();

					//cout << v_size << "\n";

					int index_size = indices.size();

					//v.emplace_back(v1);
					//v.emplace_back(v2);
					//v.emplace_back(v3);
					//v.emplace_back(v4);

					v.emplace_back(v1_scale.x); v.emplace_back(v1_scale.y); v.emplace_back(v1_scale.z);
					v.emplace_back(v2_scale.x); v.emplace_back(v2_scale.y); v.emplace_back(v2_scale.z);
					v.emplace_back(v3_scale.x); v.emplace_back(v3_scale.y); v.emplace_back(v3_scale.z);
					v.emplace_back(v4_scale.x); v.emplace_back(v4_scale.y); v.emplace_back(v4_scale.z);



					int mtl = mats.size();

					Material mat;
					mat.isLight = true;
					mat.Ke = vec3(4.0f);

					mats.emplace_back(mat);

					//Triangle tr1(v1, v3, v2, -1, -1, -1, -1, -1, -1, mtl);
					//Triangle tr2(v1, v4, v3, -1, -1, -1, -1, -1, -1, mtl);

					int v1 = index_size; //v_size;
					int v2 = index_size + 1;//v_size + 1;
					int v3 = index_size + 2;//v_size + 2;
					int v4 = index_size + 3;//v_size + 3;


					Triangle_16 tr1(v1, v3, v2, mtl);
					Triangle_16 tr2(v1, v4, v3, mtl);

					//Triangle_24 tr1(v1, v3, v2, -1, -1, -1, -1, -1, -1, mtl);
					//Triangle_24 tr2(v1, v4, v3, -1, -1, -1, -1, -1, -1, mtl);

					trs.emplace_back(tr1);
					trs.emplace_back(tr2);
				}
			}
			Ls.norm();

		}

	}

	vec3 sample_triangle(int& triangle_ind)
	{
		int ind0 = trs[triangle_ind].ind_v0;
		int ind1 = trs[triangle_ind].ind_v1;
		int ind2 = trs[triangle_ind].ind_v2;

		int v0 = 3 * indices[ind0].v;
		int v1 = 3 * indices[ind1].v;
		int v2 = 3 * indices[ind2].v;

		vec3 p0(v[v0], v[v0 + 1], v[v0 + 2]);
		vec3 p1(v[v1], v[v1 + 1], v[v1 + 2]);
		vec3 p2(v[v2], v[v2 + 1], v[v2 + 2]);

		p0 *= v_extend_i65535;
		p1 *= v_extend_i65535;
		p2 *= v_extend_i65535;

		p0 += v_min;
		p1 += v_min;
		p2 += v_min;

		float su0 = sqrt14(randf());
		float b0 = 1.0f - su0;
		float b1 = randf() * su0;

		return b0 * p0 + b1 * p1 + (1.0f - b0 - b1) * p2;
	}

	void compute_area()
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < round_light[i].size(); ++j)
			{
				int ind = round_light[i][j];

				int ind0 = trs[ind].ind_v0;
				int ind1 = trs[ind].ind_v1;
				int ind2 = trs[ind].ind_v2;

				int v0 = 3 * indices[ind0].v;
				int v1 = 3 * indices[ind1].v;
				int v2 = 3 * indices[ind2].v;

				vec3 p0(v[v0], v[v0 + 1], v[v0 + 2]);
				vec3 p1(v[v1], v[v1 + 1], v[v1 + 2]);
				vec3 p2(v[v2], v[v2 + 1], v[v2 + 2]);

				p0 *= v_extend_i65535;
				p1 *= v_extend_i65535;
				p2 *= v_extend_i65535;

				///p0 += v_min;
				///p1 += v_min;
				///p2 += v_min;

				vec3 cr((p0 - p2).cross(p1 - p2));

				//cout << "a: " << a.x << " " << a.y << " " << a.z << "\n";
				//cout << "b: " << b.x << " " << b.y << " " << b.z << "\n";
				//cout << "c: " << c.x << " " << c.y << " " << c.z << "\n";
				float area = sqrt14(cr.dot(cr)) * 0.5f;

				round_light_area[i].push_back(area);
			}
		}
	}

	bool __fastcall all_hit(node*& n, Ray& r, HitRecord& rec)
	{
		node* stack[26];

		stack[0] = n;
		float mint = 10000000.0f;
		int si = 0;
		while (si >= 0)
		{
			auto top(stack[si]);
			si--;

			float tl;

			//if (top->box.hit_axis_tl(r, tl) && tl < mint)
			if (top->box.hit_axis_tl_mint(r, tl, mint))
			{
				if (!top->leaf)
				{
					int sort_axis = r.sign[top->axis];

					stack[si + 1] = top->nodes[sort_axis];
					stack[si + 2] = top->nodes[1 - sort_axis];
					si += 2;
				}
				else
				{
					auto start = top->start, end = start + top->range;

					for (auto i = start; i < end; ++i)
					{
						int ind0 = trs[i].ind_v0;
						int ind1 = trs[i].ind_v1;
						int ind2 = trs[i].ind_v2;

						int v0 = indices[ind0].v;
						int v1 = indices[ind1].v;
						int v2 = indices[ind2].v;

						//vec3 p0 = v[v0];
						//vec3 p1 = v[v1];
						//vec3 p2 = v[v2];

						vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
						vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
						vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);

						p0 *= v_extend_i65535;
						p1 *= v_extend_i65535;
						p2 *= v_extend_i65535;

						p0 += v_min;
						p1 += v_min;
						p2 += v_min;

						vec3 e1(p1 - p0);
						vec3 e2(p2 - p0);

						hit_triangle(r, p0, e1, e2, i, &mint, rec);
						//trs[i].optimize_2(r, rec, &mint, i);
					}

				}
			}
			else
				continue;

		}
		if (mint == 10000000.0f)
		{
			return false;
		}

		return true;
	}

	vec2 __fastcall intp(const vec2& a, const vec2& b, const vec2& c, const float& u, const float& v)
	{
		return a * (1.0f - u - v) + b * u + c * v;
	}

	vec3 __fastcall intp(const vec3& a, const vec3& b, const vec3& c, const float& u, const float& v)
	{
		return a * (1.0f - u - v) + b * u + c * v;
	}

	bool __fastcall hit_color_alpha_2(node*& n, Ray& r, HitRecord& rec, vec3* color, float& transparent)
	{
		node* stack[32];

		//stack[0] = n;
		int a = r.sign[n->axis];

		stack[0] = n->nodes[a];
		stack[1] = n->nodes[1 - a];
		float mint = 10000000.0f;
		int si = 1;

		//ray_t new_ray(r.o, r.invd);

		while (si >= 0)
		{
			auto top(stack[si]);
			--si;

			float tl;


			if (top->box.hit_axis_tl(r, tl) && tl < mint)
				//if (top->box.hit_axis_tl_mint(r, tl, mint))// && tl < mint)
			{

				int sort_axis = r.sign[top->axis];

				stack[si + 1] = top->nodes[sort_axis];
				stack[si + 2] = top->nodes[1 - sort_axis];

				si += 2;

				if (top->leaf)
				{
					si -= 2;
					auto start = top->start, end = start + top->range;

					for (auto i = start; i < end; ++i)
					{
						int ind0 = trs[i].ind_v0;
						int ind1 = trs[i].ind_v1;
						int ind2 = trs[i].ind_v2;

						int v0 = indices[ind0].v;
						int v1 = indices[ind1].v;
						int v2 = indices[ind2].v;

						//vec3 p0(v[v0]);
						//vec3 p1(v[v1]);
						//vec3 p2(v[v2]);

						vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
						vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
						vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);

						p0 *= v_extend_i65535;
						p1 *= v_extend_i65535;
						p2 *= v_extend_i65535;

						p0 += v_min;
						p1 += v_min;
						p2 += v_min;

						vec3 e1(p1 - p0);
						vec3 e2(p2 - p0);

						hit_triangle(r, p0, e1, e2, i, &mint, rec);
						//trs[i].optimize_2(r, rec, &mint, i);
					}
				}
			}
			else
				continue;

		}
		if (mint < 10000000.0f)
		{
			rec.t = mint;
			uint32_t ind = rec.ti, mtl_ind = trs[ind].mtl;

			int ti = rec.ti;

			int ind0 = trs[ti].ind_v0;
			int ind1 = trs[ti].ind_v1;
			int ind2 = trs[ti].ind_v2;

			int n0 = indices[ind0].vn;
			int n1 = indices[ind1].vn;
			int n2 = indices[ind2].vn;

			//vec3 rec_normal = intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v);

			//rec.n = (vn[vn0] + vn[vn1] + vn[vn2]) / 3.0f;//intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v);

			//rec.n = rec_normal;

			//rec.n = (vn[vn0] + vn[vn1] + vn[vn2]) / 3.0f;//intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v);

			//if (vn0 < 0)
			//	cout << "negative normal index\n";

			//cout << vn0 << "\n";

			//if (!trs[ti].compute_normal)
			//	cout << "not compute normal\n";

			//rec.n = trs[ti].vn[0] > 0 ? intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v).norm() : trs[ti].n.norm();

			if (n0 >= 0)
			{
				int n0_0 = 3 * n0;
				int n0_1 = 3 * n0 + 1;
				int n0_2 = 3 * n0 + 2;

				int n1_0 = 3 * n1;
				int n1_1 = 3 * n1 + 1;
				int n1_2 = 3 * n1 + 2;

				int n2_0 = 3 * n2;
				int n2_1 = 3 * n2 + 1;
				int n2_2 = 3 * n2 + 2;

				vec3 vn0(float(vn[n0_0]) * 0.0001f, float(vn[n0_1]) * 0.0001f, float(vn[n0_2]) * 0.0001f);
				vec3 vn1(float(vn[n1_0]) * 0.0001f, float(vn[n1_1]) * 0.0001f, float(vn[n1_2]) * 0.0001f);
				vec3 vn2(float(vn[n2_0]) * 0.0001f, float(vn[n2_1]) * 0.0001f, float(vn[n2_2]) * 0.0001f);

				rec.n = intp(vn0, vn1, vn2, rec.u, rec.v).norm();
				//rec.n = intp(vn[n0], vn[n1], vn[n2], rec.u, rec.v).norm();
			}
			else
			{
				int ind0 = trs[ti].ind_v0;
				int ind1 = trs[ti].ind_v1;
				int ind2 = trs[ti].ind_v2;

				int v0 = indices[ind0].v;
				int v1 = indices[ind1].v;
				int v2 = indices[ind2].v;

				vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
				vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
				vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);

				p0 *= v_extend_i65535;
				p1 *= v_extend_i65535;
				p2 *= v_extend_i65535;

				p0 += v_min;
				p1 += v_min;
				p2 += v_min;

				vec3 e1(p1 - p0);
				vec3 e2(p2 - p0);

				rec.n = (e1.cross(e2)).norm();
			}

			/*if (trs[ti].compute_normal)
			{
			int vn0 = trs[ti].vn[0];
			int vn1 = trs[ti].vn[1];
			int vn2 = trs[ti].vn[2];

			rec.n = intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v);
			}
			else
			rec.n = trs[ti].n;
			*/
			//rec.n = trs[ti].n;

			if (mats[mtl_ind].use_texture)
			{
				int ind0 = trs[ti].ind_v0;
				int ind1 = trs[ti].ind_v1;
				int ind2 = trs[ti].ind_v2;

				int t0 = indices[ind0].vt;
				int t1 = indices[ind1].vt;
				int t2 = indices[ind2].vt;

				int t0_0 = 2 * t0;
				int t0_1 = 2 * t0 + 1;

				int t1_0 = 2 * t1;
				int t1_1 = 2 * t1 + 1;

				int t2_0 = 2 * t2;
				int t2_1 = 2 * t2 + 1;

				vec2 vt0 = vec2(vt[t0_0] * 0.0001f, vt[t0_1] * 0.0001f);
				vec2 vt1 = vec2(vt[t1_0] * 0.0001f, vt[t1_1] * 0.0001f);
				vec2 vt2 = vec2(vt[t2_0] * 0.0001f, vt[t2_1] * 0.0001f);

				vec2 t(intp(vt0, vt1, vt2, rec.u, rec.v));

				int tex_ind = mats[mtl_ind].TexInd;

				*color = texs[tex_ind].ev(t) * i255;

				//opaque cua transparent material se de xu ly rieng
				transparent = mats[mtl_ind].Tr;//mats[mtl_ind].opaque;

											   /*if (!mats[mtl_ind].alpha_material)
											   {
											   alpha = 1.0f;
											   //return true;
											   }
											   else
											   {
											   alpha = texs[tex_ind].ev_alpha(t);

											   //return true;
											   }*/

				if (mats[mtl_ind].use_alpha_texture)
				{
					//float alpha_texture = texs[tex_ind].ev_alpha(t);
					//transparent = (1.0f - alpha_texture);

					transparent = texs[tex_ind].ev_alpha(t);

					//transparent = transparent  
				}
				//if (transparent < 1.0f)
				//	transparent = 0.0f;
				/*
				if (mats[mtl_ind].use_Bump)
				{
				vec3 original_rec_n = rec.n;


				vec3 bump = texs[tex_ind].ev_bump(t);

				bump = (2.0f * bump - vec3(1.0f));

				//TBN

				int ind0 = trs[ti].ind_v0;
				int ind1 = trs[ti].ind_v1;
				int ind2 = trs[ti].ind_v2;

				int v0 = indices[ind0].v;
				int v1 = indices[ind1].v;
				int v2 = indices[ind2].v;

				vec3 p0 = v[v0];
				vec3 p1 = v[v1];
				vec3 p2 = v[v2];

				vec3 e1 = p1 - p0;
				vec3 e2 = p2 - p0;

				vec2 delta_uv1 = vt[t1] - vt[t0];
				vec2 delta_uv2 = vt[t2] - vt[t0];

				float f = 1.0f / (delta_uv1.x * delta_uv2.y - delta_uv1.y * delta_uv2.x);

				vec3 T, B, N;

				N = rec.n;

				T = f * (e1 * delta_uv2.y - e2 * delta_uv1.y);
				B = f * (e2 * delta_uv1.x - e1 * delta_uv2.x);

				T = B.cross(N).norm();
				B = T.cross(N).norm();

				rec.n = (rec.n * bump.z + T * bump.x + B * bump.y).norm();

				}
				*/

				return true;

			}

			*color = mats[mtl_ind].Kd;


			return true;
		}
		return false;

	}

	bool __fastcall hit_color_alpha(node*& n, Ray& r, HitRecord& rec, vec3* color, float& transparent)
	{
		node* stack[32];

		//stack[0] = n;
		int a = r.sign[n->axis];

		stack[0] = n->nodes[a];
		stack[1] = n->nodes[1 - a];
		float mint = 10000000.0f;
		int si = 1;
		

		//float mint = 10000000.0f;
		//int si = 0;
		//ray_t new_ray(r.o, r.invd);

		while (si >= 0)
		{
			auto top(stack[si]);
			--si;

			float tl;

			

			if (top->box.hit_axis_tl(r, tl) && tl < mint)
				//if (top->box.hit_axis_tl_mint(r, tl, mint))// && tl < mint)
			{
				int sort_axis = r.sign[top->axis];

				stack[si + 1] = top->nodes[sort_axis];
				stack[si + 2] = top->nodes[1 - sort_axis];

				si += 2;
				if (top->leaf)
				{
					si -= 2;
					auto start = top->start, end = start + top->range;

					for (auto i = start; i < end; ++i)
					{
						int ind0 = trs[i].ind_v0;
						int ind1 = trs[i].ind_v1;
						int ind2 = trs[i].ind_v2;

						/*int v0 = indices[ind0].v;
						int v1 = indices[ind1].v;
						int v2 = indices[ind2].v;

						//vec3 p0(v[v0]);
						//vec3 p1(v[v1]);
						//vec3 p2(v[v2]);

						vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
						vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
						vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);*/

						int v0 = 3 * indices[ind0].v;
						int v1 = 3 * indices[ind1].v;
						int v2 = 3 * indices[ind2].v;

						vec3 p0(v[v0], v[v0 + 1], v[v0 + 2]);
						vec3 p1(v[v1], v[v1 + 1], v[v1 + 2]);
						vec3 p2(v[v2], v[v2 + 1], v[v2 + 2]);

						p0 *= v_extend_i65535;
						p1 *= v_extend_i65535;
						p2 *= v_extend_i65535;

						p0 += v_min;
						p1 += v_min;
						p2 += v_min;

						vec3 e1(p1 - p0);
						vec3 e2(p2 - p0);

						hit_triangle(r, p0, e1, e2, i, &mint, rec);
						//trs[i].optimize_2(r, rec, &mint, i);
					}
					//continue;
				}
				//else
				//{
					
				//}
			}
			else
			//{
				
				continue;
			//}
		}
		if (mint < 10000000.0f)
		{
			rec.t = mint;
			uint32_t ind = rec.ti, mtl_ind = trs[ind].mtl;

			int ti = rec.ti;

			/*
			int ind0 = trs[ti].ind_v0;
			int ind1 = trs[ti].ind_v1;
			int ind2 = trs[ti].ind_v2;

			int n0 = indices[ind0].vn;
			int n1 = indices[ind1].vn;
			int n2 = indices[ind2].vn;
			*/

			//vec3 rec_normal = intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v);

			//rec.n = (vn[vn0] + vn[vn1] + vn[vn2]) / 3.0f;//intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v);

			//rec.n = rec_normal;

			//rec.n = (vn[vn0] + vn[vn1] + vn[vn2]) / 3.0f;//intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v);

			//if (vn0 < 0)
			//	cout << "negative normal index\n";

			//cout << vn0 << "\n";

			//if (!trs[ti].compute_normal)
			//	cout << "not compute normal\n";

			//rec.n = trs[ti].vn[0] > 0 ? intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v).norm() : trs[ti].n.norm();

			int ind0 = trs[ti].ind_v0;
			int n0 = indices[ind0].vn;

			if (n0 >= 0)
			{
				/*int n0_0 = 3 * n0;
				int n0_1 = n0_0 + 1;//3 * n0 + 1;
				int n0_2 = n0_0 + 2;//3 * n0 + 2;

				int n1_0 = 3 * n1;
				int n1_1 = n1_0 + 1;//3 * n1 + 1;
				int n1_2 = n1_0 + 2;//3 * n1 + 2;

				int n2_0 = 3 * n2;
				int n2_1 = n2_0 + 1;//3 * n2 + 1;
				int n2_2 = n2_0 + 2;//3 * n2 + 2;

				vec3 vn0(float(vn[n0_0]) * 0.0001f, float(vn[n0_1]) * 0.0001f, float(vn[n0_2]) * 0.0001f);
				vec3 vn1(float(vn[n1_0]) * 0.0001f, float(vn[n1_1]) * 0.0001f, float(vn[n1_2]) * 0.0001f);
				vec3 vn2(float(vn[n2_0]) * 0.0001f, float(vn[n2_1]) * 0.0001f, float(vn[n2_2]) * 0.0001f);
				*/

				int ind0 = trs[ti].ind_v0;
				int ind1 = trs[ti].ind_v1;
				int ind2 = trs[ti].ind_v2;

				int n0 = 3 * indices[ind0].vn;
				int n1 = 3 * indices[ind1].vn;
				int n2 = 3 * indices[ind2].vn;

				vec3 vn0(float(vn[n0]) * 0.0001f, float(vn[n0 + 1]) * 0.0001f, float(vn[n0 + 2]) * 0.0001f);
				vec3 vn1(float(vn[n1]) * 0.0001f, float(vn[n1 + 1]) * 0.0001f, float(vn[n1 + 2]) * 0.0001f);
				vec3 vn2(float(vn[n2]) * 0.0001f, float(vn[n2 + 1]) * 0.0001f, float(vn[n2 + 2]) * 0.0001f);

				rec.n = intp(vn0, vn1, vn2, rec.u, rec.v).norm();
				//rec.n = intp(vn[n0], vn[n1], vn[n2], rec.u, rec.v).norm();
			}
			else
			{
				int ind0 = trs[ti].ind_v0;
				int ind1 = trs[ti].ind_v1;
				int ind2 = trs[ti].ind_v2;

				int v0 = 3 * indices[ind0].v;
				int v1 = 3 * indices[ind1].v;
				int v2 = 3 * indices[ind2].v;

				vec3 p0(v[v0], v[v0 + 1], v[v0 + 2]);
				vec3 p1(v[v1], v[v1 + 1], v[v1 + 2]);
				vec3 p2(v[v2], v[v2 + 1], v[v2 + 2]);

				p0 *= v_extend_i65535;
				p1 *= v_extend_i65535;
				p2 *= v_extend_i65535;

				p0 += v_min;
				p1 += v_min;
				p2 += v_min;

				vec3 e1(p1 - p0);
				vec3 e2(p2 - p0);

				rec.n = (e1.cross(e2)).norm();
			}

			/*if (trs[ti].compute_normal)
			{
			int vn0 = trs[ti].vn[0];
			int vn1 = trs[ti].vn[1];
			int vn2 = trs[ti].vn[2];

			rec.n = intp(vn[vn0], vn[vn1], vn[vn2], rec.u, rec.v);
			}
			else
			rec.n = trs[ti].n;
			*/
			//rec.n = trs[ti].n;
			
			transparent = mats[mtl_ind].Tr;

			if (mats[mtl_ind].use_texture)
			{
				int ind0 = trs[ti].ind_v0;
				int ind1 = trs[ti].ind_v1;
				int ind2 = trs[ti].ind_v2;

				int t0 = indices[ind0].vt;
				int t1 = indices[ind1].vt;
				int t2 = indices[ind2].vt;

				int t0_0 = t0 << 1;//2 * t0;
				int t0_1 = t0_0 + 1;//2 * t0 + 1;

				int t1_0 = t1 << 1;//2 * t1;
				int t1_1 = t1_0 + 1;//2 * t1 + 1;

				int t2_0 = t2 << 1;//2 * t2;
				int t2_1 = t2_0 + 1;//2 * t2 + 1;

				vec2 vt0 = vec2(vt[t0_0] * 0.0001f, vt[t0_1] * 0.0001f);
				vec2 vt1 = vec2(vt[t1_0] * 0.0001f, vt[t1_1] * 0.0001f);
				vec2 vt2 = vec2(vt[t2_0] * 0.0001f, vt[t2_1] * 0.0001f);

				vec2 t(intp(vt0, vt1, vt2, rec.u, rec.v));

				int tex_ind = mats[mtl_ind].TexInd;

				*color = texs[tex_ind].ev(t) * i255;

				//opaque cua transparent material se de xu ly rieng
				//transparent = mats[mtl_ind].Tr;//mats[mtl_ind].opaque;

							 /*if (!mats[mtl_ind].alpha_material)
							 {
							 alpha = 1.0f;
							 //return true;
							 }
							 else
							 {
							 alpha = texs[tex_ind].ev_alpha(t);

							 //return true;
							 }*/

				if (mats[mtl_ind].use_alpha_texture)
				{
					//float alpha_texture = texs[tex_ind].ev_alpha(t);
					//transparent = (1.0f - alpha_texture);

					transparent = texs[tex_ind].ev_alpha(t);

					//transparent = transparent  
				}
				//if (transparent < 1.0f)
				//	transparent = 0.0f;
				
				if (mats[mtl_ind].use_Bump)
				{
					vec3 original_rec_n = rec.n;


					vec3 bump = texs[tex_ind].ev_bump(t);

					bump = (2.0f * bump - vec3(1.0f));

					//TBN

					int ind0 = trs[ti].ind_v0;
					int ind1 = trs[ti].ind_v1;
					int ind2 = trs[ti].ind_v2;

					int v0 = indices[ind0].v;
					int v1 = indices[ind1].v;
					int v2 = indices[ind2].v;

					vec3 p0 = v[v0];
					vec3 p1 = v[v1];
					vec3 p2 = v[v2];

					vec3 e1 = p1 - p0;
					vec3 e2 = p2 - p0;

					vec2 delta_uv1 = vt1 - vt0;//vt[t1] - vt[t0];
					vec2 delta_uv2 = vt2 - vt0;//vt[t2] - vt[t0];

					float f = 1.0f / (delta_uv1.x * delta_uv2.y - delta_uv1.y * delta_uv2.x);

					vec3 T, B, N;

					N = rec.n;

					T = f * (e1 * delta_uv2.y - e2 * delta_uv1.y);
					B = f * (e2 * delta_uv1.x - e1 * delta_uv2.x);

					T = B.cross(N).norm();
					B = T.cross(N).norm();

					rec.n = (rec.n * bump.z + T * bump.x + B * bump.y).norm();

				}
				

				return true;

			}

			*color = mats[mtl_ind].Kd;


			return true;
		}
		return false;

	}

	//changeed
	//bool __fastcall hit_anything_range(node*& n, Ray& r, const float& d)
	bool __fastcall hit_anything_range(node*& n, Ray& r)
	{
		node* stack[32];

		stack[0] = n->nodes[second_shadow_node];
		stack[1] = n->nodes[first_shadow_node];

		int si = 1;

		//stack[0] = n;

		while (si >= 0)
		{
			auto top(stack[si--]);

			float th;

			if (top->box.hit_shadow_th(r, th))// && th >= 0.0f)
			{
				if (!top->leaf)
				{

					stack[si + 1] = top->nodes[second_shadow_node];
					stack[si + 2] = top->nodes[first_shadow_node];


					si += 2;
				}
				else
				{
					auto start = top->start, end = start + top->range;

					for (auto i = start; i < end; ++i)
						//for (auto i = end - 1; i >= start; ++i)
					{
						int mtl = trs[i].mtl;
						//if (trs[i].hit_anything(r, d))// && !mats[mtl].isLight)
						//if(trs[i].hit_anything(r, d) && !mats[mtl].isLight)

						//if (trs[i].hit_anything(r) && !mats[mtl].isLight)
						//	return true;

						int ind0 = trs[i].ind_v0;
						int ind1 = trs[i].ind_v1;
						int ind2 = trs[i].ind_v2;

						int v0 = 3 * indices[ind0].v;
						int v1 = 3 * indices[ind1].v;
						int v2 = 3 * indices[ind2].v;

						//vec3 p0 = v[v0];
						//vec3 p1 = v[v1];
						//vec3 p2 = v[v2];

						/*vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
						vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
						vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);*/

						vec3 p0(v[v0], v[v0 + 1], v[v0 + 2]);
						vec3 p1(v[v1], v[v1 + 1], v[v1 + 2]);
						vec3 p2(v[v2], v[v2 + 1], v[v2 + 2]);

						p0 *= v_extend_i65535;
						p1 *= v_extend_i65535;
						p2 *= v_extend_i65535;

						p0 += v_min;
						p1 += v_min;
						p2 += v_min;

						vec3 e1(p1 - p0);
						vec3 e2(p2 - p0);

						if (hit_anything(r, p0, e1, e2) && !mats[mtl].isLight)
						{
							//cout << "Light Block\n";
							return true;
						}
					}
				}
			}
			else
				continue;
		}
		return false;
	}


	bool __fastcall hit_anything_alpha(node*& n, Ray& r)
	{
		node* stack[32];

		stack[0] = n->nodes[second_shadow_node];
		stack[1] = n->nodes[first_shadow_node];

		int si = 1;

		//stack[0] = n;

		while (si >= 0)
		{
			auto top(stack[si--]);

			float th;

			if (top->box.hit_shadow_th(r, th))// && th >= 0.0f)
			{
				if (!top->leaf)
				{

					stack[si + 1] = top->nodes[second_shadow_node];
					stack[si + 2] = top->nodes[first_shadow_node];


					si += 2;
				}
				else
				{
					auto start = top->start, end = start + top->range;

					for (auto i = start; i < end; ++i)
						//for (auto i = end - 1; i >= start; ++i)
					{
						int mtl = trs[i].mtl;
						//if (trs[i].hit_anything(r, d))// && !mats[mtl].isLight)
						//if(trs[i].hit_anything(r, d) && !mats[mtl].isLight)

						HitRecord rec;
						//if (trs[i].hit_anything_alpha(r, rec) && !mats[mtl].isLight)

						int ind0 = trs[i].ind_v0;
						int ind1 = trs[i].ind_v1;
						int ind2 = trs[i].ind_v2;

						int v0 = indices[ind0].v;
						int v1 = indices[ind1].v;
						int v2 = indices[ind2].v;

						//vec3 p0 = v[v0];
						//vec3 p1 = v[v1];
						//vec3 p2 = v[v2];

						vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
						vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
						vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);

						p0 *= v_extend_i65535;
						p1 *= v_extend_i65535;
						p2 *= v_extend_i65535;

						p0 += v_min;
						p1 += v_min;
						p2 += v_min;

						vec3 e1(p1 - p0);
						vec3 e2(p2 - p0);

						if (hit_triangle_alpha(r, p0, e1, e2, rec) && !mats[mtl].isLight)
						{
							if (!mats[mtl].alpha_material)
								return true;
							else
							{
								//this retrun false cause all our pain
								//this is always false and the most occluded part is litted light hell fire
								//return false;
								float alpha = mats[mtl].Tr;//mats[mtl].opaque;

								/*if (mats[mtl].use_alpha_texture)
								{
									int ind0 = trs[i].ind_v0;
									int ind1 = trs[i].ind_v1;
									int ind2 = trs[i].ind_v2;

									int t0 = indices[ind0].vt;
									int t1 = indices[ind1].vt;
									int t2 = indices[ind2].vt;

									int t0_0 = 2 * t0;
									int t0_1 = 2 * t0 + 1;

									int t1_0 = 2 * t1;
									int t1_1 = 2 * t1 + 1;

									int t2_0 = 2 * t2;
									int t2_1 = 2 * t2 + 1;

									vec2 vt0 = vec2(vt[t0_0] * 0.0001f, vt[t0_1] * 0.0001f);
									vec2 vt1 = vec2(vt[t1_0] * 0.0001f, vt[t1_1] * 0.0001f);
									vec2 vt2 = vec2(vt[t2_0] * 0.0001f, vt[t2_1] * 0.0001f);

									vec2 t(intp(vt0, vt1, vt2, rec.u, rec.v));

									int tex_ind = mats[mtl].TexInd;
									float alpha_texture = texs[tex_ind].ev_alpha(t);

									alpha = alpha_texture;//(1.0f - alpha_texture);
									if (alpha > 0.0f)
										alpha = 1.0f;
								}*/

								if (randf() > alpha)
									return true;
								else
									continue;
							}
						}
					}
				}
			}
			else
				continue;
		}
		return false;
	}

	bool __fastcall hit_anything_range_alpha(node*& n, Ray& r, float& d)
	{
		node* stack[32];

		

		//stack[0] = n->nodes[second_shadow_node];
		//stack[1] = n->nodes[first_shadow_node];

		int si = 0;

		stack[0] = n;

		while (si >= 0)
		{
			auto top(stack[si--]);

			float th;

			int sort_axis = r.sign[top->axis];

			if (top->box.hit_shadow_th(r, th))// && th >= 0.0f)
			{
				if (!top->leaf)
				{

					//stack[si + 1] = top->nodes[second_shadow_node];
					//stack[si + 2] = top->nodes[first_shadow_node];
					int sort_axis = r.sign[top->axis];

					stack[si + 1] = top->nodes[sort_axis];
					stack[si + 2] = top->nodes[1 - sort_axis];

					si += 2;
				}
				else
				{
					auto start = top->start, end = start + top->range;

					for (auto i = start; i < end; ++i)
					{
						int mtl = trs[i].mtl;

						HitRecord rec;
						//if (trs[i].hit_anything_alpha(r, rec) && rec.t <= d)

						int ind0 = trs[i].ind_v0;
						int ind1 = trs[i].ind_v1;
						int ind2 = trs[i].ind_v2;

						int v0 = indices[ind0].v;
						int v1 = indices[ind1].v;
						int v2 = indices[ind2].v;

						//vec3 p0 = v[v0];
						//vec3 p1 = v[v1];
						//vec3 p2 = v[v2];

						vec3 p0(v[3 * v0], v[3 * v0 + 1], v[3 * v0 + 2]);
						vec3 p1(v[3 * v1], v[3 * v1 + 1], v[3 * v1 + 2]);
						vec3 p2(v[3 * v2], v[3 * v2 + 1], v[3 * v2 + 2]);

						p0 *= v_extend_i65535;
						p1 *= v_extend_i65535;
						p2 *= v_extend_i65535;

						p0 += v_min;
						p1 += v_min;
						p2 += v_min;

						vec3 e1(p1 - p0);
						vec3 e2(p2 - p0);

						//if (hit_triangle_alpha(r, p0, e1, e2, rec))
						if (hit_triangle_alpha(r, p0, e1, e2, rec) && 1.00024f * rec.t < d)
						{
							if (!mats[mtl].alpha_material)
								return true;
							else
							{
								float alpha = mats[mtl].Tr;

								/*if (mats[mtl].use_alpha_texture)
								{
									int ind0 = trs[i].ind_v0;
									int ind1 = trs[i].ind_v1;
									int ind2 = trs[i].ind_v2;

									int t0 = indices[ind0].vt;
									int t1 = indices[ind1].vt;
									int t2 = indices[ind2].vt;

									int t0_0 = 2 * t0;
									int t0_1 = 2 * t0 + 1;

									int t1_0 = 2 * t1;
									int t1_1 = 2 * t1 + 1;

									int t2_0 = 2 * t2;
									int t2_1 = 2 * t2 + 1;

									vec2 vt0 = vec2(vt[t0_0] * 0.0001f, vt[t0_1] * 0.0001f);
									vec2 vt1 = vec2(vt[t1_0] * 0.0001f, vt[t1_1] * 0.0001f);
									vec2 vt2 = vec2(vt[t2_0] * 0.0001f, vt[t2_1] * 0.0001f);

									vec2 t(intp(vt0, vt1, vt2, rec.u, rec.v));

									int tex_ind = mats[mtl].TexInd;
									float alpha_texture = texs[tex_ind].ev_alpha(t);

									alpha = (1.0f - alpha_texture);
								}*/

								if (randf() > alpha)
									return true;
								else
									continue;
							}
						}
					}
				}
			}
			else
				continue;
		}
		return false;
	}

	bool __fastcall hit_anything_range_2(node*& n, Ray& r, float& d)
	{
		node* stack[32];

		//stack[0] = n->nodes[second_shadow_node];
		//stack[1] = n->nodes[first_shadow_node];

		int si = 0;

		stack[0] = n;

		while (si >= 0)
		{
			auto top(stack[si--]);

			float th;

			if (top->box.hit_shadow_th(r, th))// && th >= 0.0f)
			{
				if (!top->leaf)
				{

					//stack[si + 1] = top->nodes[second_shadow_node];
					//stack[si + 2] = top->nodes[first_shadow_node];

					int sort_axis = r.sign[top->axis];

					stack[si + 1] = top->nodes[sort_axis];
					stack[si + 2] = top->nodes[1 - sort_axis];

					si += 2;
				}
				else
				{
					auto start = top->start, end = start + top->range;

					for (auto i = start; i < end; ++i)
					{
						int mtl = trs[i].mtl;

						HitRecord rec;
						//if (trs[i].hit_anything_alpha(r, rec) && rec.t <= d)

						int ind0 = trs[i].ind_v0;
						int ind1 = trs[i].ind_v1;
						int ind2 = trs[i].ind_v2;

						int v0 = 3 * indices[ind0].v;
						int v1 = 3 * indices[ind1].v;
						int v2 = 3 * indices[ind2].v;

						//vec3 p0 = v[v0];
						//vec3 p1 = v[v1];
						//vec3 p2 = v[v2];

						vec3 p0(v[v0], v[v0 + 1], v[v0 + 2]);
						vec3 p1(v[v1], v[v1 + 1], v[v1 + 2]);
						vec3 p2(v[v2], v[v2 + 1], v[v2 + 2]);

						p0 *= v_extend_i65535;
						p1 *= v_extend_i65535;
						p2 *= v_extend_i65535;

						p0 += v_min;
						p1 += v_min;
						p2 += v_min;

						vec3 e1(p1 - p0);
						vec3 e2(p2 - p0);

						float min_t;
						unsigned int ti;

						//if (hit_triangle(r, p0, e1, e2, ti, &min_t, rec) && rec.t > eps && rec.t < d)
						if(hit_triangle_alpha(r, p0, e1, e2, rec) && rec.t > eps && 1.00000024f * rec.t <= d)
							return true;
					}
				}
			}
			else
				continue;
		}
		return false;
	}

	float __fastcall mis2(const float& pdf1, const float& pdf2)
	{
		float pdf1_2 = pdf1 * pdf1;
		return pdf1_2 / (pdf1_2 + pdf2 * pdf2);
	}


	vec3 __fastcall sibenik_tracing_artificial_light_normal_Le(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		vec3 Le(30.0f);


		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_alpha_pointer(n, new_ray, rec, &color, alpha))
			float alpha = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, alpha))
			{
				
				if (randf() < alpha)
				{
					int mtl = trs[rec.ti].mtl;
					T *= mats[mtl].Tf;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;
				}

				//bool reverse_normal = false;
				float cos_incident = new_ray.d.dot(rec.n);

				//float cos_a = new_ray.d.dot(rec.n);

				//vec3 nl;
				//if (new_ray.d.dot(rec.n) >= 0)

				vec3 original_normal = rec.n;

				//vec3 original_n = rec.n;

				if (cos_incident > 0)
				{
					//reverse_normal = true;
					rec.n = -rec.n;
				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].isLight)
				{
					if (specular)
					{
						//Le 111
						//L += T * mats[mtl].Ke;

						L += T * Le;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						//Good
						//L += T * mats[mtl].Ke * mis_weight;

						L += T * Le * mis_weight;

						//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
						//return clamp_L;

						return L;

					}
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				if (!bsdfs[mtl]->isSpecular())
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;

					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					
					float cos_mtl = light_direction.dot(rec.n);

					/*if (mats[mtl].type == RoughPlastic_type)
					{
					vec3 h = (-new_ray.d + light_direction).norm();
					cos_mtl = h.dot(light_direction);
					}*/

					if (cos_mtl > 0.0f  && !hit_anything_range(n, light_ray))
					{

						//float ilength = 1.0f / length;
						//light_direction *= ilength;

						//cos_light *= ilength;
						//cos_mtl *= ilength;

						float cos_light = abs(-light_direction.dot(Ls.normal[li]));

						//better than use cos_mtl
						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);

						//float bsdf_pdf;
						//vec3 bsdf_eval = bsdfs[mtl]->eval_and_pdf(new_ray.d, light_direction, rec.n, color, bsdf_pdf);

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						//L += T * mats[mtl].Ke * bsdf_eval * mis / pdf_light;

						L += T * Le * bsdf_eval * mis / pdf_light;

					}


					Sampling_Struct sampling;
					bool isReflect;

					//vec3 sample_normal = rec.n;

					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	sample_normal = -rec.n;

					/*vec3 original_normal = rec.n;

					if (bsdfs[mtl]->have_refraction() && reverse_normal)
					original_normal = -rec.n;
					*/

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					if (!sampling.isSpecular)
					{
						//onb local(rec.n);

						//vec3 direction = (local.u * sampling.direction.x + local.v * sampling.direction.y + local.w * sampling.direction.z).norm();
						//new_ray = Ray(new_hit_point, direction);

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = false;
					}
					else
					{

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = true;
					}
					T *= sampling.weight;
					//T *= sampling.weight;
				}
				else
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;

					Sampling_Struct sampling;
					bool isReflect;



					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	original_normal = -rec.n;


					//bool b;
					//if (!bsdfs[mtl]->isSampleVolume())
					//b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					//Quan Trọng
					//Importance
					//sampling.weight và sampling.eval là hòan tòan khácc nhau

					//vú dị: Diffuse

					//Sampling.weight = color
					//Sampling.eval = color * cos_o * ipi

					//sampling weight tính với tia sample cho bước tiếp theo
					//sampling.eval tính với tia sample light
					//nếu dùng công thức sampling.eval cho sampling.weight thì tia sẽ yếu
					//ảnh sẽ tối và russian roulette sẽ kết thúc sớm
					//dẫn đến ảnh tối hơn và kết thúc rất sớm

					T *= sampling.weight;
				}
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				L += T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

					   //vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
					   //return clamp_L;

				return L;

			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_artificial_light_normal_Ke(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;


		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_alpha_pointer(n, new_ray, rec, &color, alpha))
			float alpha = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, alpha))
			{
				if (randf() < alpha)
				{
					int mtl = trs[rec.ti].mtl;
					T *= mats[mtl].Tf;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;
				}

				//bool reverse_normal = false;
				float cos_incident = new_ray.d.dot(rec.n);

				//float cos_a = new_ray.d.dot(rec.n);

				//vec3 nl;
				//if (new_ray.d.dot(rec.n) >= 0)

				vec3 original_normal = rec.n;

				//vec3 original_n = rec.n;

				if (cos_incident > 0)
				{
					//reverse_normal = true;
					rec.n = -rec.n;
				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{
					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (specular)
					{
						//Le 111
						//L += T * mats[mtl].Ke;

						L += T * Ke;

						return L;
					}
					else
					{
						//original
						/*
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						//Good
						//L += T * mats[mtl].Ke * mis_weight;

						L += T * Ke * mis_weight;

						//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
						//return clamp_L;

						return L;
						*/
						if (mats[mtl].isLight)
						{
							int li = light_map[rec.ti];

							vec3 light_direction(new_ray.d * rec.t);

							float length2 = light_direction.length2();

							light_direction.normalize();


							float cos_light = -light_direction.dot(rec.n);

							float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

							float mis_weight = mis2(prev_pdf, pdf_light);

							//Good
							//L += T * mats[mtl].Ke * mis_weight;

							L += T * Ke * mis_weight;

							//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
							//return clamp_L;

							return L;
						}
						else
						{
							L += T * Ke;
							return L;

							//use this and no light appear
							//float pdf_light = 1.0f;
							//float mis_weight = mis2(prev_pdf, pdf_light);
							//L += T * Ke * mis_weight;

						}
					}
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				if (!bsdfs[mtl]->isSpecular())
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);

					/*if (mats[mtl].type == RoughPlastic_type)
					{
					vec3 h = (-new_ray.d + light_direction).norm();
					cos_mtl = h.dot(light_direction);
					}*/

					if (cos_mtl > 0.0f && !hit_anything_range(n, light_ray))
					{

						//float ilength = 1.0f / length;
						//light_direction *= ilength;

						//cos_light *= ilength;

						//cos_mtl *= ilength;

						float cos_light = -light_direction.dot(Ls.normal[li]);
						cos_light = abs(cos_light);

						//better than use cos_mtl
						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);

						//float bsdf_pdf;
						//vec3 bsdf_eval = bsdfs[mtl]->eval_and_pdf(new_ray.d, light_direction, rec.n, color, bsdf_pdf);

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * Ls.Ke[li] * bsdf_eval * mis / pdf_light;

					}


					Sampling_Struct sampling;
					bool isReflect;

					//vec3 sample_normal = rec.n;

					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	sample_normal = -rec.n;

					/*vec3 original_normal = rec.n;

					if (bsdfs[mtl]->have_refraction() && reverse_normal)
					original_normal = -rec.n;
					*/

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					if (!sampling.isSpecular)
					{
						//onb local(rec.n);

						//vec3 direction = (local.u * sampling.direction.x + local.v * sampling.direction.y + local.w * sampling.direction.z).norm();
						//new_ray = Ray(new_hit_point, direction);

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = false;
					}
					else
					{

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = true;
					}
					T *= sampling.weight;
					//T *= sampling.weight;
				}
				else
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;

					Sampling_Struct sampling;
					bool isReflect;



					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	original_normal = -rec.n;


					//bool b;
					//if (!bsdfs[mtl]->isSampleVolume())
					//b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					//Quan Trọng
					//Importance
					//sampling.weight và sampling.eval là hòan tòan khácc nhau

					//vú dị: Diffuse

					//Sampling.weight = color
					//Sampling.eval = color * cos_o * ipi

					//sampling weight tính với tia sample cho bước tiếp theo
					//sampling.eval tính với tia sample light
					//nếu dùng công thức sampling.eval cho sampling.weight thì tia sẽ yếu
					//ảnh sẽ tối và russian roulette sẽ kết thúc sớm
					//dẫn đến ảnh tối hơn và kết thúc rất sớm

					T *= sampling.weight;
				}
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//L += T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

				//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
				//return clamp_L;

				//return L;

				return L;
			}
		}
		//return L;
		return L;
	}

	vec3 __fastcall sibenik_tracing_artificial_light_normal_Ke_alpha(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;


		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_alpha_pointer(n, new_ray, rec, &color, alpha))
			float alpha = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, alpha))
			{
				
				if (randf() < alpha)
				{
					int mtl = trs[rec.ti].mtl;
					T *= mats[mtl].Tf;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;
				}

				//bool reverse_normal = false;
				float cos_incident = new_ray.d.dot(rec.n);

				//float cos_a = new_ray.d.dot(rec.n);

				//vec3 nl;
				//if (new_ray.d.dot(rec.n) >= 0)

				vec3 original_normal = rec.n;

				//vec3 original_n = rec.n;

				/*if (cos_incident > 0)
				{
					//reverse_normal = true;
					rec.n = -rec.n;
				}*/

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].isLight)// || mats[mtl].type == Light_Diffuse_type)
				{
					//cout << "Light here!\n";
					if (specular)
					{
						//Le 111
						L += T * mats[mtl].Ke;


						return L;
					}
					else
					{
						//if (mats[mtl].isLight)
						//{
						/*if (light_map.find(rec.ti) == light_map.end())
						{
						cout << "Light not map\n";
						cout << rec.ti << "\n";
						cout << mats[mtl].Matname << "\n";
						}*/
						int li = light_map[rec.ti];

						//cout << rec.ti << " " << li << "\n";
						//if (li != 0)
						//cout << rec.ti << " " << li << "\n";
						//	cout << li << "\n";
						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						//Good
						L += T * mats[mtl].Ke * mis_weight;

						//L += T * Ls.Ke[li] * mis_weight;


						//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
						//return clamp_L;

						return L;

						//return vec3(1.0f, 0, 0);
						/*}
						else
						{
						L += T * mats[mtl].Ke;
						return L;
						}*/
					}
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				if (!bsdfs[mtl]->isSpecular())
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					//float cos_light = -light_direction.dot(Ls.normal[li]);
					float cos_mtl = light_direction.dot(rec.n);

					/*if (mats[mtl].type == RoughPlastic_type)
					{
					vec3 h = (-new_ray.d + light_direction).norm();
					cos_mtl = h.dot(light_direction);
					}*/

					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{

						//float ilength = 1.0f / length;
						//light_direction *= ilength;

						//cos_light *= ilength;

						//cos_mtl *= ilength;

						//cos_light = abs(cos_light)

						float cos_light = abs(-light_direction.dot(Ls.normal[li]));

						//better than use cos_mtl
						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);

						//float bsdf_pdf;
						//vec3 bsdf_eval = bsdfs[mtl]->eval_and_pdf(new_ray.d, light_direction, rec.n, color, bsdf_pdf);

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						//failed at this line for dabrovic
						//work for the rest
						//if (li != 0 && li != 1)
						//	cout << li << "\n";

						//vi sao T * mats[mtl].ke co ket qua
						//con T * Ls.ke[li] va T * 4 thi sai

						//cout << Ls.Ke[li].maxc() << "\n";

						//cout << mats[mtl].Ke.maxc() << "\n";

						//L += T * mats[mtl].Ke * bsdf_eval * mis / pdf_light;
						//cout << Ls.Ke.size() << "\n";
						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// *4.0f;
																		 //L = L * 4.0f;
					}


					Sampling_Struct sampling;
					bool isReflect;

					//vec3 sample_normal = rec.n;

					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	sample_normal = -rec.n;

					/*vec3 original_normal = rec.n;

					if (bsdfs[mtl]->have_refraction() && reverse_normal)
					original_normal = -rec.n;
					*/

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					if (!sampling.isSpecular)
					{
						//onb local(rec.n);

						//vec3 direction = (local.u * sampling.direction.x + local.v * sampling.direction.y + local.w * sampling.direction.z).norm();
						//new_ray = Ray(new_hit_point, direction);

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = false;
					}
					else
					{

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = true;
					}
					T *= sampling.weight;
					//T *= sampling.weight;
				}
				else
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;

					Sampling_Struct sampling;
					bool isReflect;



					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	original_normal = -rec.n;


					//bool b;
					//if (!bsdfs[mtl]->isSampleVolume())
					//b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					if (mats[mtl].type == Transparent_type)
					{
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, sampling.direction);
						continue;
					}

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					//Quan Trọng
					//Importance
					//sampling.weight và sampling.eval là hòan tòan khácc nhau

					//vú dị: Diffuse

					//Sampling.weight = color
					//Sampling.eval = color * cos_o * ipi

					//sampling weight tính với tia sample cho bước tiếp theo
					//sampling.eval tính với tia sample light
					//nếu dùng công thức sampling.eval cho sampling.weight thì tia sẽ yếu
					//ảnh sẽ tối và russian roulette sẽ kết thúc sớm
					//dẫn đến ảnh tối hơn và kết thúc rất sớm

					T *= sampling.weight;
				}
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//L += T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

				//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
				//return clamp_L;

				//return L;

				return L;
			}
		}
		//return L;
		return L;
	}

	vec3 __fastcall sibenik_tracing_artificial_light_normal_Ke_alpha_Transparent(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;


		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_alpha_pointer(n, new_ray, rec, &color, alpha))
			float alpha = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, alpha))
			{
				//return vec3(abs(rec.n.x), abs(rec.n.y), abs(rec.n.z));
				//return rec.n;

				float cos_incident = new_ray.d.dot(rec.n);

				//float cos_a = new_ray.d.dot(rec.n);

				//vec3 nl;
				//if (new_ray.d.dot(rec.n) >= 0)

				vec3 original_normal = rec.n;

				//vec3 original_n = rec.n;

				if (cos_incident >= 0)
				{
					//reverse_normal = true;
					rec.n = -rec.n;
				}

				if (randf() < alpha)//Diffuse
				{
					/*int mtl = trs[rec.ti].mtl;
					T *= mats[mtl].Tf;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;*/
					int mtl = trs[rec.ti].mtl;
					T *= color;
				 	vec3 direction_out = cosine_hemisphere_Sampling_with_n(-rec.n);

					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, direction_out);

					prev_pdf = direction_out.dot(-rec.n) * ipi;

					continue;
				}
				/*else//Transparent
				{
					int mtl = trs[rec.ti].mtl;
					T *= mats[mtl].Tf;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;
				}*/

				//bool reverse_normal = false;
				

				//return rec.n;
				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{
					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;
					//cout << "Light here!\n";
					if (specular)
					{
						//Le 111
						L += T * Ke;//mats[mtl].Ke;


						return L;
					}
					else
					{
						//if (mats[mtl].isLight)
						//{
						/*if (light_map.find(rec.ti) == light_map.end())
						{
						cout << "Light not map\n";
						cout << rec.ti << "\n";
						cout << mats[mtl].Matname << "\n";
						}*/
						int li = light_map[rec.ti];

						//cout << rec.ti << " " << li << "\n";
						//if (li != 0)
						//cout << rec.ti << " " << li << "\n";
						//	cout << li << "\n";
						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						//Good
						L += T * Ke * mis_weight;
						//L += T * mats[mtl].Ke * mis_weight;

						//L += T * Ls.Ke[li] * mis_weight;


						//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
						//return clamp_L;

						return L;

						//return vec3(1.0f, 0, 0);
						/*}
						else
						{
						L += T * mats[mtl].Ke;
						return L;
						}*/
					}
				}

				
				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				if (!bsdfs[mtl]->isSpecular())
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					//float cos_light = -light_direction.dot(Ls.normal[li]);
					float cos_mtl = light_direction.dot(rec.n);

					/*if (mats[mtl].type == RoughPlastic_type)
					{
					vec3 h = (-new_ray.d + light_direction).norm();
					cos_mtl = h.dot(light_direction);
					}*/

					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{

						//float ilength = 1.0f / length;
						//light_direction *= ilength;

						//cos_light *= ilength;

						//cos_mtl *= ilength;

						//cos_light = abs(cos_light)

						float cos_light = abs(-light_direction.dot(Ls.normal[li]));

						//better than use cos_mtl
						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);

						//float bsdf_pdf;
						//vec3 bsdf_eval = bsdfs[mtl]->eval_and_pdf(new_ray.d, light_direction, rec.n, color, bsdf_pdf);

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						//failed at this line for dabrovic
						//work for the rest
						//if (li != 0 && li != 1)
						//	cout << li << "\n";

						//vi sao T * mats[mtl].ke co ket qua
						//con T * Ls.ke[li] va T * 4 thi sai

						//cout << Ls.Ke[li].maxc() << "\n";

						//cout << mats[mtl].Ke.maxc() << "\n";

						//L += T * mats[mtl].Ke * bsdf_eval * mis / pdf_light;
						//cout << Ls.Ke.size() << "\n";
						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// *4.0f;
																						//L = L * 4.0f;
					}


					Sampling_Struct sampling;
					bool isReflect;

					//vec3 sample_normal = rec.n;

					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	sample_normal = -rec.n;

					/*vec3 original_normal = rec.n;

					if (bsdfs[mtl]->have_refraction() && reverse_normal)
					original_normal = -rec.n;
					*/

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					if (!sampling.isSpecular)
					{
						//onb local(rec.n);

						//vec3 direction = (local.u * sampling.direction.x + local.v * sampling.direction.y + local.w * sampling.direction.z).norm();
						//new_ray = Ray(new_hit_point, direction);

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = false;
					}
					else
					{

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = true;
					}
					T *= sampling.weight;
					//T *= sampling.weight;
				}
				else
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;

					Sampling_Struct sampling;
					bool isReflect;



					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	original_normal = -rec.n;


					//bool b;
					//if (!bsdfs[mtl]->isSampleVolume())
					//b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					/*if (mats[mtl].type == Transparent_type)
					{
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, sampling.direction);
						continue;
					}*/

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					//Quan Trọng
					//Importance
					//sampling.weight và sampling.eval là hòan tòan khácc nhau

					//vú dị: Diffuse

					//Sampling.weight = color
					//Sampling.eval = color * cos_o * ipi

					//sampling weight tính với tia sample cho bước tiếp theo
					//sampling.eval tính với tia sample light
					//nếu dùng công thức sampling.eval cho sampling.weight thì tia sẽ yếu
					//ảnh sẽ tối và russian roulette sẽ kết thúc sớm
					//dẫn đến ảnh tối hơn và kết thúc rất sớm

					T *= sampling.weight;
				}
				if (i > 3)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				L += T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

				//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
				//return clamp_L;

				//return L;

				return L;
			}
		}
		//return L;
		return L;
	}

	vec3 __fastcall sibenik_tracing_artificial_light_normal_Ke_alpha_Transparent_to_BackGround(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_alpha_pointer(n, new_ray, rec, &color, alpha))
			float alpha = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, alpha))
			{
				//return vec3(abs(rec.n.x), abs(rec.n.y), abs(rec.n.z));
				//return rec.n;

				float cos_incident = new_ray.d.dot(rec.n);

				
				vec3 original_normal = rec.n;

				//vec3 original_n = rec.n;

				if (cos_incident >= 0)				
					rec.n = -rec.n;
			
				if (randf() < alpha)//Diffuse
				{
					
					int mtl = trs[rec.ti].mtl;
					//T *= color;
					vec3 direction_out = cosine_hemisphere_Sampling_with_n(-rec.n);

					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, direction_out);

					prev_pdf *= alpha;

					T *= alpha;
					
					continue;
				}
				T *= (1.0f - alpha);
				//else
				//	return T;
				/*else
				{
					vec3 direction_out = cosine_hemisphere_Sampling_with_n(-rec.n);

					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, direction_out);

					prev_pdf = direction_out.dot(rec.n) * ipi;

					continue;
				}*/
				/*else//Transparent
				{
				int mtl = trs[rec.ti].mtl;
				T *= mats[mtl].Tf;
				new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
				continue;
				}*/

				//bool reverse_normal = false;


				//return rec.n;
				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].type == Light_Transparent_type)
				{
					//cout << 1;
					if (randf() < 0.9f)//mats[mtl].Tr)
					{
						//T *= color;// *0.9f;
						
						//cout << "See Through\n";
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.02f, new_ray.d);
						prev_pdf = 1;
						continue;
					}
					else
					{					
						vec3 Ke = mats[mtl].Ke;
					
						if (specular)
						{							
							L += T * Ke;
							return L;
						}
						else
						{
							
							int li = light_map[rec.ti];
				
							vec3 light_direction(new_ray.d * rec.t);

							float length2 = light_direction.length2();

							light_direction.normalize();

							float cos_light = -light_direction.dot(rec.n);

							float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

							float mis_weight = mis2(prev_pdf, pdf_light);
						
							L += T * Ke * mis_weight;
						
							return L;

						}
					}
				}

				if (mats[mtl].type == Back_Ground_type)
				{
					cout << "BackGround\n";
					vec3 Ke = color;

					L += T * Ke;
					return L;

					/*if (specular)
					{
						L += T * Ke;
						return L;
					}
					else
					{

						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();

						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;

						return L;

					}*/
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				if (!bsdfs[mtl]->isSpecular())
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;
					/*int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					//float cos_light = -light_direction.dot(Ls.normal[li]);
					float cos_mtl = light_direction.dot(rec.n);

					

					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{

						//float ilength = 1.0f / length;
						//light_direction *= ilength;

						//cos_light *= ilength;

						//cos_mtl *= ilength;

						//cos_light = abs(cos_light)

						float cos_light = abs(-light_direction.dot(Ls.normal[li]));

						//better than use cos_mtl
						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);

						//float bsdf_pdf;
						//vec3 bsdf_eval = bsdfs[mtl]->eval_and_pdf(new_ray.d, light_direction, rec.n, color, bsdf_pdf);

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						//failed at this line for dabrovic
						//work for the rest
						//if (li != 0 && li != 1)
						//	cout << li << "\n";

						//vi sao T * mats[mtl].ke co ket qua
						//con T * Ls.ke[li] va T * 4 thi sai

						//cout << Ls.Ke[li].maxc() << "\n";

						//cout << mats[mtl].Ke.maxc() << "\n";

						//L += T * mats[mtl].Ke * bsdf_eval * mis / pdf_light;
						//cout << Ls.Ke.size() << "\n";
						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// *4.0f;
																						//L = L * 4.0f;
					}*/


					Sampling_Struct sampling;
					bool isReflect;

					//vec3 sample_normal = rec.n;

					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	sample_normal = -rec.n;

					/*vec3 original_normal = rec.n;

					if (bsdfs[mtl]->have_refraction() && reverse_normal)
					original_normal = -rec.n;
					*/

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					if (!sampling.isSpecular)
					{
						//onb local(rec.n);

						//vec3 direction = (local.u * sampling.direction.x + local.v * sampling.direction.y + local.w * sampling.direction.z).norm();
						//new_ray = Ray(new_hit_point, direction);

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = false;
					}
					else
					{

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = true;
					}
					T *= sampling.weight;
					//T *= sampling.weight;
				}
				else
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;

					Sampling_Struct sampling;
					bool isReflect;



					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	original_normal = -rec.n;


					//bool b;
					//if (!bsdfs[mtl]->isSampleVolume())
					//b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					if (mats[mtl].type == Transparent_type)
					{
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, sampling.direction);
						continue;
					}

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					//Quan Trọng
					//Importance
					//sampling.weight và sampling.eval là hòan tòan khácc nhau

					//vú dị: Diffuse

					//Sampling.weight = color
					//Sampling.eval = color * cos_o * ipi

					//sampling weight tính với tia sample cho bước tiếp theo
					//sampling.eval tính với tia sample light
					//nếu dùng công thức sampling.eval cho sampling.weight thì tia sẽ yếu
					//ảnh sẽ tối và russian roulette sẽ kết thúc sớm
					//dẫn đến ảnh tối hơn và kết thúc rất sớm

					T *= sampling.weight;
				}
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				L += 9.5 * T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

				//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
				//return clamp_L;

				//return L;

				return L;
			}
		}
		//return L;
		return L;
	}

	vec3 __fastcall sibenik_tracing_artificial_light_normal_Ke_alpha_dabrovic(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;


		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_alpha_pointer(n, new_ray, rec, &color, alpha))
			float alpha = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, alpha))
			{
				if (randf() < alpha)
				{
					int mtl = trs[rec.ti].mtl;
					T *= mats[mtl].Tf;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;
				}

				//bool reverse_normal = false;
				float cos_incident = new_ray.d.dot(rec.n);

				//float cos_a = new_ray.d.dot(rec.n);

				//vec3 nl;
				//if (new_ray.d.dot(rec.n) >= 0)

				vec3 original_normal = rec.n;

				//vec3 original_n = rec.n;

				if (cos_incident > 0)
				{
					//reverse_normal = true;
					rec.n = -rec.n;
				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{
					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (specular)
					{
						//Le 111
						//L += T * mats[mtl].Ke;

						L += T * Ke;

						return L;
					}
					else
					{
						//original
						/*
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						//Good
						//L += T * mats[mtl].Ke * mis_weight;

						L += T * Ke * mis_weight;

						//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
						//return clamp_L;

						return L;
						*/
						if (mats[mtl].isLight)
						{
							int li = light_map[rec.ti];

							vec3 light_direction(new_ray.d * rec.t);

							float length2 = light_direction.length2();

							light_direction.normalize();


							float cos_light = -light_direction.dot(rec.n);

							float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

							float mis_weight = mis2(prev_pdf, pdf_light);

							//Good
							//L += T * mats[mtl].Ke * mis_weight;

							L += T * Ke * mis_weight;

							//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
							//return clamp_L;

							return L;
						}
						else
						{
							L += T * Ke;
							return L;

							//use this and no light appear
							//float pdf_light = 1.0f;
							//float mis_weight = mis2(prev_pdf, pdf_light);
							//L += T * Ke * mis_weight;

						}
					}
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				if (!bsdfs[mtl]->isSpecular())
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);

					/*if (mats[mtl].type == RoughPlastic_type)
					{
					vec3 h = (-new_ray.d + light_direction).norm();
					cos_mtl = h.dot(light_direction);
					}*/

					//hit anything alpha = failed
					if (cos_mtl > 0.0f && !hit_anything_range(n, light_ray))
					{

						//float ilength = 1.0f / length;
						//light_direction *= ilength;

						//cos_light *= ilength;

						//cos_mtl *= ilength;

						float cos_light = abs(-light_direction.dot(Ls.normal[li]));
						//cos_light = abs(cos_light);

						//better than use cos_mtl
						float pdf_light = (length * length) / (Ls.area[li] * abs(cos_light)) * Ls.pdf[li];//cos_mtl

																										  //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);

						//float bsdf_pdf;
						//vec3 bsdf_eval = bsdfs[mtl]->eval_and_pdf(new_ray.d, light_direction, rec.n, color, bsdf_pdf);

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;

					}


					Sampling_Struct sampling;
					bool isReflect;

					//vec3 sample_normal = rec.n;

					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	sample_normal = -rec.n;

					/*vec3 original_normal = rec.n;

					if (bsdfs[mtl]->have_refraction() && reverse_normal)
					original_normal = -rec.n;
					*/

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					if (!sampling.isSpecular)
					{
						//onb local(rec.n);

						//vec3 direction = (local.u * sampling.direction.x + local.v * sampling.direction.y + local.w * sampling.direction.z).norm();
						//new_ray = Ray(new_hit_point, direction);

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = false;
					}
					else
					{

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = true;
					}
					T *= sampling.weight;
					//T *= sampling.weight;
				}
				else
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;

					Sampling_Struct sampling;
					bool isReflect;



					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	original_normal = -rec.n;


					//bool b;
					//if (!bsdfs[mtl]->isSampleVolume())
					//b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					//Quan Trọng
					//Importance
					//sampling.weight và sampling.eval là hòan tòan khácc nhau

					//vú dị: Diffuse

					//Sampling.weight = color
					//Sampling.eval = color * cos_o * ipi

					//sampling weight tính với tia sample cho bước tiếp theo
					//sampling.eval tính với tia sample light
					//nếu dùng công thức sampling.eval cho sampling.weight thì tia sẽ yếu
					//ảnh sẽ tối và russian roulette sẽ kết thúc sớm
					//dẫn đến ảnh tối hơn và kết thúc rất sớm

					T *= sampling.weight;
				}
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//L += T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

				//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
				//return clamp_L;

				//return L;

				return L;
			}
		}
		//return L;
		return L;
	}

	vec3 __fastcall sibenik_tracing_artificial_light_normal_Ke_alpha_roughness(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;


		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_alpha_pointer(n, new_ray, rec, &color, alpha))
			float alpha = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, alpha))
			{
				if (randf() < alpha)
				{
					int mtl = trs[rec.ti].mtl;
					T *= mats[mtl].Tf;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;
				}

				//bool reverse_normal = false;
				float cos_incident = new_ray.d.dot(rec.n);

				//float cos_a = new_ray.d.dot(rec.n);

				//vec3 nl;
				//if (new_ray.d.dot(rec.n) >= 0)

				vec3 original_normal = rec.n;

				//vec3 original_n = rec.n;

				if (cos_incident > 0)
				{
					//reverse_normal = true;
					rec.n = -rec.n;
				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].isLight)// || mats[mtl].type == Light_Diffuse_type)
				{
					//cout << "Light here!\n";
					if (specular)
					{
						//Le 111
						L += T * mats[mtl].Ke;


						return L;
					}
					else
					{
						//if (mats[mtl].isLight)
						//{
						/*if (light_map.find(rec.ti) == light_map.end())
						{
						cout << "Light not map\n";
						cout << rec.ti << "\n";
						cout << mats[mtl].Matname << "\n";
						}*/
						int li = light_map[rec.ti];

						//cout << rec.ti << " " << li << "\n";
						//if (li != 0)
						//cout << rec.ti << " " << li << "\n";
						//	cout << li << "\n";
						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						//Good
						L += T * mats[mtl].Ke * mis_weight;

						//L += T * Ls.Ke[li] * mis_weight;


						//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
						//return clamp_L;

						return L;

						//return vec3(1.0f, 0, 0);
						/*}
						else
						{
						L += T * mats[mtl].Ke;
						return L;
						}*/
					}
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				if (!bsdfs[mtl]->isSpecular())
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					//float cos_light = -light_direction.dot(Ls.normal[li]);
					float cos_mtl = light_direction.dot(rec.n);

					/*if (mats[mtl].type == RoughPlastic_type)
					{
					vec3 h = (-new_ray.d + light_direction).norm();
					cos_mtl = h.dot(light_direction);
					}*/

					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{

						//float ilength = 1.0f / length;
						//light_direction *= ilength;

						//cos_light *= ilength;

						//cos_mtl *= ilength;

						//cos_light = abs(cos_light)

						float cos_light = abs(-light_direction.dot(Ls.normal[li]));

						//better than use cos_mtl
						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);

						//float bsdf_pdf;
						//vec3 bsdf_eval = bsdfs[mtl]->eval_and_pdf(new_ray.d, light_direction, rec.n, color, bsdf_pdf);

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						//failed at this line for dabrovic
						//work for the rest
						//if (li != 0 && li != 1)
						//	cout << li << "\n";

						//vi sao T * mats[mtl].ke co ket qua
						//con T * Ls.ke[li] va T * 4 thi sai

						//cout << Ls.Ke[li].maxc() << "\n";

						//cout << mats[mtl].Ke.maxc() << "\n";

						//L += T * mats[mtl].Ke * bsdf_eval * mis / pdf_light;
						//cout << Ls.Ke.size() << "\n";
						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// *4.0f;
																		 //L = L * 4.0f;
					}


					Sampling_Struct sampling;
					bool isReflect;

					//vec3 sample_normal = rec.n;

					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	sample_normal = -rec.n;

					/*vec3 original_normal = rec.n;

					if (bsdfs[mtl]->have_refraction() && reverse_normal)
					original_normal = -rec.n;
					*/

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					/*if(!mats[mtl].use_Roughness)
					b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);
					else
					{
					float r = rec.r;
					//cout << r << "\n";
					b = bsdfs[mtl]->sample_roughness(new_ray.d, rec.n, original_normal, sampling, isReflect, color, r);
					}*/

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					if (!sampling.isSpecular)
					{
						//onb local(rec.n);

						//vec3 direction = (local.u * sampling.direction.x + local.v * sampling.direction.y + local.w * sampling.direction.z).norm();
						//new_ray = Ray(new_hit_point, direction);

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = false;
					}
					else
					{

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = true;
					}
					T *= sampling.weight;
					//T *= sampling.weight;
				}
				else
				{
					//if (bsdfs[mtl]->need_original_normal() && reverse_normal)
					//	rec.n = -rec.n;

					Sampling_Struct sampling;
					bool isReflect;



					//if (bsdfs[mtl]->have_refraction() && reverse_normal)
					//	original_normal = -rec.n;


					//bool b;
					//if (!bsdfs[mtl]->isSampleVolume())
					//b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					if (mats[mtl].type == Transparent_type)
					{
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, sampling.direction);
						continue;
					}

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					//Quan Trọng
					//Importance
					//sampling.weight và sampling.eval là hòan tòan khácc nhau

					//vú dị: Diffuse

					//Sampling.weight = color
					//Sampling.eval = color * cos_o * ipi

					//sampling weight tính với tia sample cho bước tiếp theo
					//sampling.eval tính với tia sample light
					//nếu dùng công thức sampling.eval cho sampling.weight thì tia sẽ yếu
					//ảnh sẽ tối và russian roulette sẽ kết thúc sớm
					//dẫn đến ảnh tối hơn và kết thúc rất sớm

					T *= sampling.weight;
				}
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//L += T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

				//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
				//return clamp_L;

				//return L;

				return L;
			}
		}
		//return L;
		return L;
	}


	vec3 __fastcall sibenik_tracing_enviroment(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;


		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_alpha_pointer(n, new_ray, rec, &color, alpha))
			float alpha = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, alpha))
			{
				if (randf() < alpha)
				{
					int mtl = trs[rec.ti].mtl;
					T *= mats[mtl].Tf;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				//float cos_a = new_ray.d.dot(rec.n);

				//vec3 nl;
				//if (new_ray.d.dot(rec.n) >= 0)

				if (cos_incident > 0)
					rec.n = -rec.n;

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
				int mtl = trs[rec.ti].mtl;

				/*if (mats[mtl].isLight)
				{
				float u, v;

				Env.DirectionToUV(new_ray.d, u, v);

				vec3 env_value = Env.UVToValue(u, v);

				if (specular)
				{
				//Le 111
				//L += T * mats[mtl].Ke;
				L += T * env_value;
				return L;
				}
				else
				{
				int li = light_map[rec.ti];

				vec3 light_direction(new_ray.d * rec.t);

				float length2 = light_direction.length2();

				light_direction.normalize();


				float cos_light = -light_direction.dot(rec.n);

				float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

				float mis_weight = mis2(prev_pdf, pdf_light);

				//Good
				L += T * env_value * mis_weight;

				return L;

				}
				}*/

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				if (!bsdfs[mtl]->isSpecular())
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					Ray light_ray(hit_point, light_direction);


					
					float cos_mtl = light_direction.dot(rec.n);

					if (cos_mtl > 0.0f  && !hit_anything_range(n, light_ray))
						//if (cos_mtl > 0.0f && !hit_anything_range(n, light_ray))
					{
						//if (use_eviroment)
						//{

						float cos_light = abs(-light_direction.dot(Ls.normal[li]));

						float ilength =1.0f / length;
						light_direction *= ilength;

						cos_light *= ilength;
						cos_mtl *= ilength;

						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						/*
						float u, v;

						Env.DirectionToUV(new_ray.d, u, v);

						vec3 env_value = Env.UVToValue(u, v);
						float env_pdf = Env.UVToPdf(u, v);

						//too bright
						//float mis = mis2(prev_pdf, env_pdf);
						//L += T * env_value * mis;

						//noisy
						//Enviroment_1
						//float mis = mis2(pdf_light * env_pdf, bsdf_pdf * cont_prob);
						//L += T * env_value * bsdf_eval * mis / (pdf_light * env_pdf);
						*/


						float u, v;

						Env.DirectionToUV(new_ray.d, u, v);

						u += rotation_u;

						vec3 env_value = Env.UVToValue(u, v);

						//good
						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);
						L += T * mats[mtl].Tf * env_value * bsdf_eval * mis / (pdf_light);


						//L += T * Le * bsdf_eval * mis / (pdf_light * env_pdf);
						//}
						/*else
						{
						float ilength = 1.0f / length;
						light_direction *= ilength;

						cos_light *= ilength;
						cos_mtl *= ilength;


						//better than use cos_mtl
						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

						//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * Le * bsdf_eval * mis / pdf_light;
						}*/
					}


					Sampling_Struct sampling;
					bool isReflect;

					bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					if (!sampling.isSpecular)
					{
						//onb local(rec.n);

						//vec3 direction = (local.u * sampling.direction.x + local.v * sampling.direction.y + local.w * sampling.direction.z).norm();

						vec3 direction = sampling.direction;

						new_ray = Ray(new_hit_point, direction);

						specular = false;
					}
					else
					{

						new_ray = Ray(new_hit_point, sampling.direction);
						specular = true;
					}
					T *= sampling.weight;
					//T *= sampling.weight;
				}
				else
				{
					Sampling_Struct sampling;
					bool isReflect;

					bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					//Quan Trọng
					//Importance
					//sampling.weight và sampling.eval là hòan tòan khácc nhau

					//vú dị: Diffuse

					//Sampling.weight = color
					//Sampling.eval = color * cos_o * ipi

					//sampling weight tính với tia sample cho bước tiếp theo
					//sampling.eval tính với tia sample light
					//nếu dùng công thức sampling.eval cho sampling.weight thì tia sẽ yếu
					//ảnh sẽ tối và russian roulette sẽ kết thúc sớm
					//dẫn đến ảnh tối hơn và kết thúc rất sớm

					T *= sampling.weight;
				}
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{
				//also correct
				if (use_eviroment)
				{

					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					//if (specular)
					//{
					//	L += T * env_value;
					//	return L;
					//}
					//else
					//{
					float env_pdf = Env.UVToPdf(u, v);

					float mis = mis2(prev_pdf, env_pdf);

					L += T * env_value * mis;

					return L;
					//}
				}
				else
				{
					L += T;
					return L;
				}

				//correct
				/*
				float u, v;

				Env.DirectionToUV(new_ray.d, u, v);

				u += rotation_u;

				vec3 env_value = Env.UVToValue(u, v);

				if (specular)
				{
				L += T * env_value;
				return L;
				}
				else
				{
				float env_pdf = Env.UVToPdf(u, v);

				float mis = mis2(prev_pdf, env_pdf);

				L += T * env_value * mis;

				return L;
				}
				*/
			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		//vec3 Ke = 

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				//return color;
				//int mtl = trs[rec.ti].mtl;

				//if (mats[mtl].type == Light_Diffuse_type)
					//return rec.n;
				//else
				//	return vec3(0.0f);

				//if (Tr < 1.0f)
				//	Tr = 0.0f;
				//if (Tr > 0.0f && Tr != 1.0f)
				//	cout << Tr << "\n";
				//if (Tr != 0.0f && Tr != 1.0f)
				//	cout << Tr << "\n";
				if (randf() <= Tr)
				{
					//int mtl = trs[rec.ti].mtl;
					//T *= mats[mtl].Tf;
					//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);

					//* Tr lam mat cai rem
					//T *= Tr;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					continue;
				}
				/*if (Tr == 1.0f)
				{
					//cout << "a\n";
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					continue;
				}*/



				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;


				if (cos_incident > 0)
					rec.n = -rec.n;


				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].type == Back_Ground_type)
				{
					L += T * 4 * color;
					return L;
				}

				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{
				
					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;
					//if (mats[mtl].isLight)
					//{
						//mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{
						//float scale_dist = rec.t * randf() * 0.01;

						L += T * Ke * 6;//(1.0f + expf(scale_dist));

						return L;
					}

						if (specular || mats[mtl].type == Light_Diffuse_type)
						{

							L += T * Ke;// *0.1f;

							return L;
						}
						else
						{
							int li = light_map[rec.ti];

							vec3 light_direction(new_ray.d * rec.t);

							float length2 = light_direction.length2();

							light_direction.normalize();


							float cos_light = -light_direction.dot(rec.n);

							float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

							float mis_weight = mis2(prev_pdf, pdf_light);

							//Good
							//L += T * mats[mtl].Ke * mis_weight;

							L += T * Ke * mis_weight;// *0.1f;


							return L;

						}
					//}
					/*else
					{
						//return vec3(0.0f);
						L += T * mats[mtl].Kd;// *0.1f;// *0.0000001f;
						return L;
					}*/
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				//sample direct lighting				

				//if (!bsdfs[mtl]->isSpecular())
				//{
					if (use_directional_light)
					{
						vec3 sun_direction = inifinite_light_direction;
						if (blur_directional_light)
						{
							//from hit point move 1 unit in light direction
							//create a sphere of radius equal to blur radius
							//sample random point on sphere, minus these point with hitpoint
							//this will be the new sun_direction

							vec3 blur_center = hit_point + sun_direction;

							float u = randf();
							float v = randf();

							float theta = tau * u;
							float phi = acosf(2.0f * v - 1.0f);

							float sin_phi, cos_phi;

							FTA::sincos(phi, &sin_phi, &cos_phi);

							float sin_theta, cos_theta;

							FTA::sincos(theta, &sin_theta, &cos_theta);

							float x = directional_blur_radius * sin_phi * cos_theta;
							float y = directional_blur_radius * sin_phi * sin_theta;
							float z = directional_blur_radius * cos_phi;

							//original

							vec3 sample_point = blur_center + vec3(x, y, z);
							sun_direction = (sample_point - hit_point).norm();

							//wrong
							//sun_direction = vec3(x, y, z);
							//sun_direction.normalize();

						}
						//Ray directional_lighting(hit_point, inifinite_light_direction);
						Ray directional_lighting(hit_point, sun_direction);

						//if (!hit_anything_range(n, directional_lighting))
						if (!hit_anything_alpha(n, directional_lighting))
						{
							vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

							//MIS
							//float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
							//float sun_pdf = 1.0f;//there is only one way to sample this direction

							//float mis = mis2(sun_pdf, directional_pdf * cont_prob);


							//L += T * sun_color * directional_bsdf * mis;

							//No MIS
							//L += T * mats[mtl].Tf * sun_color * directional_bsdf;

						}
					}
					if (use_area_light)
					{
						int li = Ls.sample();
						vec3 light_point(Ls.sample_light(li));

						vec3 light_direction(light_point - hit_point);

						float length = light_direction.length();

						float ilength = 1.0f / length;
						light_direction *= ilength;

						Ray light_ray(hit_point, light_direction);


						float cos_mtl = light_direction.dot(rec.n);


						if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
						{
							float cos_light = abs(-light_direction.dot(Ls.normal[li]));


							float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																										 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

							vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
							float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


							float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

							//cout << li << "\n";

							//cout << Ls.Ke.size() << "\n";
							L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

						}
					}

					Sampling_Struct sampling;
					bool isReflect;


					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					new_ray = Ray(new_hit_point, sampling.direction);

					specular = sampling.isSpecular;


					T *= sampling.weight;
				//}
				/*else
				{

					Sampling_Struct sampling;
					bool isReflect;

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					T *= sampling.weight;
				}*/
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//L += T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

				//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
				//return clamp_L;

				//return L;

				//if use enviroment compute this value

				//compute sky color will create better result

				//this is from Ray Tracing In One Weekend


				//if (i == 0)
				//{	
					T *= color;
					if (compute_sky_color)
					{
						/*vec3 sky_direction = new_ray.d;

						float t = 0.5f * (sky_direction.y + 1.0f);

						//vec3 sky_color = (1.0f - t) * vec3(1.0f) + t * vec3(0.5f, 0.7f, 4.0f);

						vec3 sky_color = (1.0f - t) * vec3(1.0f) + t * vec3(0.2f, 0.4f, 6.0f);

						L += T * sky_color;*/

						L += T * 12.0f;

						return L;
					}
					else
					{
						float u, v;

						Env.DirectionToUV(new_ray.d, u, v);

						u += rotation_u;

						vec3 env_value = Env.UVToValue(u, v);

						L += T * env_value;
						return L;

						/*//if (specular)
						//{
						//	L += T * env_value;
						//	return L;
						//}
						//else
						//{
						float env_pdf = Env.UVToPdf(u, v);

						float mis = mis2(prev_pdf, env_pdf);

						L += T * env_value * mis;

						return L;*/
					}
				//}
				/*else
				{
					//return color;
					T *= color;
					L += 2 * T;// *0.1f;// *0.00001f;
					return L;
				}*/
				/*if (specular)
				{

				L += T * sky_color;

				return L;
				}
				else
				{
				int li = light_map[rec.ti];

				vec3 light_direction(new_ray.d * rec.t);

				float length2 = light_direction.length2();

				light_direction.normalize();


				float cos_light = -light_direction.dot(rec.n);

				float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

				float mis_weight = mis2(prev_pdf, pdf_light);

				//Good
				//L += T * mats[mtl].Ke * mis_weight;

				L += T * sky_color * mis_weight;

				return L;

				}*/

			}
		}
		//return L;
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_2(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{	
				int mtl = trs[rec.ti].mtl;

				//if (mtl == 341)
				//{
				//	return L += T * 4;
				//}

				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				
				/*if (Tr > 0)
				{
					if (randf() <= Tr)
					{
						T *= Tr;
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{
						L += 0.2 * T;
						return L;
					}
				}*/

				if (Tr > 0 && randf() <= Tr)
				{
					
					T *= Tr;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					continue;
					
				}
				

				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;
					
					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{
						//float scale_dist = rec.t * randf() * 0.01;

						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}
					
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				
				if (use_directional_light)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{
						
						vec3 blur_center = hit_point + sun_direction;

						float u = randf();
						float v = randf();

						float theta = tau * u;
						float phi = acosf(2.0f * v - 1.0f);

						float sin_phi, cos_phi;

						FTA::sincos(phi, &sin_phi, &cos_phi);

						float sin_theta, cos_theta;

						FTA::sincos(theta, &sin_theta, &cos_theta);

						float blur_radius = directional_blur_radius, scale_sun = 1.0f;

						if (randf() < 0.1f)
						{
							blur_radius = small_blur_radius;
							scale_sun = 6.0f;
						}

						/*float x = directional_blur_radius * sin_phi * cos_theta;
						float y = directional_blur_radius * sin_phi * sin_theta;
						float z = directional_blur_radius * cos_phi;*/

						float x = blur_radius * sin_phi * cos_theta;
						float y = blur_radius * sin_phi * sin_theta;
						float z = blur_radius * cos_phi;
						
						vec3 sample_point = blur_center + vec3(x, y, z);
						sun_direction = (sample_point - hit_point).norm();

					//}
					
					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					//HitRecord new_rec;

					/*if (hit_color_alpha(n, directional_lighting, new_rec, &color, Tr))
					{
						
						if (randf() <= Tr)
						{
							vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

							float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
							float sun_pdf = 1.0f;//there is only one way to sample this direction

							float mis = mis2(sun_pdf, directional_pdf * cont_prob);


							L += T * sun_color * directional_bsdf * mis;

							//L += T * sun_color * directional_bsdf;
							//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
							//continue;
						}
					}*/
					///*else
					//vec3 new_color;
					//int new_mtl = trs[new_rec.ti].mtl;
					//if (!hit_color_alpha(n, directional_lighting, new_rec, &new_color, Tr))
					if(!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

						
						L += T * sun_color * directional_bsdf * scale_sun;

					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				//sum up light

				int ind = randf() * round_light[0].size();

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;
				
				if (i > 40)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{
					

					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}
				

			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_circle(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				int mtl = trs[rec.ti].mtl;

				
				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;
			
				if (Tr > 0 && randf() <= Tr)
				{

					T *= Tr;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					continue;

				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{
					

						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);


				if (use_directional_light)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					/*float x = directional_blur_radius * sin_phi * cos_theta;
					float y = directional_blur_radius * sin_phi * sin_theta;
					float z = directional_blur_radius * cos_phi;*/

					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					//HitRecord new_rec;

					/*if (hit_color_alpha(n, directional_lighting, new_rec, &color, Tr))
					{

					if (randf() <= Tr)
					{
					vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

					float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
					float sun_pdf = 1.0f;//there is only one way to sample this direction

					float mis = mis2(sun_pdf, directional_pdf * cont_prob);


					L += T * sun_color * directional_bsdf * mis;

					//L += T * sun_color * directional_bsdf;
					//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					//continue;
					}
					}*/
					///*else
					//vec3 new_color;
					//int new_mtl = trs[new_rec.ti].mtl;
					//if (!hit_color_alpha(n, directional_lighting, new_rec, &new_color, Tr))
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);


						L += T * sun_color * directional_bsdf * scale_sun;

					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				//sum up light

				//int ind = randf() * round_light[0].size();

			

				//if (hit_point.y <= 332.914)
				//{
					int min_distance = (hit_point - light_center[0]).length2();
					int min_light = 0;

					for (int j = 1; j < 4; ++j)
					{
						float distance = (hit_point - light_center[j]).length2();

						if (distance < min_distance)
						{
							min_distance = distance;
							min_light = j;
						}
					}

					//sample min_light
					//if (min_distance <= 810000)//distance 1000//400
					//{
						//cout << "sample";

						//clock_t a = clock();

						//a = (clock() - a);

						//float seed = ((double)a) / CLOCKS_PER_SEC;

						//srand(seed);

						int ind = randf() * round_light[min_light].size();

						int light_ind = round_light[min_light][ind];
						
						//cout << ind << " " << light_ind << "\n";

						vec3 light_point = sample_triangle(light_ind);

						vec3 light_direction = (light_point - hit_point);

						float length2 = light_direction.length2();

						float length = sqrt14(length2);
					

						float inv_length = 1.0f / length;

		

						light_direction *= inv_length;

						Ray light_ray(hit_point, light_direction);

						//length += 40.0f;
						//if (!hit_anything_alpha(n, light_ray))
						//if(!hit_anything_range_alpha(n, light_ray, length))
						if(!hit_anything_range_2(n, light_ray, length))
						{
							//cout << "Not Hit\n";
							vec3 norm(0, -1, 0);

							float cos_light = -light_direction.dot(norm);

							cos_light = abs(cos_light);

							//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

																										 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

							vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
							//float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


							//float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

							float inv_length2 = 1.0f / length2;
							L += T * round_light_Ke * bsdf_eval * inv_length2;// *mis / pdf_light;// * 0.1f
						}
					//}
				//}


				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_sphere(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				if (Tr > 0 && randf() <= Tr)
				{

					T *= Tr;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					continue;

				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);


				if (use_directional_light)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					/*float x = directional_blur_radius * sin_phi * cos_theta;
					float y = directional_blur_radius * sin_phi * sin_theta;
					float z = directional_blur_radius * cos_phi;*/

					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					//HitRecord new_rec;

					/*if (hit_color_alpha(n, directional_lighting, new_rec, &color, Tr))
					{

					if (randf() <= Tr)
					{
					vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

					float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
					float sun_pdf = 1.0f;//there is only one way to sample this direction

					float mis = mis2(sun_pdf, directional_pdf * cont_prob);


					L += T * sun_color * directional_bsdf * mis;

					//L += T * sun_color * directional_bsdf;
					//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					//continue;
					}
					}*/
					///*else
					//vec3 new_color;
					//int new_mtl = trs[new_rec.ti].mtl;
					//if (!hit_color_alpha(n, directional_lighting, new_rec, &new_color, Tr))
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);


						L += T * sun_color * directional_bsdf * scale_sun;

					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				//sum up light

				//int ind = randf() * round_light[0].size();

				//sample hemisphere
				/*vec3 sphere_norm(0, -1, 0);
				vec3 coord = uniform_hemisphere_Sampling(sphere_norm);

				vec3 sphere_point = local_sphere.u * coord.x + local_sphere.v * coord.y + local_sphere.w * coord.z;

				sphere_point = light_center[0] - vec3(0, 120, 0) + sphere_point * sphere_radius;

				//sphere_point *= sphere_radius;

				vec3 light_direction(sphere_point - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				float inv_length = 1.0f / length;

				light_direction *= inv_length;*/

				vec3 light_position;//ham nay co trong scene ma bo
				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();
				
				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))
				if (!hit_anything_range_2(n, light_ray, length))
				//if(!hit_anything_range_alpha(n, light_ray, length))
				{
					//cout << "Not Hit\n";
					//vec3 norm(0, -1, 0);

					//float cos_light = -light_direction.dot(norm);

					//cos_light = abs(cos_light);

					//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

					//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

					float pdf_light = 1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					//float inv_length2 = 1.0f / length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis / pdf_light;// *mis / pdf_light;// * 0.1f
				}
				//}
				//}


				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_circle_2(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (Tr > 0)
				{
					//T *= color;
					//Good
					/*
					if (randf() < 0.5f)//0.5
					{
						if (randf() <= Tr)
						{
							T *= 0.99f;
							//T *= Tr;
							//T *= expf(-randf() * 0.1f);

							//if (hit_point.y > 2237.5f)
							//	T *= 0.95f;
							//else
							//	T *= 0.8f;
							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
							continue;
						}
						else
						{
							T *= 0.99f;
							vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);//Reflect(new_ray.d, rec.n);


							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, direction_out);
							continue;
						}
					}
					else
					{
						float ior = 1.9f;//1.36

						float F = Dielectric_Reflectance(ior, cos_incident);

						//T *= color;

						if (randf() < F)
						{
							T *= F;
							vec3 direction_out //= Reflect(new_ray.d, rec.n);
											   = cosine_hemisphere_Sampling_with_n(rec.n);

							float cos_o = direction_out.dot(rec.n);

							//T *= cos_o;

							new_ray = Ray(hit_point, direction_out);
							continue;
						}
						else
						{

							T *= (1.0f - F);
							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
							continue;
						}
					}
					*/


					//Good2
					/*
					T *= 2 * expf(-cos_incident);
					if (randf() < 0.5f)//0.5
					{

						if (randf() <= Tr)
						{
							T *= 0.99f;
							//T *= Tr;
							//T *= expf(-randf() * 0.1f);

							//if (hit_point.y > 2237.5f)
							//	T *= 0.95f;
							//else
							//	T *= 0.8f;
							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
							continue;
						}
						else
						{
							T *= 0.99f;
							vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);//Reflect(new_ray.d, rec.n);


							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, direction_out);
							continue;
						}
					}
					else
					{
						float ior = 1.9f;//1.36

						float F = Dielectric_Reflectance(ior, cos_incident);

						//T *= color;

						if (randf() < F)
						{
							T *= F;
							vec3 direction_out //= Reflect(new_ray.d, rec.n);
								= cosine_hemisphere_Sampling_with_n(rec.n);

							float cos_o = direction_out.dot(rec.n);

							//T *= cos_o;

							new_ray = Ray(hit_point, direction_out);
							continue;
						}
						else
						{

							T *= (1.0f - F);
							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
							continue;
						}
					}
					*/

					//Good diffuse
/*if (randf() < 0.8f)
{

	//float blur_radius = 0.02f;
	//vec3 blur_point;

	//random_point_on_sphere(blur_radius, blur_point);

	//vec3 hit_point(new_ray.o + new_ray.d * rec.t * 1.0002f);
	//vec3 direction_point = hit_point + new_ray.d + blur_point;
	//vec3 direction_out = (direction_point - hit_point).norm();


	//optimize

	//vec3 direction_out = (new_ray.d + blur_point).norm();

	new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
	continue;
}
else
{

	vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);

	float cos_o = direction_out.dot(rec.n);


	float sin_o = sqrt14(1.0f - cos_o * cos_o);
	float horizon_scatter = powf(sin_o, 0.5);

	vec3 R(0.5f);

	vec3 weight = color * cos_o + horizon_scatter * R * cos_incident * ipi;

	//T *= color * direction_out.dot(rec.n);

	T *= weight;

	new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
	continue;
}*/

					
					//Good

					if (randf() < 0.5f)
					{
						
						//float blur_radius = 0.02f;
						//vec3 blur_point;

						//random_point_on_sphere(blur_radius, blur_point);

						//vec3 hit_point(new_ray.o + new_ray.d * rec.t * 1.0002f);
						//vec3 direction_point = hit_point + new_ray.d + blur_point;
						//vec3 direction_out = (direction_point - hit_point).norm();


						//optimize

						//vec3 direction_out = (new_ray.d + blur_point).norm();

						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{

						/*vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);

						float cos_o = direction_out.dot(rec.n);
						

						float sin_o = sqrt14(1.0f - cos_o * cos_o);
						float horizon_scatter = powf(sin_o, 0.5);

						vec3 R(0.5f);

						vec3 weight = color * cos_o + horizon_scatter * R * cos_incident * ipi;

						//T *= color * direction_out.dot(rec.n);

						T *= weight;

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;*/

						Sampling_Struct sampling;
						bool isReflect;

						bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

						vec3 direction_out = sampling.direction;
						vec3 weight = sampling.weight;
						float pdf = sampling.pdf;

						T *= weight;
						prev_pdf = pdf;

						//T *= Tr;

						if (isReflect)
						
							new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						
						else
							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, direction_out);
						
						continue;
					}
					

/*float F = Dielectric_Reflectance(1.1, cos_incident);
if (randf() < 1.0f - F)
{
	T *= Tr;
	//float blur_radius = 0.02f;
	//vec3 blur_point;

	//random_point_on_sphere(blur_radius, blur_point);

	//vec3 hit_point(new_ray.o + new_ray.d * rec.t * 1.0002f);
	//vec3 direction_point = hit_point + new_ray.d + blur_point;
	//vec3 direction_out = (direction_point - hit_point).norm();


	//optimize

	//vec3 direction_out = (new_ray.d + blur_point).norm();

	new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
	continue;
}
else
{

	vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);

	//T *= color * direction_out.dot(rec.n);

	new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
	continue;
}*/

				
					//no too good
					/*
					float F = Dielectric_Reflectance(1.4, cos_incident);

					T *= Tr;
					if (randf() < F)
					{
						T *= color;
						vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);//Reflect(new_ray.d, rec.n);

						prev_pdf = direction_out.dot(rec.n);


						new_ray = Ray(hit_point, direction_out);

						continue;
					}
					else
					{
						if (randf() < Tr)
						{
							//T *= Tr;
							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
							continue;
						}
						else
						{
							//T *= expf(-cos_incident);
							vec3 direction_out = cosine_hemisphere_Sampling_with_n(-rec.n);
							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, direction_out);
							continue;
						}
					}
					*/


					
					//failed
					/*T *= 1.5f * expf(-cos_incident);
					if (randf() < 0.5f)//0.5
					{

						if (randf() <= Tr)
						{
							T *= 0.99f;
		
							//vec3 direction_out = cosine_hemisphere_Sampling_with_n(-rec.n);

							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
							continue;
						}
						else
						{
							T *= 0.99f;
							vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);//Reflect(new_ray.d, rec.n);


							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, direction_out);
							continue;
						}
					}
				else
				{
					float ior = 1.1f;//1.36

					float F = Dielectric_Reflectance(ior, cos_incident);

					
					if (randf() < F)
					{
						T *= F;
						vec3 direction_out = Reflect(new_ray.d, rec.n);
						//= cosine_hemisphere_Sampling_with_n(rec.n);

						

						new_ray = Ray(hit_point, direction_out);
						continue;
					}
					else
					{

						T *= (1.0f - F);
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
				}*/
					//failed
					/*T *= 2 * expf(-cos_incident);
					if (randf() <= Tr)
					{
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{
						vec3 direction_out = Reflect(new_ray.d, rec.n);
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, direction_out);
						continue;
					}
					*/


					

					/*if (randf() <= Tr)
					{
						T *= Tr;
						//T *= expf(-randf() * 0.1f);

						//if (hit_point.y > 2237.5f)
						//	T *= 0.95f;
						//else
						//	T *= 0.8f;
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{
						vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);//Reflect(new_ray.d, rec.n);
						
						
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, direction_out);
						continue;
					}
					*/
					/*
					float ior = 1.36f;

					float F = Dielectric_Reflectance(ior, cos_incident);

					T *= color;

					if (randf() < F)
					{
						T *= F;
						vec3 direction_out = Reflect(new_ray.d, rec.n);//cosine_hemisphere_Sampling_with_n(rec.n);
						new_ray = Ray(hit_point, direction_out);
						continue;
					}
					else
					{
						T *= (1.0f - F);
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					*/
					/*float ior = 1.9f;

					float F = Dielectric_Reflectance(ior, cos_incident);

					if (randf() < F)
					{
						T *= F;
						vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);
						new_ray = Ray(hit_point, direction_out);
						continue;
					}
					else
					{
						T *= (1.0f - F);
						//vec3 direction_out = cosine_hemisphere_Sampling_with_n(-rec.n);
						//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, direction_out);

						vec3 blur_center = hit_point - rec.n;

						float u = randf();
						float v = randf();

						float theta = tau * u;
						float phi = acosf(2.0f * v - 1.0f);

						float sin_phi, cos_phi;

						FTA::sincos(phi, &sin_phi, &cos_phi);

						float sin_theta, cos_theta;

						FTA::sincos(theta, &sin_theta, &cos_theta);

						float blur_radius = 0.02f;

					

						
						float x = blur_radius * sin_phi * cos_theta;
						float y = blur_radius * sin_phi * sin_theta;
						float z = blur_radius * cos_phi;

						vec3 sample_point = blur_center + vec3(x, y, z);
						vec3 direction_out = (sample_point - hit_point).norm();

						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, direction_out);

						continue;
					}*/

				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					/*float x = directional_blur_radius * sin_phi * cos_theta;
					float y = directional_blur_radius * sin_phi * sin_theta;
					float z = directional_blur_radius * cos_phi;*/

					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					//HitRecord new_rec;

					/*if (hit_color_alpha(n, directional_lighting, new_rec, &color, Tr))
					{

					if (randf() <= Tr)
					{
					vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

					float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
					float sun_pdf = 1.0f;//there is only one way to sample this direction

					float mis = mis2(sun_pdf, directional_pdf * cont_prob);


					L += T * sun_color * directional_bsdf * mis;

					//L += T * sun_color * directional_bsdf;
					//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					//continue;
					}
					}*/
					///*else
					//vec3 new_color;
					//int new_mtl = trs[new_rec.ti].mtl;
					//if (!hit_color_alpha(n, directional_lighting, new_rec, &new_color, Tr))
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);


						L += T * sun_color * directional_bsdf * scale_sun;

					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				//sum up light

				//int ind = randf() * round_light[0].size();

				//sample hemisphere
				/*vec3 sphere_norm(0, -1, 0);
				vec3 coord = uniform_hemisphere_Sampling(sphere_norm);

				vec3 sphere_point = local_sphere.u * coord.x + local_sphere.v * coord.y + local_sphere.w * coord.z;

				sphere_point = light_center[0] - vec3(0, 120, 0) + sphere_point * sphere_radius;

				//sphere_point *= sphere_radius;

				vec3 light_direction(sphere_point - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				float inv_length = 1.0f / length;

				light_direction *= inv_length;*/

				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))
				
				if (length <= 1000 && hit_point.y < 2798.2f && !hit_anything_range_2(n, light_ray, length))
				
				//if (length < 2000 && !hit_anything_range_2(n, light_ray, length))

					//if(!hit_anything_range_alpha(n, light_ray, length))
				{
					//cout << "Not Hit\n";
					vec3 norm(0, -1, 0);

					float cos_light = -light_direction.dot(norm);

					//cos_light = abs(cos_light);

					//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

					//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}
				

				vec3 random_sphere_point;

				random_point_on_sphere(pic_light_radius, random_sphere_point);

				vec3 picture_light_point = picture_light_position + random_sphere_point;

				vec3 cylinder_light_direction(picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				if (cylinder_length <= 1000 && hit_point.y < 2437.5 && !hit_anything_range_2(n, cylinder_light_ray, cylinder_length2))
				{
					//cout << "Hit Cylinder\n";
					float pdf_light = pdf_light_cylinder;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = cylinder_light_direction.dot(rec.n);

					float inv_length2 = 1.0f / cylinder_length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light;// *100000 * inv_length2;

				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_circle_blur(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				//return color;
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (Tr > 0)
				{					
					//Good diffuse
					if (randf() < 0.8f)
					{

					//float blur_radius = 0.02f;
					//vec3 blur_point;

					//random_point_on_sphere(blur_radius, blur_point);

					//vec3 hit_point(new_ray.o + new_ray.d * rec.t * 1.0002f);
					//vec3 direction_point = hit_point + new_ray.d + blur_point;
					//vec3 direction_out = (direction_point - hit_point).norm();


					//optimize

					//vec3 direction_out = (new_ray.d + blur_point).norm();

						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{

						/*vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);

						float cos_o = direction_out.dot(rec.n);


						float sin_o = sqrt14(1.0f - cos_o * cos_o);
						float horizon_scatter = powf(sin_o, 0.5);

						vec3 R(0.5f);

						vec3 weight = color * cos_o + horizon_scatter * R * cos_incident * ipi;

						//T *= color * direction_out.dot(rec.n);

						T *= weight;

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;*/

						vec3 direction_out = Reflect(new_ray.d, rec.n);

						float blur_radius = 0.4f;
						vec3 blur_point;

						random_point_on_sphere(blur_radius, blur_point);

						direction_out = (direction_out + blur_point).norm();

						float cos_o = direction_out.dot(rec.n);
						T *= powf(abs(cos_o), cos_incident);

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;

					}


					

				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light && mats[mtl].type != Glass_type || mats[mtl].type != ThinDielectric_type)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					/*float x = directional_blur_radius * sin_phi * cos_theta;
					float y = directional_blur_radius * sin_phi * sin_theta;
					float z = directional_blur_radius * cos_phi;*/

					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					//HitRecord new_rec;

					/*if (hit_color_alpha(n, directional_lighting, new_rec, &color, Tr))
					{

					if (randf() <= Tr)
					{
					vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

					float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
					float sun_pdf = 1.0f;//there is only one way to sample this direction

					float mis = mis2(sun_pdf, directional_pdf * cont_prob);


					L += T * sun_color * directional_bsdf * mis;

					//L += T * sun_color * directional_bsdf;
					//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					//continue;
					}
					}*/
					///*else
					//vec3 new_color;
					//int new_mtl = trs[new_rec.ti].mtl;
					//if (!hit_color_alpha(n, directional_lighting, new_rec, &new_color, Tr))
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

						//if(mats[mtl].type != Glass_type || mats[mtl].type != ThinDielectric_type)
							L += T * sun_color * directional_bsdf * scale_sun;
						//else
						//	L += T * vec3(0.4941f, 0.449f, 0.2196f) * directional_bsdf * scale_sun;

						//0.9882, 0.8980, 0.4392


					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				//sum up light

				//int ind = randf() * round_light[0].size();

				//sample hemisphere
				/*vec3 sphere_norm(0, -1, 0);
				vec3 coord = uniform_hemisphere_Sampling(sphere_norm);

				vec3 sphere_point = local_sphere.u * coord.x + local_sphere.v * coord.y + local_sphere.w * coord.z;

				sphere_point = light_center[0] - vec3(0, 120, 0) + sphere_point * sphere_radius;

				//sphere_point *= sphere_radius;

				vec3 light_direction(sphere_point - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				float inv_length = 1.0f / length;

				light_direction *= inv_length;*/

				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))

				if (length <= 1000 && hit_point.y < 2798.2f && !hit_anything_range_2(n, light_ray, length))

					//if (length < 2000 && !hit_anything_range_2(n, light_ray, length))

					//if(!hit_anything_range_alpha(n, light_ray, length))
				{
					//cout << "Not Hit\n";
					vec3 norm(0, -1, 0);

					float cos_light = -light_direction.dot(norm);

					//cos_light = abs(cos_light);

					//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

					//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}


				vec3 random_sphere_point;

				random_point_on_sphere(pic_light_radius, random_sphere_point);

				vec3 picture_light_point = picture_light_position + random_sphere_point;

				vec3 cylinder_light_direction(picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				if (cylinder_length <= 1000 && hit_point.y < 2437.5 && !hit_anything_range_2(n, cylinder_light_ray, cylinder_length2))
				{
					//cout << "Hit Cylinder\n";
					float pdf_light = pdf_light_cylinder;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = cylinder_light_direction.dot(rec.n);

					float inv_length2 = 1.0f / cylinder_length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light;// *100000 * inv_length2;

				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}



	vec3 __fastcall sibenik_tracing_directional_lighting_Rough_Plastic(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (Tr > 0)
				{
					//Good diffuse
					if (randf() < 0.8f)
					{

						//float blur_radius = 0.02f;
						//vec3 blur_point;

						//random_point_on_sphere(blur_radius, blur_point);

						//vec3 hit_point(new_ray.o + new_ray.d * rec.t * 1.0002f);
						//vec3 direction_point = hit_point + new_ray.d + blur_point;
						//vec3 direction_out = (direction_point - hit_point).norm();


						//optimize

						//vec3 direction_out = (new_ray.d + blur_point).norm();

						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{

						

						vec3 direction_out = Reflect(new_ray.d, rec.n);

						float blur_radius = 0.4f;
						vec3 blur_point;

						random_point_on_sphere(blur_radius, blur_point);

						direction_out = (direction_out + blur_point).norm();

						float cos_o = direction_out.dot(rec.n);
						T *= powf(abs(cos_o), cos_incident);

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;

					}




				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					/*float x = directional_blur_radius * sin_phi * cos_theta;
					float y = directional_blur_radius * sin_phi * sin_theta;
					float z = directional_blur_radius * cos_phi;*/

					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					//HitRecord new_rec;

					/*if (hit_color_alpha(n, directional_lighting, new_rec, &color, Tr))
					{

					if (randf() <= Tr)
					{
					vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

					float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
					float sun_pdf = 1.0f;//there is only one way to sample this direction

					float mis = mis2(sun_pdf, directional_pdf * cont_prob);


					L += T * sun_color * directional_bsdf * mis;

					//L += T * sun_color * directional_bsdf;
					//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					//continue;
					}
					}*/
					///*else
					//vec3 new_color;
					//int new_mtl = trs[new_rec.ti].mtl;
					//if (!hit_color_alpha(n, directional_lighting, new_rec, &new_color, Tr))
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);


						L += T * sun_color * directional_bsdf * scale_sun;

					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				//sum up light

				//int ind = randf() * round_light[0].size();

				//sample hemisphere
				/*vec3 sphere_norm(0, -1, 0);
				vec3 coord = uniform_hemisphere_Sampling(sphere_norm);

				vec3 sphere_point = local_sphere.u * coord.x + local_sphere.v * coord.y + local_sphere.w * coord.z;

				sphere_point = light_center[0] - vec3(0, 120, 0) + sphere_point * sphere_radius;

				//sphere_point *= sphere_radius;

				vec3 light_direction(sphere_point - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				float inv_length = 1.0f / length;

				light_direction *= inv_length;*/

				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))

				if (length <= 1000 && hit_point.y < 2798.2f && !hit_anything_range_2(n, light_ray, length))

					//if (length < 2000 && !hit_anything_range_2(n, light_ray, length))

					//if(!hit_anything_range_alpha(n, light_ray, length))
				{
					//cout << "Not Hit\n";
					vec3 norm(0, -1, 0);

					float cos_light = -light_direction.dot(norm);

					//cos_light = abs(cos_light);

					//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

					//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}


				vec3 random_sphere_point;

				random_point_on_sphere(pic_light_radius, random_sphere_point);

				vec3 picture_light_point = picture_light_position + random_sphere_point;

				vec3 cylinder_light_direction(picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				if (cylinder_length <= 1000 && hit_point.y < 2437.5 && !hit_anything_range_2(n, cylinder_light_ray, cylinder_length2))
				{
					//cout << "Hit Cylinder\n";
					float pdf_light = pdf_light_cylinder;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = cylinder_light_direction.dot(rec.n);

					float inv_length2 = 1.0f / cylinder_length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light;// *100000 * inv_length2;

				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_SSS(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (Tr > 0)
				{
					//Good diffuse
					if (randf() < 0.8f)
					{					
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{

						float distance = 0.6f / cos_incident;

						if (isScatter(distance))
						{
							vec3 hit_point(new_ray.o + new_ray.d * distance);

							vec3 Transmitance(vec3(expf(-sigma_t.x * distance), expf(-sigma_t.y * distance), expf(-sigma_t.z * distance)));

							vec3 density = sigma_t * Transmitance;

							float atmos_pdf = (density.x + density.y + density.z) / 3;

							T *= Transmitance * sigma_s / atmos_pdf;

							vec3 direction_out;
							float p = HG_Sample(-new_ray.d, direction_out);

							new_ray = Ray(hit_point, direction_out);

							prev_pdf = p;

							specular = false;
						}
						else
						{
							new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
							continue;
						}
						

					}




				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					
					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();


					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);


						L += T * sun_color * directional_bsdf * scale_sun;

					}
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				

				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))

				if (length <= 1000 && hit_point.y < 2798.2f && !hit_anything_range_2(n, light_ray, length))

					//if (length < 2000 && !hit_anything_range_2(n, light_ray, length))

					//if(!hit_anything_range_alpha(n, light_ray, length))
				{
					//cout << "Not Hit\n";
					vec3 norm(0, -1, 0);

					float cos_light = -light_direction.dot(norm);

					//cos_light = abs(cos_light);

					//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

					//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}


				vec3 random_sphere_point;

				random_point_on_sphere(pic_light_radius, random_sphere_point);

				vec3 picture_light_point = picture_light_position + random_sphere_point;

				vec3 cylinder_light_direction(picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				if (cylinder_length <= 1000 && hit_point.y < 2437.5 && !hit_anything_range_2(n, cylinder_light_ray, cylinder_length2))
				{
					//cout << "Hit Cylinder\n";
					float pdf_light = pdf_light_cylinder;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = cylinder_light_direction.dot(rec.n);

					float inv_length2 = 1.0f / cylinder_length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light;// *100000 * inv_length2;

				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_Fall_Off_Mask(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (Tr > 0)
				{
					
					//float floatOut = acos(1.0 - (1 - (cos_incident * -1))) / (M_PI * 0.5);
					
					float floatOut = acosf(-cos_incident) * i2pi * 2.0f;

					
					//if (randf() < (1.0f - cos_incident) * 6.0f)

					if(randf() < floatOut)
					{						
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{
						
						vec3 direction_out = Reflect(new_ray.d, rec.n);

						float blur_radius = 0.4f;
						vec3 blur_point;

						random_point_on_sphere(blur_radius, blur_point);

						direction_out = (direction_out + blur_point).norm();

						float cos_o = direction_out.dot(rec.n);
						T *= powf(abs(cos_o), cos_incident);

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;

					}




				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					/*float x = directional_blur_radius * sin_phi * cos_theta;
					float y = directional_blur_radius * sin_phi * sin_theta;
					float z = directional_blur_radius * cos_phi;*/

					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					//HitRecord new_rec;

					/*if (hit_color_alpha(n, directional_lighting, new_rec, &color, Tr))
					{

					if (randf() <= Tr)
					{
					vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

					float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
					float sun_pdf = 1.0f;//there is only one way to sample this direction

					float mis = mis2(sun_pdf, directional_pdf * cont_prob);


					L += T * sun_color * directional_bsdf * mis;

					//L += T * sun_color * directional_bsdf;
					//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					//continue;
					}
					}*/
					///*else
					//vec3 new_color;
					//int new_mtl = trs[new_rec.ti].mtl;
					//if (!hit_color_alpha(n, directional_lighting, new_rec, &new_color, Tr))
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);


						L += T * sun_color * directional_bsdf * scale_sun;

					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				

				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))

				if (length <= 1000 && hit_point.y < 2798.2f && !hit_anything_range_2(n, light_ray, length))

					//if (length < 2000 && !hit_anything_range_2(n, light_ray, length))

					//if(!hit_anything_range_alpha(n, light_ray, length))
				{
					//cout << "Not Hit\n";
					vec3 norm(0, -1, 0);

					float cos_light = -light_direction.dot(norm);

					//cos_light = abs(cos_light);

					//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

					//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}


				vec3 random_sphere_point;

				random_point_on_sphere(pic_light_radius, random_sphere_point);

				vec3 picture_light_point = picture_light_position + random_sphere_point;

				vec3 cylinder_light_direction(picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				if (cylinder_length <= 1000 && hit_point.y < 2437.5 && !hit_anything_range_2(n, cylinder_light_ray, cylinder_length2))
				{
					//cout << "Hit Cylinder\n";
					float pdf_light = pdf_light_cylinder;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = cylinder_light_direction.dot(rec.n);

					float inv_length2 = 1.0f / cylinder_length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light;// *100000 * inv_length2;

				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_transparent_Mask(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (Tr > 0)
				{

					//float floatOut = acos(1.0 - (1 - (cos_incident * -1))) / (M_PI * 0.5);

					float floatOut = acosf(-cos_incident) * i2pi * 2.0f;


					//if (randf() < (1.0f - cos_incident) * 6.0f)

					if (randf() < floatOut)
					{
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{

						vec3 direction_out = Reflect(new_ray.d, rec.n);

						float blur_radius = 0.4f;
						vec3 blur_point;

						random_point_on_sphere(blur_radius, blur_point);

						direction_out = (direction_out + blur_point).norm();

						float cos_o = direction_out.dot(rec.n);
						T *= powf(abs(cos_o), cos_incident);

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;

					}




				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					/*float x = directional_blur_radius * sin_phi * cos_theta;
					float y = directional_blur_radius * sin_phi * sin_theta;
					float z = directional_blur_radius * cos_phi;*/

					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					//HitRecord new_rec;

					/*if (hit_color_alpha(n, directional_lighting, new_rec, &color, Tr))
					{

					if (randf() <= Tr)
					{
					vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

					float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
					float sun_pdf = 1.0f;//there is only one way to sample this direction

					float mis = mis2(sun_pdf, directional_pdf * cont_prob);


					L += T * sun_color * directional_bsdf * mis;

					//L += T * sun_color * directional_bsdf;
					//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					//continue;
					}
					}*/
					///*else
					//vec3 new_color;
					//int new_mtl = trs[new_rec.ti].mtl;
					//if (!hit_color_alpha(n, directional_lighting, new_rec, &new_color, Tr))
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);


						L += T * sun_color * directional_bsdf * scale_sun;

					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}



				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))

				if (length <= 1000 && hit_point.y < 2798.2f && !hit_anything_range_2(n, light_ray, length))

					//if (length < 2000 && !hit_anything_range_2(n, light_ray, length))

					//if(!hit_anything_range_alpha(n, light_ray, length))
				{
					//cout << "Not Hit\n";
					vec3 norm(0, -1, 0);

					float cos_light = -light_direction.dot(norm);

					//cos_light = abs(cos_light);

					//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

					//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}


				vec3 random_sphere_point;

				random_point_on_sphere(pic_light_radius, random_sphere_point);

				vec3 picture_light_point = picture_light_position + random_sphere_point;

				vec3 cylinder_light_direction(picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				if (cylinder_length <= 1000 && hit_point.y < 2437.5 && !hit_anything_range_2(n, cylinder_light_ray, cylinder_length2))
				{
					//cout << "Hit Cylinder\n";
					float pdf_light = pdf_light_cylinder;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = cylinder_light_direction.dot(rec.n);

					float inv_length2 = 1.0f / cylinder_length2;
					L += T * vec3(4, 1.6, 0) * bsdf_eval * mis * cos_mtl / pdf_light;// *100000 * inv_length2;

				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_circle_blur_final(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				//return color;
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (Tr > 0)
				{
					//Good diffuse
					if (randf() < 0.8f)
					{

						//float blur_radius = 0.02f;
						//vec3 blur_point;

						//random_point_on_sphere(blur_radius, blur_point);

						//vec3 hit_point(new_ray.o + new_ray.d * rec.t * 1.0002f);
						//vec3 direction_point = hit_point + new_ray.d + blur_point;
						//vec3 direction_out = (direction_point - hit_point).norm();


						//optimize

						//vec3 direction_out = (new_ray.d + blur_point).norm();

						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{

						/*vec3 direction_out = cosine_hemisphere_Sampling_with_n(rec.n);

						float cos_o = direction_out.dot(rec.n);


						float sin_o = sqrt14(1.0f - cos_o * cos_o);
						float horizon_scatter = powf(sin_o, 0.5);

						vec3 R(0.5f);

						vec3 weight = color * cos_o + horizon_scatter * R * cos_incident * ipi;

						//T *= color * direction_out.dot(rec.n);

						T *= weight;

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;*/

						vec3 direction_out = Reflect(new_ray.d, rec.n);

						float blur_radius = 0.4f;
						vec3 blur_point;

						random_point_on_sphere(blur_radius, blur_point);

						direction_out = (direction_out + blur_point).norm();

						float cos_o = direction_out.dot(rec.n);
						T *= powf(abs(cos_o), cos_incident);

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;

					}




				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light && mats[mtl].type != Glass_type || mats[mtl].type != ThinDielectric_type)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					/*float x = directional_blur_radius * sin_phi * cos_theta;
					float y = directional_blur_radius * sin_phi * sin_theta;
					float z = directional_blur_radius * cos_phi;*/

					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					//HitRecord new_rec;

					/*if (hit_color_alpha(n, directional_lighting, new_rec, &color, Tr))
					{

					if (randf() <= Tr)
					{
					vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

					float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
					float sun_pdf = 1.0f;//there is only one way to sample this direction

					float mis = mis2(sun_pdf, directional_pdf * cont_prob);


					L += T * sun_color * directional_bsdf * mis;

					//L += T * sun_color * directional_bsdf;
					//new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
					//continue;
					}
					}*/
					///*else
					//vec3 new_color;
					//int new_mtl = trs[new_rec.ti].mtl;
					//if (!hit_color_alpha(n, directional_lighting, new_rec, &new_color, Tr))
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

						//if(mats[mtl].type != Glass_type || mats[mtl].type != ThinDielectric_type)
						L += T * sun_color * directional_bsdf * scale_sun;
						//else
						//	L += T * vec3(0.4941f, 0.449f, 0.2196f) * directional_bsdf * scale_sun;

						//0.9882, 0.8980, 0.4392


					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_alpha(n, light_ray))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				//sum up light

				//int ind = randf() * round_light[0].size();

				//sample hemisphere
				/*vec3 sphere_norm(0, -1, 0);
				vec3 coord = uniform_hemisphere_Sampling(sphere_norm);

				vec3 sphere_point = local_sphere.u * coord.x + local_sphere.v * coord.y + local_sphere.w * coord.z;

				sphere_point = light_center[0] - vec3(0, 120, 0) + sphere_point * sphere_radius;

				//sphere_point *= sphere_radius;

				vec3 light_direction(sphere_point - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				float inv_length = 1.0f / length;

				light_direction *= inv_length;*/

				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))

				if (length <= 1000 && hit_point.y < 2406.7f && !hit_anything_range_2(n, light_ray, length))

					//if (length < 2000 && !hit_anything_range_2(n, light_ray, length))

					//if(!hit_anything_range_alpha(n, light_ray, length))
				{
					//cout << "Not Hit\n";
					vec3 norm(0, -1, 0);

					float cos_light = -light_direction.dot(norm);

					//cos_light = abs(cos_light);

					//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

					//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(1.0f, 1.0f, 0.2f) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}


				vec3 random_sphere_point;

				random_point_on_sphere(pic_light_radius, random_sphere_point);

				vec3 picture_light_point = picture_light_position +random_sphere_point;

				vec3 cylinder_light_direction(picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				vec3 vertical_norm(0, 1, 0);

				float cos_view = cylinder_light_direction.dot(vertical_norm);

				if (cylinder_length <= 1000 && hit_point.y < 2437.5 && !hit_anything_range_2(n, cylinder_light_ray, cylinder_length))
				{
					//cout << "Hit Cylinder\n";
					float pdf_light = 1.0f;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = cylinder_light_direction.dot(rec.n);

					//float inv_length2 = 1.0f / cylinder_length2;
					L += T * vec3(4.0f, 4.0f, 0.2f) * bsdf_eval * mis * cos_mtl / pdf_light;// *100000 * inv_length2;

				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	//chay tot
	vec3 __fastcall sibenik_tracing_directional_lighting_circle_blur_final_2(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				//return color;
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (Tr > 0)
				{
					//Good diffuse
					if (randf() < 0.8f)
					{						
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{

						vec3 direction_out = Reflect(new_ray.d, rec.n);

						float blur_radius = 0.4f;
						vec3 blur_point;

						random_point_on_sphere(blur_radius, blur_point);

						direction_out = (direction_out + blur_point).norm();

						float cos_o = direction_out.dot(rec.n);
						T *= powf(abs(cos_o), cos_incident);

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;

					}




				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light && mats[mtl].type != Glass_type || mats[mtl].type != ThinDielectric_type)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					/*float x = directional_blur_radius * sin_phi * cos_theta;
					float y = directional_blur_radius * sin_phi * sin_theta;
					float z = directional_blur_radius * cos_phi;*/

					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;
					
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

						//if(mats[mtl].type != Glass_type || mats[mtl].type != ThinDielectric_type)
						L += T * sun_color * directional_bsdf * scale_sun;
						//else
						//	L += T * vec3(0.4941f, 0.449f, 0.2196f) * directional_bsdf * scale_sun;

						//0.9882, 0.8980, 0.4392


					}//*/
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_range_alpha(n, light_ray, length))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				//sum up light

				//int ind = randf() * round_light[0].size();

				//sample hemisphere
				/*vec3 sphere_norm(0, -1, 0);
				vec3 coord = uniform_hemisphere_Sampling(sphere_norm);

				vec3 sphere_point = local_sphere.u * coord.x + local_sphere.v * coord.y + local_sphere.w * coord.z;

				sphere_point = light_center[0] - vec3(0, 120, 0) + sphere_point * sphere_radius;

				//sphere_point *= sphere_radius;

				vec3 light_direction(sphere_point - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				float inv_length = 1.0f / length;

				light_direction *= inv_length;*/

				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))

				if (length <= 1000 && hit_point.y < 2406.7f && !hit_anything_range_2(n, light_ray, length))

					//if (length < 2000 && !hit_anything_range_2(n, light_ray, length))

					//if(!hit_anything_range_alpha(n, light_ray, length))
				{
					//cout << "Not Hit\n";
					//vec3 norm(0, -1, 0);

					//float cos_light = -light_direction.dot(norm);

					//cos_light = abs(cos_light);

					//float pdf_light = length2 / (round_light_area[min_light][light_ind] * cos_light) * round_light_pdf[min_light];//cos_mtl

					//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(1.0f, 1.0f, 0.2f) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}

				//sphere
				vec3 random_sphere_point;
				random_point_on_sphere(pic_light_radius, random_sphere_point);
				vec3 picture_light_point = picture_light_position + random_sphere_point;

				//circle
				//float rx, rz;
				//random_point_on_circle(pic_light_radius, rx, rz);
				//vec3 picture_light_point = picture_light_position + vec3(rx, 0, rz);

				vec3 cylinder_light_direction = (picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				//vec3 vertical_norm(0, 1, 0);

				//float cos_view = cylinder_light_direction.dot(vertical_norm);

				//float cos_view = cylinder_light_direction.dot(rec.n);

				//vec3 inv_vertical_norm(0.0f, -1.0f, 0.0f);

				float cos_view = cylinder_light_direction.y;//-cylinder_light_direction.dot(inv_vertical_norm);

				//if (cylinder_length <= 2600.0f && cos_view >= 0.6f && hit_point.y <= 2290.875732f && !hit_anything_range_alpha(n, cylinder_light_ray, cylinder_length))
				//if (cylinder_length <= 2600.0f && mats[mtl].Matname == "Rectangle2106816060" && hit_point.y <= limit_y && !hit_anything_range_alpha(n, cylinder_light_ray, cylinder_length))
				//if(mats[mtl].Matname == "Rectangle2106816060" && hit_point.y <= limit_y)//2290.875732f)
				
				//if(cylinder_length <= 2000.0f && cos_view >= 0.6f && hit_point.y <= 2990.0f)//limit_y)
				//if (cylinder_length <= 2500.0f && cos_view >= 0.68f && hit_point.y <= 2420.0f)
			
				if (cylinder_length <= 2200.0f && cos_view >= 0.88f && hit_point.y <= 2480.0f)
				{
					//cout << "Hit Cylinder\n";
					float pdf_light = pdf_light_radius;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = cylinder_light_direction.dot(rec.n);

					float inv_length2 = 1.0f / cylinder_length2;
					//L += T * vec3(6.0f, 0.0f, 0.0f) * bsdf_eval * mis * cos_mtl;// / pdf_light;// *100000 * inv_length2;

					//255 165 0 orange
					//1 0.647 0
					//L += T * vec3(4.0f, 4.0f, 0.8f) * bsdf_eval * mis * cos_mtl * inv_length2 * 1000000.0f;

					//L += T * 4.0f * vec3(1.0f, 0.647f, 0.0f) * bsdf_eval;// *cos_mtl;// *inv_length2 * 1000000.0f;

					//L += T * 6000000.0f * vec3(1.0f, 0.513f, 0.0f) * cos_mtl * bsdf_eval * inv_length2;

					float real_length = length * 0.004f;

					float inv_real_length2 = 1.0f / (real_length * real_length);

					L += T * 16 * vec3(1.0f, 0.513f, 0.0f) * inv_real_length2;// *cos_mtl;

					//float length = 

					//L += T * vec3(1.0f, 1.0f, 0.2f) * bsdf_eval;// *mis / pdf_light;// l * cos_mtl * inv_length2 * 1000000.0f;
				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_circle_blur_final_version(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				//return color;
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

			
				if (Tr > 0)
				{
					
					if (randf() < 0.8f)
					{
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{

						vec3 direction_out = Reflect(new_ray.d, rec.n);

						float blur_radius = 0.4f;
						vec3 blur_point;

						random_point_on_sphere(blur_radius, blur_point);

						direction_out = (direction_out + blur_point).norm();

						float cos_o = direction_out.dot(rec.n);
						T *= powf(abs(cos_o), cos_incident);

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;

					}




				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				
				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light && mats[mtl].type != Glass_type || mats[mtl].type != ThinDielectric_type)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

				
					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					
					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;

					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);
						
						L += T * sun_color * directional_bsdf * scale_sun;					
					}
				}
				if (use_area_light)
				{
					int li = Ls.sample();
					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);

					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);


					float cos_mtl = light_direction.dot(rec.n);


					if (cos_mtl > 0.0f && !hit_anything_range_alpha(n, light_ray, length))
					{
						float cos_light = abs(-light_direction.dot(Ls.normal[li]));


						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																									 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * mats[mtl].Tf * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

					}
				}

				
				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				
				if (length <= 1000 && hit_point.y < 2406.7f && !hit_anything_range_2(n, light_ray, length))
				{
					
					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(1.0f, 1.0f, 0.2f) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}

				//sphere
				vec3 random_sphere_point;
				random_point_on_sphere(pic_light_radius, random_sphere_point);
				vec3 picture_light_point = picture_light_position + random_sphere_point;

				//circle
				//float rx, rz;
				//random_point_on_circle(pic_light_radius, rx, rz);
				//vec3 picture_light_point = picture_light_position + vec3(rx, 0, rz);

				vec3 cylinder_light_direction = (picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				float cos_view = cylinder_light_direction.y;
				float cos_mtl = cylinder_light_direction.dot(rec.n);

				if (cylinder_length <= 2400.0f && cos_view >= 0.6f && cos_mtl >= 0.0f && hit_point.y <= 2420.0f)
				{
					
					float pdf_light = pdf_light_radius;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					

					float inv_length2 = 1.0f / cylinder_length2;
					
					//float real_length = length * 0.004f;

					//float inv_real_length2 = 1.0f / (real_length * real_length);

					//L += T * 16 * vec3(1.0f, 0.513f, 0.0f) * inv_real_length2;

					//L += T * 16 * vec3(1.0f, 0.513f, 0.0f) * bsdf_eval * cos_mtl;

					L += T * 4 * vec3(1.0f, 0.513f, 0.0f) * bsdf_eval;
				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall cone_light(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color(1.0f, 1.0f, 1.0f);

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float sigma = 0.1f;

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				//return color;
				int mtl = trs[rec.ti].mtl;


				if (mats[mtl].type == Back_Ground_type)
				{
					//L += 19.5 * T * color;
					L += 29.5 * T * color;
					return L;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;

				cos_incident = abs(cos_incident);

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (Tr > 0)
				{
					//Good diffuse
					if (randf() < 0.8f)
					{
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, new_ray.d);
						continue;
					}
					else
					{

						vec3 direction_out = Reflect(new_ray.d, rec.n);

						float blur_radius = 0.4f;
						vec3 blur_point;

						random_point_on_sphere(blur_radius, blur_point);

						direction_out = (direction_out + blur_point).norm();

						float cos_o = direction_out.dot(rec.n);
						T *= powf(abs(cos_o), cos_incident);

						new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, direction_out);
						continue;

					}




				}


				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{

					vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

					if (i == 0 && mats[mtl].type == Light_Diffuse_type)
					{


						L += T * Ke;// *6;//(1.0f + expf(scale_dist));

						return L;
					}

					if (specular || mats[mtl].type == Light_Diffuse_type)
					{

						L += T * Ke;

						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						float length2 = light_direction.length2();

						light_direction.normalize();


						float cos_light = -light_direction.dot(rec.n);

						float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

						float mis_weight = mis2(prev_pdf, pdf_light);

						L += T * Ke * mis_weight;


						return L;

					}

				}

				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				if (use_directional_light && mats[mtl].type != Glass_type || mats[mtl].type != ThinDielectric_type)
				{
					vec3 sun_direction = inifinite_light_direction;
					//if (blur_directional_light)
					//{

					vec3 blur_center = hit_point + sun_direction;

					float u = randf();
					float v = randf();

					float theta = tau * u;
					float phi = acosf(2.0f * v - 1.0f);

					float sin_phi, cos_phi;

					FTA::sincos(phi, &sin_phi, &cos_phi);

					float sin_theta, cos_theta;

					FTA::sincos(theta, &sin_theta, &cos_theta);

					float blur_radius = directional_blur_radius, scale_sun = 1.0f;

					if (randf() < 0.1f)
					{
						blur_radius = small_blur_radius;
						scale_sun = 6.0f;
					}

					
					float x = blur_radius * sin_phi * cos_theta;
					float y = blur_radius * sin_phi * sin_theta;
					float z = blur_radius * cos_phi;

					vec3 sample_point = blur_center + vec3(x, y, z);
					sun_direction = (sample_point - hit_point).norm();

					//}

					Ray directional_lighting(hit_point, sun_direction);

					float Tr = 0.0f;

					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

						//if(mats[mtl].type != Glass_type || mats[mtl].type != ThinDielectric_type)
						L += T * sun_color * directional_bsdf * scale_sun;
						

					}
				}
				
				float x, z;

				int light_ind = 0;

				float max_length2 = (light_center[0] - hit_point).length2();

				for (int l = 1; l < 4; ++l)
				{
					float current_length2 = (light_center[l] - hit_point).length2();

					if (current_length2 < max_length2)
					{
						light_ind = l;
						max_length2 = current_length2;
					}
				}

				random_point_on_circle(light_radius, x, z);

				vec3 light_position = light_center[light_ind] + vec3(x, 0, z);

				vec3 light_direction(light_position - hit_point);

				float length2 = light_direction.length2();

				float length = sqrt14(length2);

				light_direction = light_direction.norm();

				Ray light_ray(hit_point, light_direction);

				//length += 40.0f;
				//if (!hit_anything_alpha(n, light_ray))
				//if(!hit_anything_range_alpha(n, light_ray, length))

				if (length <= 1000 && hit_point.y < 2406.7f && !hit_anything_range_2(n, light_ray, length))

				{
					
					float pdf_light = pdf_light_radius;//0.25f * ipi;//1;
					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = light_direction.dot(rec.n);

					float inv_length2 = 1.0f / length2;
					L += T * vec3(1.0f, 1.0f, 0.2f) * bsdf_eval * mis * cos_mtl / pdf_light * 100000 * inv_length2;// *mis / pdf_light;// * 0.1f
				}

				//sphere
				//vec3 random_sphere_point;
				//random_point_on_sphere(pic_light_radius, random_sphere_point);
				//vec3 picture_light_point = picture_light_position + random_sphere_point;

				//2603.1689f

				//if (picture_light_point.y >= 2603.1689f)
				//	picture_light_point.y = 2.0f * 2603.1689f - picture_light_point.y;

				//circle
				float rx, rz;
				random_point_on_circle(pic_light_radius, rx, rz);
				vec3 picture_light_point = picture_light_position + vec3(rx, 0, rz);

				vec3 cylinder_light_direction = (picture_light_point - hit_point);

				float cylinder_length2 = cylinder_light_direction.length2();

				float cylinder_length = sqrt14(length2);

				cylinder_light_direction = cylinder_light_direction.norm();

				Ray cylinder_light_ray(hit_point, cylinder_light_direction);

				float upper_limit = 2603.1689f - 5.0f;

				vec3 inv_normal(0.0f, -1.0f, 0.0f);

				float cos_view = -cylinder_light_direction.dot(inv_normal);

				//picture_light_point.y <= upper_limit

				if (cylinder_length <= 2200.0f && cos_view >= 0.7f && hit_point.y <= 2460.0f)////3000.0f)
				{
					//cout << "Hit Cylinder\n";
					float pdf_light = pdf_light_radius;

					vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, cylinder_light_direction, rec.n, original_normal, color);
					float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, cylinder_light_direction, rec.n, original_normal);


					float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

					float cos_mtl = cylinder_light_direction.dot(rec.n);

					//float inv_length2 = 1.0f / cylinder_length2;
					//L += T * vec3(6.0f, 0.0f, 0.0f) * bsdf_eval * mis * cos_mtl;// / pdf_light;// *100000 * inv_length2;

					//255 165 0 orange
					//1 0.647 0
					//L += T * vec3(4.0f, 4.0f, 0.8f) * bsdf_eval * mis * cos_mtl * inv_length2 * 1000000.0f;

					//L += T * 4.0f * vec3(1.0f, 0.647f, 0.0f) * bsdf_eval;// *cos_mtl;// *inv_length2 * 1000000.0f;

					//L += T * 6000000.0f * vec3(1.0f, 0.513f, 0.0f) * cos_mtl * bsdf_eval * inv_length2;

					float real_length = cylinder_length * 0.001f;

					float inv_real_length2 = 1.0f / (real_length * real_length);

					L += T * 4 * vec3(1.0f, 0.513f, 0.0f) * bsdf_eval;// *inv_real_length2;// *cos_mtl;

																			 //float length = 

																			 //L += T * vec3(1.0f, 1.0f, 0.2f) * bsdf_eval;// *mis / pdf_light;// l * cos_mtl * inv_length2 * 1000000.0f;
				}

				Sampling_Struct sampling;
				bool isReflect = true;


				bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

				if (!b)
					return vec3(0.0f);

				rec.n = isReflect ? rec.n : -rec.n;

				vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

				//if (sampling.ray_shift == 0.0f)
				///	new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;
				//else
				//	new_hit_point = new_ray.o + new_ray.d * rec.t + sampling.direction * sampling.ray_shift;

				prev_pdf = sampling.pdf;

				new_ray = Ray(new_hit_point, sampling.direction);

				specular = sampling.isSpecular;


				T *= sampling.weight;

				if (i > 6)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//T *= color;
				if (compute_sky_color)
				{


					//L += T * 19.5;// *12.0f;

					L += T * 29.5;

					return L;
				}
				else
				{
					float u, v;

					Env.DirectionToUV(new_ray.d, u, v);

					u += rotation_u;

					vec3 env_value = Env.UVToValue(u, v);

					L += T * env_value;
					return L;
				}


			}
		}
		return L;
	}

	vec3 __fastcall sibenik_tracing_directional_lighting_sun_mirror(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		//vec3 Ke = 

		for (int i = 0; i < 60; ++i)
		{
			float Tr = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, Tr))
			{
				if (randf() <= Tr)
				{
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;
				}

				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;


				if (cos_incident > 0)
				{

					rec.n = -rec.n;
				}


				int mtl = trs[rec.ti].mtl;

				if ((mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type))
				{
					if (use_area_light)
					{
						vec3 Ke = mats[mtl].isLight ? mats[mtl].Ke : mats[mtl].Kd;

						if (specular)
						{

							L += T * Ke;// *0.1f;

							return L;
						}
						else
						{
							int li = light_map[rec.ti];

							vec3 light_direction(new_ray.d * rec.t);

							float length2 = light_direction.length2();

							light_direction.normalize();


							float cos_light = -light_direction.dot(rec.n);

							float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

							float mis_weight = mis2(prev_pdf, pdf_light);

							//Good
							//L += T * mats[mtl].Ke * mis_weight;

							L += T * Ke * mis_weight;// *0.1f;


							return L;

						}
					}
					else
					{
						//return vec3(0.0f);
						L += T;// *0.1f;// *0.0000001f;
						return L;
					}
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				//float cont_prob = max(min(T.maxc(), 1.0f), 0.01f);

				float cont_prob = min(T.maxc(), 1.0f);

				//sample direct lighting				

				if (use_directional_light)
				{
					vec3 sun_direction = inifinite_light_direction;
					if (blur_directional_light)
					{
						//from hit point move 1 unit in light direction
						//create a sphere of radius equal to blur radius
						//sample random point on sphere, minus these point with hitpoint
						//this will be the new sun_direction

						vec3 blur_center = hit_point + sun_direction;

						float u = randf();
						float v = randf();

						float theta = tau * u;
						float phi = acosf(2.0f * v - 1.0f);

						float sin_phi, cos_phi;

						FTA::sincos(phi, &sin_phi, &cos_phi);

						float sin_theta, cos_theta;

						FTA::sincos(theta, &sin_theta, &cos_theta);

						float x = directional_blur_radius * sin_phi * cos_theta;
						float y = directional_blur_radius * sin_phi * sin_theta;
						float z = directional_blur_radius * cos_phi;

						//original

						vec3 sample_point = blur_center + vec3(x, y, z);
						sun_direction = (sample_point - hit_point).norm();

						//wrong
						//sun_direction = vec3(x, y, z);
						//sun_direction.normalize();

					}
					//Ray directional_lighting(hit_point, inifinite_light_direction);
					Ray directional_lighting(hit_point, sun_direction);

					//if (!hit_anything_range(n, directional_lighting))
					if (!hit_anything_alpha(n, directional_lighting))
					{
						vec3 directional_bsdf = bsdfs[mtl]->eval(new_ray.d, sun_direction, rec.n, original_normal, color);

						//MIS
						//float directional_pdf = bsdfs[mtl]->pdf(new_ray.d, inifinite_light_direction, rec.n, original_normal);
						//float sun_pdf = 1.0f;//there is only one way to sample this direction

						//float mis = mis2(sun_pdf, directional_pdf * cont_prob);


						//L += T * sun_color * directional_bsdf * mis;

						//No MIS
						//use a different light color for specular material

						vec3 light_color = !bsdfs[mtl]->isSpecular() ? sun_color : vec3(0.03, 0.02, 0.01);
						L += T * light_color * directional_bsdf;

					}
				}

				if (!bsdfs[mtl]->isSpecular())
				{

					if (use_area_light)
					{
						int li = Ls.sample();
						vec3 light_point(Ls.sample_light(li));

						vec3 light_direction(light_point - hit_point);

						float length = light_direction.length();

						float ilength = 1.0f / length;
						light_direction *= ilength;

						Ray light_ray(hit_point, light_direction);


						float cos_mtl = light_direction.dot(rec.n);


						if (cos_mtl > 0.0f && !hit_anything_range(n, light_ray))
						{
							float cos_light = abs(-light_direction.dot(Ls.normal[li]));


							float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

																										 //float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

							vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
							float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


							float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

							//cout << li << "\n";

							//cout << Ls.Ke.size() << "\n";
							L += T * Ls.Ke[li] * bsdf_eval * mis / pdf_light;// * 0.1f

						}
					}

					Sampling_Struct sampling;
					bool isReflect;


					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;


					prev_pdf = sampling.pdf;

					new_ray = Ray(new_hit_point, sampling.direction);

					specular = sampling.isSpecular;


					T *= sampling.weight;
				}
				else
				{

					Sampling_Struct sampling;
					bool isReflect;

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;

					T *= sampling.weight;
				}
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{

				//L += T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

				//vec3 clamp_L = clamp_vector(L, 0.0f, 1.0f);
				//return clamp_L;

				//return L;

				//if use enviroment compute this value

				//compute sky color will create better result

				//this is from Ray Tracing In One Weekend


				if (i == 0)
				{
					if (compute_sky_color)
					{
						vec3 sky_direction = new_ray.d;

						float t = 0.5f * (sky_direction.y + 1.0f);

						//vec3 sky_color = (1.0f - t) * vec3(1.0f) + t * vec3(0.5f, 0.7f, 4.0f);

						vec3 sky_color = (1.0f - t) * vec3(1.0f) + t * vec3(0.2f, 0.4f, 6.0f);

						L += T * sky_color;

						return L;
					}
					else
					{
						float u, v;

						Env.DirectionToUV(new_ray.d, u, v);

						u += rotation_u;

						vec3 env_value = Env.UVToValue(u, v);

						L += T * env_value;
						return L;

						/*//if (specular)
						//{
						//	L += T * env_value;
						//	return L;
						//}
						//else
						//{
						float env_pdf = Env.UVToPdf(u, v);

						float mis = mis2(prev_pdf, env_pdf);

						L += T * env_value * mis;

						return L;*/
					}
				}
				else
				{
					L += T;// *0.1f;// *0.00001f;
					return L;
				}
				/*if (specular)
				{

				L += T * sky_color;

				return L;
				}
				else
				{
				int li = light_map[rec.ti];

				vec3 light_direction(new_ray.d * rec.t);

				float length2 = light_direction.length2();

				light_direction.normalize();


				float cos_light = -light_direction.dot(rec.n);

				float pdf_light = length2 / (Ls.area[li] * cos_light) * Ls.pdf[li];

				float mis_weight = mis2(prev_pdf, pdf_light);

				//Good
				//L += T * mats[mtl].Ke * mis_weight;

				L += T * sky_color * mis_weight;

				return L;

				}*/

			}
		}
		//return L;
		return L;
	}

	vec3 __fastcall sibenik_tracing_point_light_only(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		vec3 Ke = delta_light_Ke;

		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_alpha_pointer(n, new_ray, rec, &color, alpha))
			float alpha = 0.0f;
			if (hit_color_alpha(n, new_ray, rec, &color, alpha))
			{
				if (randf() < alpha)
				{
					int mtl = trs[rec.ti].mtl;
					T *= mats[mtl].Tf;
					new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.000002f, new_ray.d);
					continue;
				}


				float cos_incident = new_ray.d.dot(rec.n);

				vec3 original_normal = rec.n;

				if (cos_incident > 0)
					rec.n = -rec.n;


				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].isLight || mats[mtl].type == Light_Diffuse_type)
				{
					return vec3(0.0f);
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				if (!bsdfs[mtl]->isSpecular())
				{

					vec3 light_point = delta_light_position;

					vec3 light_direction(light_point - hit_point);



					float length = light_direction.length();

					float ilength = 1.0f / length;
					light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);

					float cos_mtl = light_direction.dot(rec.n);

					//length *= 1.11f;

					if (cos_mtl > 0.0f && !hit_anything_range_2(n, light_ray, length))
					{
						/*float cos_light = -light_direction.dot(Ls.normal[li]);

						float pdf_light = (length * length) / (Ls.area[li] * cos_light) * Ls.pdf[li];//cos_mtl

						//float pdf_light = (length * length) / (Ls.area[li] * cos_mtl) * Ls.pdf[li];

						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);


						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * Ls.Ke[li] * bsdf_eval * mis / pdf_light;
						*/

						//MIS version
						/*
						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);
						float bsdf_pdf = bsdfs[mtl]->pdf(new_ray.d, light_direction, rec.n, original_normal);

						float pdf_light = 1.0f;

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);

						L += T * Ke * bsdf_eval * mis / pdf_light;
						*/



						vec3 bsdf_eval = bsdfs[mtl]->eval(new_ray.d, light_direction, rec.n, original_normal, color);


						//Fall off
						float length2 = length * length;
						L += T * mats[mtl].Tf * Ke * bsdf_eval / length2;

						//Exponential
						//L += T * Ke * bsdf_eval * expf(-length * 0.5f);

						//if (i == 0)
						//	return L;
					}


					Sampling_Struct sampling;
					bool isReflect;


					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					rec.n = isReflect ? rec.n : -rec.n;

					vec3 new_hit_point = new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f;

					prev_pdf = sampling.pdf;

					new_ray = Ray(new_hit_point, sampling.direction);

					specular = sampling.isSpecular;

					T *= sampling.weight;
				}
				else
				{

					Sampling_Struct sampling;
					bool isReflect;

					bool b = bsdfs[mtl]->sample(new_ray.d, rec.n, original_normal, sampling, isReflect, color);

					if (!b)
						return vec3(0.0f);

					if (mats[mtl].type == Transparent_type)
					{
						new_ray = Ray(new_ray.o + new_ray.d * rec.t * 1.0002f, sampling.direction);
						continue;
					}

					rec.n = isReflect ? rec.n : -rec.n;

					new_ray = Ray(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f, sampling.direction);
					prev_pdf = sampling.pdf;
					specular = true;


					T *= sampling.weight;
				}
				if (i > 2)
				{
					prev_pdf *= cont_prob;
					if (cont_prob < randf())
					{
						//L += T;
						//return L;
						break;
					}
					float inv_cont_prob = 1.0f / cont_prob;
					T *= inv_cont_prob;

				}
			}
			else
			{
				//L += T;// *vec3(4.575004577636719, 3.5907576084136963, 1.5497708320617676);

				return L;
			}
		}
		//cheat because back face cause many dark place
		//L += T;
		return L;
	}

	vec3 Uniform_Hemisphere_Sampling()
	{
		float cos_theta = randf();

		float sin_theta = sqrt14(1.0f - cos_theta * cos_theta);

		float phi = 2.0f * pi * randf();

		float x = cosf(phi) * sin_theta;

		float y = sinf(phi) * sin_theta;

		float z = cos_theta;

		return vec3(x, y, z);
	}

	bool __fastcall hit_anything_range_test(node*& n, Ray& r, const float& d, const int& first)
	{
		node* stack[26];



		stack[0] = n->nodes[0];
		stack[1] = n->nodes[1];

		int si = 1;

		while (si >= 0)
		{
			auto top(stack[si--]);

			float th;

			if (top->box.hit_shadow_th(r, th) && th > 0.0f)
			{
				if (!top->leaf)
				{
					stack[si + 1] = top->nodes[1 - first];
					stack[si + 2] = top->nodes[first];

					si += 2;
				}
				else
				{

					auto start = top->start, end = start + top->range;

					for (auto i = start; i < end; ++i)
					{
						int mtl = trs[i].mtl;
						//if (trs[i].hit_anything(r, d) && !mats[mtl].isLight)
						//if (trs[i].hit_anything(r) && !mats[mtl].isLight)

						int ind0 = trs[i].ind_v0;
						int ind1 = trs[i].ind_v1;
						int ind2 = trs[i].ind_v2;

						int v0 = indices[ind0].v;
						int v1 = indices[ind1].v;
						int v2 = indices[ind2].v;

						vec3 p0 = v[v0];
						vec3 p1 = v[v1];
						vec3 p2 = v[v2];

						vec3 e1(p1 - p0);
						vec3 e2(p2 - p0);

						if (hit_anything(r, p0, e1, e2))
							return true;
					}
				}
			}
			else
				continue;
		}
		return false;
	}

	vec3 __fastcall ray_casting(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		//float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float alpha;

		if (hit_color_alpha(n, new_ray, rec, &color, alpha))
		{
			//return rec.n;
			//vec3 v = 0.5f * rec.n + vec3(0.5f);
			//return v;
			//int mtl = trs[rec.ti].mtl;
			//if (mats[mtl].Matname == "metal_gold")
			return color;
			//return vec3(0.0f);

		}
		return vec3(1.0f);
	}

	vec3 __fastcall ray_casting_normal(node*& n, Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		//float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float alpha;

		if (hit_color_alpha(n, new_ray, rec, &color, alpha))
		{
			return color;
			//return 0.5f * (rec.n + vec3(1.0f));
		}
		return vec3(0.0f);
	}

	void __fastcall ray_casting_color_normal(node*& n, Ray& r, vec3& col, vec3& norm)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		//float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;

		float alpha = 1.0f;

		if (hit_color_alpha(n, new_ray, rec, &color, alpha))
		{
			col = color;
			norm = rec.n;
		}
		else
		{
			col = vec3(1.0f);
			norm = vec3(0.0f);
		}
	}


};
#endif // !_VEC3_H_


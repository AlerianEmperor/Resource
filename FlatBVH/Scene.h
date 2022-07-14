#ifndef _SCENE_H_
#define _SCENE_H_
#define _CRT_SECURE_NO_WARNINGS
#include "AreaLight.h"
//#include "EnvLight.h"
#include "Triangle.h"
#include "material.h"
#include "onb.h"
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <string>
#include <fstream>
//#include <omp.h>

#define max_depth 5
#define max_bin 60
#define imax_bin 1 / max_bin
#define ct 1//traversal cost, or ray bbox intersection cost
#define ci 4//intersect triangle cost

#define alpha 0.95
#define minus_alpha 1.0f - alpha

//ct  ci   leaf
//1   4    2      67.3
//1   2    2      68.3
//1   4    4      68.7
//1   8    4      68.8
using namespace std;

#define pi 3.1415926f
#define tau 6.2831853f
//#define SWAP(x, y) do { typeof(x) SWAP = x; x = y; y = SWAP; } while (0)
#define SWAP(a, b) swap_detail(&(a), &(b), (char[(sizeof(a) == sizeof(b)) ? (ptrdiff_t)sizeof(a) : -1]){0}, sizeof(a))




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

	node() {}
	node(int s, int r) : start(s), range(r) {}
};

struct FlatNode
{
	BBox box;

	unsigned start;
	int16_t range;
	int8_t axis;
	bool leaf;

	//unsigned int start : 25;
	//unsigned int range : 3;
	//unsigned int axis  : 2;
	//unsigned int leaf  : 1;
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
	//delete n;

	//n = NULL;
	//free(n);
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

void get_face_index(char*& t, const int& vs, const int& vts, vector<Face>& faces)
{
	string s = t;
	int length = s.find_last_of("0123456789");
	s = s.substr(0, length + 1);

	int sign = 1;
	int count = 0;
	vector<int> index;
	//vector<int> indices;
	int face_type = Single_Line;
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

	if (size == 9)
	{
		for (int i = 0; i < 9; i += 3)
		{
			//int ind0 = index[i];
			//int ind1 = index[i + 1];
			fixIndex(index[i], vs);
			fixIndex(index[i + 1], vts);

			//cout << index[i] << "\n";
			//indices.push_back(index[i]);
			//indices.push_back(index[i + 1]);
		}
		Face f(index[0], index[1], index[3], index[4], index[6], index[7]);
		faces.emplace_back(f);
	}
	else if (size == 12)
	{
		if (face_type == Single_Line)
		{
			for (int i = 0; i < 12; i += 3)
			{
				fixIndex(index[i], vs);
				fixIndex(index[i + 1], vts);
				//cout << index[i] << "\n";
			}
			//0 1 2  3 4 5  6 7 8  9 10 11 
			Face f1(index[0], index[1], index[3], index[4], index[6], index[7]);
			Face f2(index[0], index[1], index[6], index[7], index[9], index[10]);

			faces.emplace_back(f1);
			faces.emplace_back(f2);
		}
		if (face_type == Double_Line)
		{
			for (int i = 0; i < 12; i += 2)
			{
				fixIndex(index[i], vs);
				fixIndex(index[i + 1], vts);
			}

			Face f1(index[0], index[1], index[2], index[3], index[4], index[5]);
			Face f2(index[0], index[1], index[4], index[5], index[6], index[7]);
			Face f3(index[0], index[1], index[6], index[7], index[8], index[9]);
			Face f4(index[0], index[1], index[8], index[9], index[10], index[11]);

			faces.emplace_back(f1);
			faces.emplace_back(f2);
			faces.emplace_back(f3);
			faces.emplace_back(f4);

			//          01      10  11
			//
			//     23           89
			//
			//        45     67

		}
	}
	else if (size == 6)
	{
		//cout << "a";
		for (int i = 0; i < 6; i += 2)
		{
			fixIndex(index[i], vs);
			fixIndex(index[i + 1], vts);
			//indices.emplace_back(index[i]);
			//indices.emplace_back(index[i + 1]);
		}

		Face f(index[0], index[1], index[2], index[3], index[4], index[5]);
		faces.emplace_back(f);
	}
	else if (size == 8)
	{
		for (int i = 0; i < 8; i += 2)
		{
			fixIndex(index[i], vs);
			fixIndex(index[i + 1], vts);
		}

		Face f1(index[0], index[1], index[2], index[3], index[4], index[5]);
		Face f2(index[0], index[1], index[4], index[5], index[6], index[7]);

		faces.emplace_back(f1);
		faces.emplace_back(f2);

		//indices.emplace_back(index[0]);
		//indices.emplace_back(index[1]);

		//indices.emplace_back(index[2]);
		//indices.emplace_back(index[3]);

		//indices.emplace_back(index[4]);
		//indices.emplace_back(index[5]);

		//indices.emplace_back(index[0]);
		//indices.emplace_back(index[1]);

		//indices.emplace_back(index[4]);
		//indices.emplace_back(index[5]);

		//indices.emplace_back(index[6]);
		//indices.emplace_back(index[7]);
	}
	else if (size == 10)
	{
		for (int i = 0; i < 10; i += 2)
		{
			fixIndex(index[i], vs);
			fixIndex(index[i + 1], vts);
		}
		Face f1(index[0], index[1], index[2], index[3], index[4], index[5]);
		Face f2(index[0], index[1], index[4], index[5], index[6], index[7]);
		Face f3(index[0], index[1], index[6], index[7], index[8], index[9]);

		faces.emplace_back(f1);
		faces.emplace_back(f2);
		faces.emplace_back(f3);
	}

}

struct Scene
{
	Scene() {}
	Scene(const string& filename, const string& enviroment_path = "", const float& rotation = 0)
	{
		ReadObj(filename);	
	}
	
	vector<FlatNode> flat_nodes;

	vector<Triangle> trs;
	vector<Material> mats;
	vector<vec3> vt;
	vector<Texture> texs;

	AreaLight Ls;

	unordered_map<int, int> light_map;

	vector<vec3> v;
	vector<vec3> vn;
	
	int leaf = 4;
	bool use_eviroment;

	string direction;

	void delete_before_render()
	{
		vector<vec3>().swap(v);
		vector<vec3>().swap(vn);
	}

	void delete_struct()
	{
		vector<vec3>().swap(v);
		vector<vec3>().swap(vt);
		vector<vec3>().swap(vn);
		//vector<Face>().swap(fs);
		vector<Triangle>().swap(trs);
		vector<Material>().swap(mats);
		for (auto& v : texs)
		{
			vector<float>().swap(v.c);
			//stbi_image_free(v.c);
		}

		vector<Texture>().swap(texs);

		vector<FlatNode>().swap(flat_nodes);
		//vector<uint32_t>().swap(ti);

	}

	void ReadMtl(const string& filename, unordered_map<string, int>& mtlMap)//mtlMap chi dung cho giai doan doc file, ender thi khong can)
	{
		cout << "texture path:" << direction + filename << "\n";
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
					mats[countMaterial].type = Mirror;
				if (strncmp(type, "Glass", 5) == 0)
					mats[countMaterial].type = Glass;
				if (strncmp(type, "Conductor", 9) == 0)
				{
					mats[countMaterial].type = Conductor;
					t += 9;
					char complex_ior_name[64];
					sscanf_s("%s", complex_ior_name);

					string ior_name = complex_ior_name;

					for (int i = 0; i < 40; ++i)
					{
						if (ior_name == Ior_List[i].name)
						{
							mats[countMaterial].complex_ior = Ior_List[i];
							break;
						}
					}
				}
			}
			else if (t[0] == 'N')
			{
				if (t[1] == 'i')
				{
					sscanf_s(t += 3, "%f", &mats[countMaterial].Ni);
					float Ni = mats[countMaterial].Ni;
					float R0 = (Ni - 1) / (Ni + 1);
					mats[countMaterial].R0 = R0 * R0;
				}
				else if (t[1] == 's')
					sscanf_s(t += 3, "%f", &mats[countMaterial].Ns);
			}

			else if (t[0] == 'd')
				sscanf_s(t += 2, "%d", &mats[countMaterial].d);
			else if (t[0] == 'T')
			{
				if (t[1] == 'r')
					sscanf_s(t += 3, "%f", &mats[countMaterial].Tr);
				else if (t[1] == 'f')
					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Tf.x, &mats[countMaterial].Tf.y, &mats[countMaterial].Tf.z);
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
				if (t[1] == 'a')
					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Ka.x, &mats[countMaterial].Ka.y, &mats[countMaterial].Ka.z);
				else if (t[1] == 'd')
				{
					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Kd.x, &mats[countMaterial].Kd.y, &mats[countMaterial].Kd.z);
					//mats[countMaterial].type = Diffuse;
				}
				else if (t[1] == 's')
				{
					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Ks.x, &mats[countMaterial].Ks.y, &mats[countMaterial].Ks.z);
					//mats[countMaterial].type = Diffuse;
				}
				else if (t[1] == 'e')
				{
					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].Ke.x, &mats[countMaterial].Ke.y, &mats[countMaterial].Ke.z);
					if (mats[countMaterial].Ke.maxc() > 0)
					{
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

					int w, h, n;
					//tex.c = stbi_load((direction + realname).c_str(), &w, &h, &n, 3);

					unsigned char* c = stbi_load((direction + realname).c_str(), &w, &h, &n, 3);

					//float* value = new float[3 * w * h];
					tex.c.resize(3 * w * h);

					//powf(c[i] / 255.0f, 2.2f), powf(c[i + 1] / 255.0f, 2.2f), powf(c[i + 2] / 255.0f, 2.2f)

					for (int e = 0; e < 3 * w * h; ++e)
					{
						tex.c[e] = powf(c[e] / 255.0f, 2.2f);
					}
					
					//vector<float>().swap(tex.c);

					//tex.c = value;

					free(c);
					//free(value);

					tex.w = w;
					tex.h = h;
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
				mats[countMaterial].use_alpha = true;
				sscanf_s(t += 6, "%s", texturename);
				string realname;

				getfilename(string(texturename), realname);

				int w2, h2, n;
				//stbi_grey_alpha fail
				//3 fail
				//0 fail
				//stbi_grey fail
				//4 

				//tex.a = stbi_load((direction + realname).c_str(), &w2, &h2, &n, 3);
				//tex.w2 = w2;
				//tex.h2 = h2;
				texs[countTexture] = tex;
			}
		}
	}

	void ReadObj(const string& p)
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
		int li = -1;
		bool isLight = false;
		int count = 0;
		unordered_map<string, int> mtlMap;
		bool read_mtl = false;

		int count_light = 0;
		int num_v = 0;
		int num_vt = 0;
		int num_vn = 0;

		while (f.getline(line, 1024))
		{
			char* t = line;
			int prev_space = strspn(t, " \t");
			t += prev_space;

			if (strncmp(t, "v", 1) == 0)
			{
				t += 1;
				float x, y, z;
				if (strncmp(t, " ", 1) == 0)
				{
					int post_space = strspn(t, " \t");
					t += post_space;
					sscanf_s(t, "%f %f %f", &x, &y, &z);
					v.emplace_back(vec3(x, y, z));
					++num_v;
				}
				else if (strncmp(t, "t", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					sscanf_s(t, "%f %f", &x, &y);
					vt.emplace_back(vec3(x, y, 0.0f));
					++num_vt;
				}
				else if (strncmp(t, "n", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					sscanf_s(t, "%f %f %f", &x, &y, &z);
					vn.emplace_back(vec3(x, y, z));
					++num_vn;
				}
			}
			else if (t[0] == 'm' && !read_mtl)//mtllib
			{
				sscanf_s(t += 7, "%s", mtlib);

				string realname = (string)mtlib;
				ReadMtl(realname, mtlMap);
				read_mtl = true;
			}
		}
		//cout << v.size() << "\n";
		//cout << vt.size() << "\n";
		f.clear();
		f.seekg(0, ios::beg);
		int fi = -1;
		int line_ind = 0;
		while (f.getline(line, 1024))
		{
			char* t = line;
			int prev_space = strspn(t, " \t");
			t += prev_space;
			if (strncmp(t, "f", 1) == 0)
			{
				//t += 1;
				int post_space = strspn(t + 1, " \t");

				t += post_space + 1;


				int face_type;
				vector<Face> faces;
				get_face_index(t, num_v, num_vt, faces);
				//cout << faces.size() << "\n";
				//for (int i = 0; i < faces.size(); ++i)
				//	fs.emplace_back(faces[i]);

				for (int i = 0; i < faces.size(); ++i)
				{
					++fi;
					int v0 = faces[i].v[0];
					int v1 = faces[i].v[1];
					int v2 = faces[i].v[2];
					vec3 p0(v[v0]);
					vec3 p1(v[v1]);
					vec3 p2(v[v2]);

					int vt0 = faces[i].vt[0];
					int vt1 = faces[i].vt[1];
					int vt2 = faces[i].vt[2];

					//if (!isLight)
					//{
					Triangle tri(p0, p1, p2, vt0, vt1, vt2, mtl);
					trs.emplace_back(tri);
					//}
					//else
					if (isLight)
					{
						Ls.add_light(p0, p1, p2);
						Ls.Ke.emplace_back(mats[mtl].Ke);
					}
				}
			}

			else if (strncmp(t, "usemtl", 6) == 0)//(t[0] == 'u')
			{
				sscanf_s(t += 7, "%s", name);
				string realname = (string)name;

				mtl = mtlMap.find(realname)->second;

				vec3 Ke = mats[mtl].Ke;

				if (Ke.maxc() > 0)
					isLight = true;
				else
					isLight = false;
			}
		}

		Ls.norm();
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
			if (best_split == 0 || best_split == range)//turn into leaf	//42.232s 42.393s
			{
				best_split = range / 2;			
			}
			n->axis = best_axis;

			sort(simp.begin() + start, simp.begin() + start + range, [best_axis](const SimpleTriangle& s1, const SimpleTriangle& s2)
			{
				return s1.c[best_axis] < s2.c[best_axis];
			});

			sah(n->nodes[0], start, best_split, boxes, simp);
			//sah(n->left, start + best_split, range - best_split, boxes, simp);
			sah(n->nodes[1], start + best_split, range - best_split, boxes, simp);
		}
	}

	void flat_bvh(node*& n, vector<FlatNode>& flat_nodes)
	{
		queue<node*> queue_node;

		queue_node.emplace(n);
		int current_index = 0;
		int left_index = 0;

		while (!queue_node.empty())
		{
			node* front = queue_node.front();

			queue_node.pop();

			if (front->leaf)
			{
				FlatNode flatnode;

				flatnode.box = front->box;
				flatnode.start = front->start;
				flatnode.range = front->range;
				flatnode.axis = front->axis;
				flatnode.leaf = true;

				//flat_nodes[current_index++] = flatnode;
				flat_nodes.emplace_back(flatnode);
				//++left_index;
			}
			else
			{
				FlatNode flatnode;

				flatnode.box = front->box;
				flatnode.start = ++left_index;
				flatnode.range = 0;
				flatnode.axis = front->axis;
				flatnode.leaf = false;

				//flat_nodes[current_index++] = flatnode;

				flat_nodes.emplace_back(flatnode);

				queue_node.emplace(front->nodes[0]);
				queue_node.emplace(front->nodes[1]);

				++left_index;
			}
		}

		//delete_bvh(&n);
	}

	void build_bvh(node*& root, const int& l)
	{
		leaf = l;
		int s = trs.size();

		vector<BBox> boxes;
		vector<SimpleTriangle> simp;
		boxes.resize(s);
		simp.resize(s);
		//ti.resize(s);
		for (int i = 0; i < s; ++i)
		{
			vec3 p0 = trs[i].p0;

			boxes[i].expand(p0);
			boxes[i].expand(p0 + trs[i].e1);
			boxes[i].expand(p0 + trs[i].e2);
			simp[i] = { boxes[i].c(), i };
		}
		//int depth = 0;
		sah(root, 0, s, boxes, simp);

		vector<FlatNode> flat;
		flat_bvh(root, flat);

		flat_nodes = flat;

		

		//delete_bvh(&root);
		//int left = get_tree_height(n->nodes[0]);
		//int right = get_tree_height(n->nodes[1]);

		//first_axis = left > right ? 0 : 1;

		vector<Triangle> new_trs;
		//vector<Face> new_fs;
		//new_fs.resize(fs.size());
		new_trs.resize(s);

		for (int i = 0; i < s; ++i)
		{			
			new_trs[i] = trs[simp[i].i];			
		}
		trs = new_trs;
		
		
		vector<SimpleTriangle>().swap(simp);
		
		
		vector<Triangle>().swap(new_trs);
	}

	bool __fastcall all_hit_flat_nodes(Ray& r, HitRecord& rec)
	{
		int stack[26];

		stack[0] = 0;
		float mint = 10000000.0f;
		int si = 0;
		while (si >= 0)
		{
			int top_index = stack[si--];
			auto top(flat_nodes[top_index]);
		

			float tl;

			if (top.box.hit_axis_tl(r, tl) && tl < mint)
			{
				if (!top.leaf)
				{
					int sort_axis = r.sign[top.axis];//1

					int start = top.start;
					stack[si + 1] = start;//top->nodes[sort_axis];//1
					stack[si + 2] = start + 1;//top->nodes[1 - sort_axis];//0

					if (sort_axis)
						swap(stack[si + 1], stack[si + 2]);

					si += 2;
				}
				else
				{
					//auto
					unsigned int start = top.start, end = start + top.range;//top->start + top->num;
																	  //HitRecord temp;

					for (unsigned int i = start; i < end; ++i)
					{
						trs[i].mooler_tumbor_pointer(r, rec, &mint, i);
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
		rec.t = mint;
		return true;
	}

	bool __fastcall hit_color_flat_node(Ray& r, HitRecord& rec, vec3* color)
	{
		int stack[26];

		//stack[0] = n;

		int a = r.sign[flat_nodes[0].axis];

		//int left = flat_nodes[a].start;
		stack[0] = 1 + a;//left;
		stack[1] = 2 - a;//left + 1;

		/*if (a)
		{
			stack[0] = 2;
			stack[1] = 1;
		}*/
		float mint = 10000000.0f;
		int si = 1;
		while (si >= 0)
		{
			int top_index = stack[si--];
			auto top(flat_nodes[top_index]);


			float tl;

			if (top.box.hit_axis_tl(r, tl) && tl < mint)
				//if(top->box.hit_axis_tl_mint(r, tl, mint))
			{
				//if (!top->leaf)
				//{
				int sort_axis = r.sign[top.axis];//1

												 //stack[si + 1] = top->nodes[sort_axis] ? top->nodes[sort_axis] : 0;//1  //32.2
												 //stack[si + 2] = top->nodes[1 - sort_axis] ? top->nodes[1 - sort_axis] : 0;//0

												 //32.49s
				stack[si + 1] = top.start + sort_axis;//top.nodes[sort_axis];//1
				stack[si + 2] = top.start + 1 - sort_axis;//top.nodes[1 - sort_axis];//0//32.49

				//if (sort_axis)
				//	swap(stack[si + 1], stack[si + 2]);
				//39.293s
				//stack[si + 1] = top->nodes[1 - sort_axis];
				//stack[si + 2] = top->nodes[sort_axis];

				//stack[si + 1] = !top->leaf ? top->nodes[sort_axis] : 0;   /32.6
				//stack[si + 2] = !top->leaf ? top->nodes[1 - sort_axis] : 0;

				//top->nodes[0]->box


				si += 2;

				//}
				//else
				if (top.range)
				{
					si -= 2;
					//auto start = top.start, end = start + top.range;//top->start + top->num;
																	//HitRecord temp;
					unsigned int start = top.start, end = start + top.range;
					for (unsigned int i = start; i < end; ++i)
					{
						trs[i].mooler_tumbor_pointer(r, rec, &mint, i);
					}
				}
			}
			else
				continue;

		}
		if (mint != 10000000.0f)
		{
			rec.t = mint;
			uint32_t ind = rec.ti, mtl_ind = trs[ind].mtl;

			if (mats[mtl_ind].use_texture)
			{

				int t0 = trs[ind].vt0,
					t1 = trs[ind].vt1,
					t2 = trs[ind].vt2;

				vec3 t(intp(vt[t0], vt[t1], vt[t2], rec.u, rec.v));

				int tex_ind = mats[mtl_ind].TexInd;

				*color = texs[tex_ind].ev(t);

				return true;

			}

			*color = (mats[mtl_ind].Kd.maxc() > 0 ? mats[mtl_ind].Kd : mats[mtl_ind].Ks);

			return true;
		}
		return false;

	}

	bool __fastcall hit_anything_range_flat_node(Ray& r, const float& d)
	{
		int stack[26];

		//int a = r.sign[n->axis];

		//stack[0] = n->nodes[a];
		//stack[1] = n->nodes[1 - a];



		stack[0] = 1;//n->nodes[0];
		stack[1] = 2;// n->nodes[1];

		int si = 1;
		//stack[si] = n;
		//int ti_ = Ls.ti[li];

		//int si = 0;
		//stack[0] = n;
		while (si >= 0)
		{
			int top_index = stack[si--];
			auto top(flat_nodes[top_index]);

			float th;

			if (top.box.hit_shadow_th(r, th) && th > 0.0f)
				//if(top->box.hit_shadow(r))
			{
				if (!top.leaf)
				{
					int start = top.start;
					//34s
					stack[si + 1] = start;
					stack[si + 2] = start + 1;

					//32s
					//stack[si + 1] = top->nodes[1];
					//stack[si + 2] = top->nodes[0];


					si += 2;
				}
				else
					//if(top->leaf)
				{
					//si -= 2;
					//auto start = top.start, end = start + top.range;//top->start + top->num;

																	//HitRecord temp;
					unsigned int start = top.start, end = start + top.range;
					for (auto i = start; i < end; ++i)
					{
						int mtl = trs[i].mtl;
						if (trs[i].hit_anything(r, d) && !mats[mtl].isLight)
							return true;
					}
				}
			}
			else
				continue;
		}
		return false;
	}

	vec3 __fastcall intp(const vec3& a, const vec3& b, const vec3& c, const float& u, const float& v)
	{
		return a * (1.0f - u - v) + b * u + c * v;
		//return c * (1.0f - u - v) + a * u + b * v;
		//return a + b * u + c * v;
	}
	
	float __fastcall mis2(const float& pdf1, const float& pdf2)
	{
		float pdf1_2 = pdf1 * pdf1;
		return pdf1_2 / (pdf1_2 + pdf2 * pdf2);
	}

	vec3 __fastcall path_tracing_flat_node(Ray& r)
	{
		vec3 L(0.0f);
		vec3 T(1.0f);
		HitRecord rec;
		vec3 color;

		float prev_pdf = 1.0f;

		Ray new_ray(r);
		bool specular = true;
		vec3 prev_normal;
		//vec3 Le(2.0f);
		for (int i = 0; i < 60; ++i)
		{
			//if (hit_color_pointer(n, new_ray, rec, &color, alpha))
			//if (hit_color(n, new_ray, rec, &color))
			if(hit_color_flat_node(new_ray, rec, &color))
			{
				int ti = rec.ti;

				vec3 rec_normal(trs[ti].n);

				float cos_incident = new_ray.d.dot(rec_normal);

				bool reverse_normal = false;
				if (cos_incident > 0)
				{
					reverse_normal = true;
					rec_normal = -rec_normal;
				}
				prev_normal = rec_normal;
				//vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
				int mtl = trs[rec.ti].mtl;

				if (mats[mtl].isLight)
				{
					if (specular)
					{
						//Le 111
						L += T * mats[mtl].Ke;
						//L += T * Le;
						return L;
					}
					else
					{
						int li = light_map[rec.ti];

						vec3 light_direction(new_ray.d * rec.t);

						//float length2 = light_direction.length2();

						//light_direction.normalize();


						float cos_light = -light_direction.dot(rec_normal);

						float cos_mtl = light_direction.dot(prev_normal);

						//float pdf_light = 1.0f / (Ls.area[li] * abs(cos_light * cos_mtl)) * Ls.pdf[li];
						float pdf_light = Ls.pdf[li] / (Ls.area[li] * abs(cos_light * cos_mtl));


						//float pdf_mtl = cos_mtl * ipi;

						float mis_weight = mis2(prev_pdf, pdf_light);

						//L += T * mats[mtl].Ke * prev_pdf / (prev_pdf + pdf_light);

						L += T * mats[mtl].Ke * mis_weight;

						//Le 111
						//L += T * mats[mtl].Ke * mis_weight / prev_pdf;

						//L += T * mats[mtl].Ke / prev_pdf;

						//L += T * mats[mtl].Ke / (prev_pdf + pdf_light);

						return L;

					}
				}

				vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec_normal * 0.0002f);

				int type = mats[mtl].type;

				float cont_prob = min(T.maxc(), 1.0f);

				switch (type)
				{
					//if (type == Diffuse
				case Diffuse:
				{

					//T *= color * ipi;

					//int num_light = 32;
					//vector<float> pdf_light_list;
					//vector<vec3> position;

					//pdf_light_list.resize(32);
					//position.resize(32);

					int li = Ls.sample();


					vec3 light_point(Ls.sample_light(li));

					vec3 light_direction(light_point - hit_point);


					float length = light_direction.length();



					//float ilength = 1.0f / length;

					//light_direction *= ilength;

					Ray light_ray(hit_point, light_direction);

					//float cos_mtl = abs(light_direction.dot(rec.n));
					//float cos_light = abs(-light_direction.dot(Ls.normal[li]));



					//if (!hit_anything_range(n, light_ray, length * 0.999f))
					if(!hit_anything_range_flat_node(light_ray, length * 0.999f))
					{
						float cos_mtl = (light_direction.dot(rec_normal));
						float cos_light = (-light_direction.dot(Ls.normal[li]));

						//cos_light *= ilength;

						//float pdf_light = (length * length) / (Ls.area[li] * abs(cos_light * cos_mtl)) * Ls.pdf[li];
						//vec3 bsdf_eval(color * ipi);
						//float bsdf_pdf = abs(cos_mtl) * ipi;

						float pdf_light = (length * length) / (Ls.area[li] * abs(cos_light * cos_mtl)) * Ls.pdf[li];

						//float pdf_light = Ls.pdf[li] / (Ls.area[li] * abs(cos_light * cos_mtl));

						//Original Diffuse Sampling
						vec3 bsdf_eval(color * ipi);
						float bsdf_pdf = abs(cos_mtl * ipi);

						//Simpler Sampling
						//vec3 bsdf_eval(color * ipi);
						//float bsdf_pdf = ipi;

						//vec3 bsdf_eval = color;
						//float bsdf_pdf = cos_mtl * ipi;

						float mis = mis2(pdf_light, bsdf_pdf * cont_prob);
						//float mis = mis2(pdf_light, bsdf_pdf);

						L += T * Ls.Ke[li] * bsdf_eval * mis / pdf_light;

						//return L;
						//return true;
					}


					onb local(rec_normal);


					BSDF_Sample sample = SampleDiffuse(rec_normal);

					vec3 new_dir(sample.direction.x * local.u + sample.direction.y * local.v + sample.direction.z * local.w);

					new_ray = Ray(hit_point, new_dir);

					prev_pdf = sample.pdf;
					specular = sample.is_specular;

					T *= color;
					break;
				}
				//else if (type == Mirror)
				case Mirror:
				{
					BSDF_Sample sample = SampleMirror(new_ray, rec_normal);

					new_ray = Ray(hit_point, sample.direction);
					prev_pdf = 1.0f;
					specular = true;
					break;
				}
				//else if (type == Glass)
				case Glass:
				{
					float ior = 1.6f;
					//float sign = cos_incident >= 0 ? 1.0f : -1.0f;
					bool is_Reflect;

					BSDF_Sample sample = SampleGlass(new_ray, rec_normal, ior, is_Reflect);


					vec3 intersect_point = is_Reflect ? new_ray.o + new_ray.d * rec.t + rec_normal * 0.0002f
						: new_ray.o + new_ray.d * rec.t - rec_normal * 0.0002f;

					new_ray = Ray(intersect_point, sample.direction);
					prev_pdf = 1.0f;
					specular = true;
					break;
				}
				case Conductor:
				{
					BSDF_Sample sample = SampleConductor(new_ray.d, rec_normal, mats[mtl].complex_ior);

					new_ray = Ray(hit_point, sample.direction);
					prev_pdf = 1.0f;
					specular = true;
					T *= color * sample.color;
				}
				default:
					break;
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
				return L + T;
				//return vec3(0.0f, 1.0f, 0.0f);
				//return L;// +T;
				//break;
			}
		}
		return L;
	}

};


#endif // !_VEC3_H_

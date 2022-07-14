#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define i255 1 / 255

int __fastcall clamp_int_texture(const int& i, const int& a, const int& b)
{
	//bad, cause error!!!
	//return max(a, min(i, b));



	return min(max(a, i), b);
	/*if (i <= a)
	return a;
	else if (i >= b)
	return b;
	return i;*/
}

struct Texture
{
	//string Texname = "";
	int w;
	int h;
	//int w2;
	//int h2;
	//vector<float> cs;
	//unsigned char* c;
	vector<float> c;
	
	//unsigned char* a;
	//vector<float> cs;
	//unsigned char* alpha;
	vec3 __fastcall ev(const vec3& t) const
	{
		//float u = t.x - floorf(t.x),//-floorf de loai integer : 1.099 - 1 = 0.99, !=[0, 1] ->[0, 1]
		//      v = t.y - floorf(t.y);
		

		//float u = t.x -//xs_ToInt(t.x),//xs_RoundToInt(t.x),//xs_CRoundToInt(t.x),
		//	v = t.y - //xs_ToInt(t.y);//xs_RoundToInt(t.y);//xs_CRoundToInt(t.y);

		//int x = clamp_int_texture(int(u * w), 0, w - 1),
		//    y = clamp_int_texture(int(v * h), 0, h - 1);

		int x = (t.x - floorf(t.x)) * (w - 1);
		int y = (t.y - floorf(t.y)) * (h - 1);
		//int i = x + ( h - 1 - y) * w;
		
		int i = 3 * (x + y * w);

		return vec3(c[i], c[i + 1], c[i + 2]);

		//return vec3(c[i], c[i + 1], c[i + 2]);
		//return vec3(cs[i], cs[i + 1], cs[i + 2]);
		
		
		//return vec3(powf(c[i] / 255.0f, 2.2f), powf(c[i + 1] / 255.0f, 2.2f), powf(c[i + 2] / 255.0f, 2.2f)); //*/
		
																											  
		//return vec3(powf(c[i] * 0.003921f, 2.2f), powf(c[i + 1] * 0.003921f, 2.2f), powf(c[i + 2] * 0.003921f, 2.2f));
		//return vec3(c[i] / 255.0f, c[i + 1] / 255.0f, c[i + 2] / 255.0f);
	}
	/*float __fastcall ev_alpha(const vec3& t)
	{
		float u = t.x - floorf(t.x), 
			  v = t.y - floorf(t.y);
		int x = clamp_int_texture(int(u * w2), 0, w2 - 1),
			y = clamp_int_texture(int(v * h2), 0, h2 - 1);
		//int i = x + ( h - 1 - y) * w;

		int i = 3 * (x + y * w);
		
		return a[i];//vec3(a[i]);
	}*/

	/*vec4 __fastcall ev_alpha(const vec3& t)
	{
		float u = t.x - floorf(t.x), 
			  v = t.y - floorf(t.y);
		

		int x = clamp_int_texture(int(u * w), 0, w - 1),
			y = clamp_int_texture(int(v * h), 0, h - 1);
		
		int i = 3 * (x + y * w);
		int j = 4 * (x + y * w);
		return vec4(powf(c[i] / 255.0f, 2.2f), powf(c[i + 1] / 255.0f, 2.2f), powf(c[i + 2] / 255.0f, 2.2f), alpha[j]); 
	}//*/

	/*void loadTexture(const string& p)
	{
		int n;
		c = stbi_load(p.c_str(), &w, &h, &n, 3);
	}*/
};

#endif // !_TEXTURE_H_

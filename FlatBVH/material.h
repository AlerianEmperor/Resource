#ifndef _MATERIAL_H_
#define _MATERIAL_H_
#include "Rnd.h"
#include "Texture.h"
#include "Complex_Ior.h"

enum MaterialType
{
	Diffuse,
	Mirror,
	Glass,
	Conductor
};

struct Material
{
	vec3 Kd;
	vec3 Ks;	
	vec3 Ke;
	
	
	ComplexIor complex_ior;

	//string Texname = "";
	int TexInd = -1;
	int type = Diffuse;
	bool use_texture = false;
	bool isLight = false;
	//float wd = 1.0f;
	//float ws = 1.0f;
	//float iwd = 1.0f;
	//float iws = 1.0f;
	bool use_alpha = false;
	
	float Tr = -1;
	float Ni = -1;
	//float iNi = -1;
	float Ns = -1;
	int d = -1;
	int illum;

	float R0;
	string Matname = "";

	vec3 Ka;
	vec3 Tf;

	//int mapKd;

	//Texture* mapKd = NULL;
	//int mapKa;
	//int mapKd;
	//int mapBump;
	//float ax = -1;
	//float ay = -1;
	Material() {}

};

struct BSDF_Sample
{
	vec3 direction;
	float pdf;

	vec3 color;
	bool is_specular;
	bool NULL_Sample = false;

	BSDF_Sample() {}
	BSDF_Sample(vec3 direction_, float pdf_, vec3 color_, bool is_specular_) : direction(direction_), pdf(pdf_), color(color_), is_specular(is_specular_) {}
	BSDF_Sample(vec3 direction_, float pdf_, bool is_specular_) : direction(direction_), pdf(pdf_), is_specular(is_specular_) {}
};

//---------DIFFUSE-------------

vec3 __fastcall cosine_hemisphere()
{
	float u = randf(), v = randf();
	float r = sqrt14(u);
	float theta = 2.0f * pi * v;

	vec3 direction(r * cosf(theta), r * sinf(theta), sqrt14(1.0f - u));
	return direction.norm();
}

BSDF_Sample __fastcall SampleDiffuse(vec3& n)
{	
	vec3 direction(cosine_hemisphere());

	float pdf = direction.dot(n) * ipi;

	//sample.direction = direction;
	//sample.pdf = pdf;
	//sample.is_specular = false;

	return BSDF_Sample(direction, pdf, false);
}

//---------MIRROR--------------
BSDF_Sample __fastcall SampleMirror(Ray& r, vec3& n)
{
	//vec3 direction(r.d - 2.0f * r.d.dot(n) * n);
	//sample.direction = direction;
	

	vec3 direction((r.d - 2.0f * r.d.dot(n) * n).norm());
	//sample.pdf = 1.0f;
	//sample.is_specular = true;

	return BSDF_Sample(direction, 1.0f, true);
}

//---------GLASS---------------

vec3 __fastcall Reflect(vec3& d, vec3& n)
{
	return (d - 2.0f * d.dot(n) * n).norm();
}

vec3 __fastcall Refract(vec3& d, vec3& n, float& eta, float& cos_t)
{
	float cos_i = d.dot(n);

	return (eta * d + (eta * cos_i - cos_t) * n).norm();

}

void __fastcall Dielectric_Reflectance(const float& eta, const float& cos_i, float& reflectance, float& cos_t)
{
	float sin_t2 = eta * eta * (1.0f - cos_i * cos_i);

	if (sin_t2 > 1.0f)
	{
		reflectance = 1.0f;
		cos_t = 0.0f;
		return;
	}

	cos_t = sqrt14(1.0f - sin_t2);

	float r_s = (eta * cos_i - cos_t) / (eta * cos_i + cos_t);
	float r_p = (eta * cos_t - cos_i) / (eta * cos_t + cos_i);

	reflectance = (r_s * r_s + r_p * r_p) * 0.5f;
}

BSDF_Sample __fastcall SampleGlass(Ray& r, vec3& n, float& ior, bool& is_Reflect)
{
	float cos_i = r.d.dot(n);

	float eta = cos_i >= 0 ? ior : 1.0f / ior;

	cos_i = abs(cos_i);

	float reflectance, cos_t;

	Dielectric_Reflectance(eta, cos_i, reflectance, cos_t);

	if (randf() < reflectance)
	{
		//sample.direction = reflect(r.d, n);
		//sample.pdf = reflectance;
		//sample.is_specular = true;
		//is_Reflect = true;

		vec3 direction(Reflect(r.d, n));
		is_Reflect = true;

		return BSDF_Sample(direction, reflectance, true);
	}
	else
	{
		//sample.direction = refract(r.d, -n, eta, cos_t);
		//sample.pdf = 1.0f - reflectance;
		//sample.is_specular = true;
		//is_Reflect = false;

		vec3 direction(Refract(r.d, -n, eta, cos_t));
		is_Reflect = false;
		return BSDF_Sample(direction, 1.0f - reflectance, true);
	}
}

//Rough Dielectric
void __fastcall RoughDielectric()
{

}

//Plastic
/*
BSDF_Sample __fastcall SamplePlastic(vec3& direction_in, vec3& n, const float& ior)
{
	float cos_i = -direction_in.dot(n);

	if (cos_i <= 0)
	{
		BSDF_Sample sample;
		sample.NULL_Sample = true;
		return sample;
	}

	float eta = 1.0f / ior;

	float reflectance, cos_t;
	Dielectric_Reflectance(eta, cos_i, reflectance, cos_t);

	float specular_prob = reflectance;

	if (randf() < specular_prob)
	{
		vec3 direction(Reflect(direction_in, n));
		
		return BSDF_Sample(direction, specular_prob, vec3(1.0f), true);
	}
	else
	{
		vec3 direction(cosine_hemisphere());
		return BSDF_Sample(direction, )
	}
}
*/
//Conductor

float __fastcall conductor_reflectance(float& eta, float& k, float& cos_i)
{
	float cos_i2 = cos_i * cos_i;
	float sin_i2 = 1.0f - cos_i2;
	float sin_i4 = sin_i2 * sin_i2;

	float x = eta * eta - k * k - sin_i2;
	float a2_b2 = sqrt14(x * x + 4.0f * eta * eta * k * k);

	float a = sqrt14((a2_b2 + x) * 0.5f);

	float r_s = ((a2_b2 + cos_i2) - 2.0f * a * cos_i) / ((a2_b2 + cos_i2) + (2.0f * a * cos_i));
	float r_p = ((cos_i2 * a2_b2 + sin_i4) - (2.0f * a * cos_i * sin_i2)) / ((cos_i2 * a2_b2 + sin_i4) + (2.0f * a * cos_i * sin_i2));

	return (r_s * r_s + r_p * r_p) * 0.5f;
}

vec3 __fastcall conductor_reflectance_rgb(ComplexIor& ior, float& cos_i)
{
	float x = conductor_reflectance(ior.eta.x, ior.k.x, cos_i);
	float y = conductor_reflectance(ior.eta.y, ior.k.y, cos_i);
	float z = conductor_reflectance(ior.eta.z, ior.k.z, cos_i);

	return vec3(x, y, z);
}



BSDF_Sample __fastcall SampleConductor(vec3& direction_in, vec3& n, ComplexIor& ior)
{
	float cos_i = -direction_in.dot(n);

	if (cos_i <= 0.0f)
	{
		BSDF_Sample sample;
		sample.NULL_Sample = true;
		return sample;
	}

	vec3 direction(Reflect(direction_in, n).norm());
	
	vec3 color(conductor_reflectance_rgb(ior, cos_i));

	return BSDF_Sample(direction, 1.0f, color, true);
}







//Smooth Coat
/*
BSDF_Sample __fastcall SampleSmoothCoat(vec3& direction_in, vec3& n, vec3& scale_sigma, float& ior)
{
	float eta = 1.0f / ior;

	float cos_i = -direction_in.dot(n);

	//fi = reflectance probability
	//cos_ti = cos refraction
	float fi, cos_ti;
	Dielectric_Reflectance(eta, cos_i, fi, cos_ti);

	float average_transmitance = exp(-2.0f * scale_sigma.average());

	float sub_weight = average_transmitance * (1.0f - fi);

	float specular_weight = fi;

	float specular_prob = specular_weight / (specular_weight + sub_weight);

	if (randf() < specular_prob)
	{
		//sample.direction = reflect(direction_in, n);
		//sample.pdf = specular_prob;
		//sample.color = vec3(fi + sub_weight);
		//sample.is_specular = true;

		vec3 direction(reflect(direction_in, n));

		return BSDF_Sample(direction, vec3(fi + sub_weight), specular_prob, true);
	}

	vec3 direction_in_sub(direction_in.x * eta, -cos_ti, direction_in.z * eta);

	float sub_sample
}
*/
#endif // !_MATERIAL_H_
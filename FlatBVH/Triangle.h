#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include "BBox.h"
#define eps 0.000000001f
struct HitRecord
{
	//HitRecord() {}
	vec3 n;
	float t;
	float u;
	float v;
	int ti;
};

struct Triangle
{
	
	Triangle() {}

	Triangle(vec3 p0, vec3 p1, vec3 p2, int vt0_, int vt1_, int vt2_, int m) : p0(p0), e1(p1 - p0), e2(p2 - p0), vt0(vt0_), vt1(vt1_), vt2(vt2_), n(e1.cross(e2).norm()), mtl(m) {}
	//Triangle(vec3 p0, vec3 p1, vec3 p2, int vt0_, int vt1_, int vt2_, int m) : p0(p0), e1(p1 - p0), e2(p2 - p0), n(e1.cross(e2).norm()), mtl(m) { vt[0] = vt0_, vt[1] = vt1_, vt[2] = vt2_; }

	bool __fastcall mooler_tumbor_pointer(const Ray& r, HitRecord& rec, float* mint, unsigned int& ti) const
	{
		/*
		vec3 pvec(r.d.cross(e2));
		float det = pvec.dot(e1);
		if (det > -0.000002f && det < 0.000002f)
		//if(abs(det) < 0.000002f)
			return false;
		float inv_det = 1.0f / det;

		vec3 tvec(r.o - p0);
		vec3 qvec(tvec.cross(e1));
	
		float u = tvec.dot(pvec) * inv_det;
			
		if (u <= eps || u >= 1.0f)
			return false;

		float v = qvec.dot(r.d) * inv_det;
		if (v <= eps || u + v >= 1.0f)
			return false;
		
		//float u = tvec.dot(pvec) * inv_det;
		//float v = qvec.dot(r.d) * inv_det;
		//if (u <= eps || u >= 1.0f || v <= eps || u + v >= 1.0f)
		//	return false;
		float t = e2.dot(qvec) * inv_det;
		if (t >= *mint || t < eps)
			return false;

		*mint = t;
		//rec.t = t;
		rec.n = n;
		rec.u = u;
		rec.v = v;
		rec.ti = ti;
		return true;
		*/

		vec3 pvec(r.d.cross(e2));

		//float pvec_[3] = { pvec.x, pvec.y, pvec.z };
		//float e1_[3] = { e1.x, e1.y, e1.z };

		//float det = avx_scalarproduct(pvec_, e1_, 3);
		float det = pvec.dot(e1);
		//if (det > -0.000002f && det < 0.000002f)
		if(abs(det) < 0.00000002f)
			return false;
		float inv_det = 1.0f / det;

		vec3 tvec(r.o - p0);
		vec3 qvec(tvec.cross(e1));

		float u = tvec.dot(pvec) * inv_det;

		if (u <= eps || u >= 1.0f)
			return false;

		float v = qvec.dot(r.d) * inv_det;
		if (v <= eps || u + v >= 1.0f)
			return false;

		//float u = tvec.dot(pvec) * inv_det;
		//float v = qvec.dot(r.d) * inv_det;
		//if (u <= eps || u >= 1.0f || v <= eps || u + v >= 1.0f)
		//	return false;
		float t = e2.dot(qvec) * inv_det;
		if (t >= *mint || t <= eps)
			return false;

		*mint = t;
		rec.t = t;
		rec.n = n;
		rec.u = u;
		rec.v = v;
		rec.ti = ti;
		return true;
	}

	bool __fastcall hit_anything(const Ray& r, const float& d) const
	{
		vec3 pvec(r.d.cross(e2));
		float det = pvec.dot(e1);
		//if (det > -0.000002f && det < 0.000002f)
		if(abs(det) < 0.00000002f)
			return false;
		float inv_det = 1.0f / det;

		vec3 tvec(r.o - p0);
		vec3 qvec(tvec.cross(e1));

		
		//if (u < 0.0f || u > 1.0f)

		float u = tvec.dot(pvec) * inv_det;	
		if (u < 0.0f || u >= 1.0f)
			return false;


		float v = qvec.dot(r.d) * inv_det;
		if (v <= 0.0f || u + v >= 1.0f)
			return false;

		//float t = e2.dot(qvec) * inv_det;
		//return t <= d && t >= 0.0001f;

		//float t = e2.dot(qvec) * inv_det;
		
		//if (t > d || t < 0.00001f)
		//	return false;
		//return true;

		//return t <= d && t >= eps;//>=0.00002f
		
		float t = e2.dot(qvec) * inv_det;

		//if (t > d || t < 0.00002f)
		//	return false;

		return t <= d && t >= 0.0000002f;
		//return true;
	}

	

	vec3 e1;
	vec3 e2;
	vec3 p0;
	//int vt[3];
	
	int mtl;

	int vt0;
	int vt1;
	int vt2;
	

	vec3 n;

};


/*
struct Triangle
{
	Triangle() {}
	Triangle(vec3 p0_, vec3 p1, vec3 p2, int f, int m) : p0(p0_), e1(p0 - p1), e2(p0 - p2), n(e1.cross(e2).norm()), ti(f), mtl(m) {}

	bool __fastcall mooler_tumbor_pointer(const Ray& r, HitRecord& rec, float* mint) const
	{
		//vec3 ABC = e1
		//vec3 DEF = e2
		//vec3 GHI = r.d

		vec3 JKL(p0 - r.o);

		vec3 DEF_GHI(e2.cross(r.d));

		//float denom = e1.dot(DEF_cross_GHI);

		float inv_denom = 1.0f / (e1.dot(DEF_GHI));

		float beta = JKL.dot(DEF_GHI) * inv_denom;

		if (beta <= 0.0001f || beta >= 1.0f)
			return false;

		vec3 ABC_JKL(e1.cross(JKL));

		float gamma = r.d.dot(ABC_JKL) * inv_denom;

		if (gamma <= 0.0001f || beta + gamma >= 1.0f)
			return false;

		float t = -e2.dot(ABC_JKL) * inv_denom;

		if (t >= 0.0001f && t < *mint)
		{
			*mint = t;
			rec.t = t;
			rec.u = beta;
			rec.v = gamma;
			rec.n = n;
			rec.ti = ti;
			return true;
		}

		return false;
		//beta = u
		//gamma = v

	}

	bool __fastcall hit_anything(const Ray& r, const float& d) const
	{
		//vec3 ABC = e1
		//vec3 DEF = e2
		//vec3 GHI = r.d

		vec3 JKL(p0 - r.o);

		vec3 DEF_GHI(e2.cross(r.d));

		float denom = e1.dot(DEF_GHI);

		if (denom > -0.000002f && denom < 0.000002f)
			return false;

		float inv_denom = 1.0f / (denom);

		float beta = JKL.dot(DEF_GHI) * inv_denom;

		if (beta <= 0.0001f || beta >= 1.0f)
			return false;

		vec3 ABC_JKL(e1.cross(JKL));

		float gamma = r.d.dot(ABC_JKL) * inv_denom;

		if (gamma <= 0.0001f || beta + gamma >= 1.0f)
			return false;

		float t = -e2.dot(ABC_JKL) * inv_denom;

		return t >= 0.0001f && t <= d;

	}

	vec3 e1;
	vec3 e2;
	vec3 p0;
	vec3 n;
	int ti;
	int mtl;
};
*/

#endif // !_TRIANGLE_H_


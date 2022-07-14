#ifndef _VEC3_H_
#define _VEC3_H_
#include <iostream>

//_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

//_mm_setcsr(_mm_getcsr() | 0x8040);
//#include "Distribution.h"
/*
#define min(x, y) x < y ? x : y
#define max(x, y) x > y ? x : y
#define square(x) x * x
const float eps = 0.0001f;
const float inf = 1e20f;
#define pi 3.1415f
#define ipi 0.3183f
#define i2pi 0.1591f
#define i180 0.0055f
#define abs(x) ((x)<0 ? -(x) : (x))
*/
#define min(x, y) x < y ? x : y
#define max(x, y) x > y ? x : y
//#define max3(x, y, z) max(max(x, y), z)
//#define min3(x, y, z) min(min(x, y), z)
#define abs(x) ((x)<0 ? -(x) : (x))
//#define square(x) x * x
const float eps = 0.0000001f;
//const float inf = 1000000000.0f;
//#define delta 1e-4
#define inf 1e20
const float pi = 3.1415926f;
const float ipi = 0.3183098f;
#define i2pi 0.1591f
#define i180 0.0055f


double inline __declspec (naked) __fastcall sqrt14(double n)
{
	_asm fld qword ptr[esp + 4]
		_asm fsqrt
	_asm ret 8
}

//max speed 14.121 //14.619 in bad condition, improve min max function// 14.569
//new max speed 13.69 robust bvh traversal
//13.69 = 2% leaf sorting + 1.5 % robust bvh traversal 
/*
struct vec3
{
	 vec3() : x(0), y(0), z(0) {}
	 vec3(float v) : x(v), y(v), z(v) {}
	 vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
	
	
	float x, y, z;
	
	
	float operator[](const int& i) const { return (&x)[i]; }


	friend vec3& __fastcall operator+(const vec3& a, const vec3& b)
	{
		return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
	}
    friend vec3& __fastcall operator-(const vec3& a, const vec3& b)
	{
		
		return vec3(a.x - b.x, a.y - b.y, a.z - b.z);
	}
	friend vec3& __fastcall operator*(const vec3& a, const vec3& b)
	{
		
		return vec3(a.x * b.x, a.y * b.y, a.z * b.z); 
	}
	 friend vec3& operator*=(const vec3& a, const float& v)
	{
		//float ax = a.x, ay = a.y, az = a.z;
		return vec3(a.x * v, a.y * v, a.z * v);
	}

	friend vec3& __fastcall operator*(const vec3& a, const float& v) { return vec3(a.x * v, a.y * v, a.z * v); }
	friend vec3& __fastcall operator*(const float& v, const vec3& a) { return vec3(a.x * v, a.y * v, a.z * v); }
	void operator+=(const vec3& v) { x += v.x; y += v.y; z += v.z; }
	void operator*=(const float& value) { x *= value; y *= value; z *= value; }

	float maxc() const { float d = max(x, y); return max(d, z); }//{ return max(max(x, y), z); }//
	float minc() const { float d = min(x, y); return min(d, z); }//{ return min(min(x, y), z); }//

	friend vec3& operator/(const vec3& a, const vec3& b) { return vec3(a.x / b.x, a.y / b.y, a.z / b.z); }
	friend vec3& operator/=(const vec3& a, const float& v) { return vec3(a.x / v, a.y / v, a.z / v); }
	friend vec3& operator-(const vec3& a) { return vec3(-a.x, -a.y, -a.z); }

	friend vec3& operator/(const vec3& a, const float& v) { return vec3(a.x / v, a.y / v, a.z / v); }
	friend vec3& operator/(const vec3& a, const int& v) { return vec3(a.x / v, a.y / v, a.z / v); }
	

	vec3& __fastcall norm() const
	{
		float l = 1 / sqrt14(x*x + y*y + z*z); return *this * l;
	}
	 float _fastcall dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }

	 float _fastcall length() const { return sqrt14(x * x + y * y + z * z); }
	 float _fastcall length2() const { return x * x + y * y + z * z; }
	
};
*/


//float

struct vec3
{
	vec3() : x(0), y(0), z(0) {}//{ v[0] = v[1] = v[2] = 0; }
	vec3(float v_) : x(v_), y(v_), z(v_) {}//{ v[0] = v[1] = v[2] = v_; }
	vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}//{ v[0] = x_; v[1] = y_; v[2] = z_; }

	//float a : 16;
	union
	{
		struct
		{
			float x, y, z;
		};
		struct
		{
			float v[3];
		};
	};

		vec3 __fastcall norm() const
		{
			float l = 1 / sqrt14(x*x + y*y + z*z); return *this * l;
		}

		void __fastcall normalize() const
		{
			float l = 1 / sqrt14(x*x + y*y + z*z); *this *= l;
		}

		float _fastcall dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }

		
		/*float __fastcall dot(const vec3& vec) const 
		{
			float c[3];

			for (int i = 0; i < 3; i +=3)
			{
				c[i] = v[i] * vec[i];
				c[i + 1] = v[i + 1] * vec[i + 1];
				c[i + 2] = v[i + 2] * vec[i + 2];
			}
			return c[0] + c[1] + c[2];
			//return x * v.x + y * v.y + z * v.z;
		}*/
		vec3 __fastcall cross(const vec3& b) const
		{
			//return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
		
			float bx = b.x, by = b.y, bz = b.z;
			return{ y * bz - z * by, z * bx - x * bz, x * by - y * bx };

			//return *this - this->dot(b) * *this;

		}
		float _fastcall length() const { return sqrt14(x * x + y * y + z * z); }
		float _fastcall length2() const { return x * x + y * y + z * z; }

	float operator[](const int& i) const { return (&x)[i]; }


	friend vec3 __fastcall operator+(const vec3& a, const vec3& b)
	{
		return{ a.x + b.x, a.y + b.y, a.z + b.z };
	}
	friend vec3 __fastcall operator-(const vec3& a, const vec3& b)
	{

		return{ a.x - b.x, a.y - b.y, a.z - b.z };
	}
	friend vec3 __fastcall operator*(const vec3& a, const vec3& b)
	{

		return{ a.x * b.x, a.y * b.y, a.z * b.z };
	}
	friend vec3 operator*=(const vec3& a, const float& v)
	{
		//float ax = a.x, ay = a.y, az = a.z;
		return{ a.x * v, a.y * v, a.z * v };
	}

	friend vec3 __fastcall operator*(const vec3& a, const float& v) { return{ a.x * v, a.y * v, a.z * v }; }
	friend vec3 __fastcall operator*(const float& v, const vec3& a) { return{ a.x * v, a.y * v, a.z * v }; }
	vec3 operator+=(const vec3& v) { x += v.x; y += v.y; z += v.z;  return *this; }
	void operator*=(const float& value) { x *= value; y *= value; z *= value; }
	void operator*=(const vec3& value) { x *= value.x; y *= value.y; z *= value.z; }

	float maxc() const { float d = max(x, y); return max(d, z); }//{ return max(max(x, y), z); }//
	float minc() const { float d = min(x, y); return min(d, z); }//{ return min(min(x, y), z); }//

	float average()
	{
		return (x + y + z) / 3.0f;
	}

	friend vec3 operator/(const vec3& a, const vec3& b) { return{ a.x / b.x, a.y / b.y, a.z / b.z }; }
	friend vec3 operator/=(const vec3& a, const float& v) { return{ a.x / v, a.y / v, a.z / v }; }
	friend vec3 operator-(const vec3& a) { return{ -a.x, -a.y, -a.z }; }

	friend vec3 operator/(const vec3& a, const float& v) { return{ a.x / v, a.y / v, a.z / v }; }
	friend vec3 operator/(const vec3& a, const int& v) { return{ a.x / v, a.y / v, a.z / v }; }


	
	
	bool all_zero()
	{
		return x == y == z == 0.0f;
	}
};

/*
struct vec4
{
	float x;
	float y;
	float z;
	float w;
	vec4() : x(0), y(0), z(0), w(0) {}
	vec4(float x_, float y_, float z_, float w_) : x(x_), y(y_), z(z_), w(w_) {}
	vec4(float v) : x(v), y(v), z(v), w(v) {}

};
*/

float avx_scalarproduct(float * array1, float * array2, size_t length) {
	__m256 vsum = _mm256_setzero_ps();
	size_t i = 0;
	if (length >= 8 * 8) {
		__m256 vsum1 = _mm256_setzero_ps();
		__m256 vsum2 = _mm256_setzero_ps();
		__m256 vsum3 = _mm256_setzero_ps();
		__m256 vsum4 = _mm256_setzero_ps();
		__m256 vsum5 = _mm256_setzero_ps();
		__m256 vsum6 = _mm256_setzero_ps();
		__m256 vsum7 = _mm256_setzero_ps();
		__m256 vsum8 = _mm256_setzero_ps();
		for (; i + 8 * 8 - 1 < length; i += 8 * 8) { // could unroll further
			vsum1 = _mm256_fmadd_ps(_mm256_loadu_ps(array1 + i), _mm256_loadu_ps(array2 + i), vsum1);
			vsum2 = _mm256_fmadd_ps(_mm256_loadu_ps(array1 + i + 8), _mm256_loadu_ps(array2 + i + 8), vsum2);
			vsum3 = _mm256_fmadd_ps(_mm256_loadu_ps(array1 + i + 16), _mm256_loadu_ps(array2 + i + 16), vsum3);
			vsum4 = _mm256_fmadd_ps(_mm256_loadu_ps(array1 + i + 24), _mm256_loadu_ps(array2 + i + 24), vsum4);
			vsum5 = _mm256_fmadd_ps(_mm256_loadu_ps(array1 + i + 32), _mm256_loadu_ps(array2 + i + 32), vsum5);
			vsum6 = _mm256_fmadd_ps(_mm256_loadu_ps(array1 + i + 40), _mm256_loadu_ps(array2 + i + 40), vsum6);
			vsum7 = _mm256_fmadd_ps(_mm256_loadu_ps(array1 + i + 48), _mm256_loadu_ps(array2 + i + 48), vsum7);
			vsum8 = _mm256_fmadd_ps(_mm256_loadu_ps(array1 + i + 56), _mm256_loadu_ps(array2 + i + 56), vsum8);
		}
		vsum1 = _mm256_add_ps(vsum1, vsum2);
		vsum3 = _mm256_add_ps(vsum3, vsum4);
		vsum5 = _mm256_add_ps(vsum5, vsum6);
		vsum7 = _mm256_add_ps(vsum7, vsum8);
		vsum1 = _mm256_add_ps(vsum1, vsum3);
		vsum5 = _mm256_add_ps(vsum5, vsum7);
		vsum = _mm256_add_ps(vsum1, vsum5);
	}
	for (; i + 7 < length; i += 8) { // could unroll further
		vsum = _mm256_fmadd_ps(_mm256_loadu_ps(array1 + i), _mm256_loadu_ps(array2 + i), vsum);
	}
	float buffer[8];
	_mm256_storeu_ps(buffer, vsum);
	float sum = buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];
	for (; i < length; ++i) {
		sum += array1[i] * array2[i];
	}
	return sum;
}


int dot(const vec3& a, const vec3& b)
{
	//float ax = a.x, ay = a.y, az = a.z, bx = b.x, by = b.y, bz = b.z;
	//return ax * bx + ay * by + az * bz;
	return a.x *b.x + a.y *b.y + a.z*b.z;
}

vec3 cross(const vec3& a, const vec3& b)
{
	//return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
	float ax = a.x, ay = a.y, az = a.z, bx = b.x, by = b.y, bz = b.z;
	return{ ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx };
}



#endif // !_MATH_H_


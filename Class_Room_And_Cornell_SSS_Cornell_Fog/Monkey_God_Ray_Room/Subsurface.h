#ifndef _SubSurface_Scattering_H_
#define _SubSurface_Scattering_H_

#include "BSDF.h"

//Apendix A
double C1(const double& n)
{
	double r;
	if (n > 1.0)
		r = -9.23372 + n * (22.2272 + n * (-20.9292 + n * (10.2291 + n * (-2.54396 + 0.254913 * n))));
	else
		r = 0.919317 + n * (-3.4793 + n * (6.75335 + n *  (-7.80989 + n *(4.98554 - 1.36881 * n))));
	return r * 0.5;
}
//Apendix A
double C2(const double& n)
{
	double r = -1641.1 + n * (1213.67 + n * (-568.556 + n * (164.798 + n * (-27.0181 + 1.91826 * n))));
	r += (((135.926 / n) - 656.175) / n + 1376.53) / n;
	return r / 3.0;
}

double min3(const double& x, const double& y, const double& z)
{
	const double r = x < y ? x : y;

	return r < z ? r : z;
}

class DirPole
{
public:
	//index of refraction : IOR
	const double eta = 1.4;

	//const value
	const vec3 sigma_s = vec3(0.68, 0.70, 0.55);
	const vec3 sigma_a = vec3(0.0024, 0.0090, 0.12);
	const vec3 g = vec3(0.0, 0.0, 0.0);

	//https://slidetodoc.com/subsurface-scattering-jaroslav-kivnek-ksvi-mff-uk-jaroslav/
	//Approx 1: Principal of Similarity
	const vec3 sigma_t = sigma_a + sigma_s;
	const vec3 sigma_sp = sigma_s * (1.0f - g);
	const vec3 sigma_tp = sigma_sp + sigma_a;
	const vec3 albedo_p = sigma_sp / sigma_tp;

	//solving equation diffusion
	const vec3 D = 1.0f / (3.0f * sigma_tp);
	const vec3 sigma_tr = sqrt(sigma_a / D);

	//http://people.compute.dtu.dk/jerf/papers/dirsss.pdf
	//reference 21
	const vec3 de = 2.131 * D / sqrt(albedo_p);

	//Appendex A:
	const double Cp = (1.0 - 2.0 * C1(eta)) / 4.0;
	const double Ce = (1.0 - 3.0 * C2(eta)) / 2.0;

	//The factor (4Cφ(1/η))−1 is the normalization factor also used by d’Eon and Irving
	//page 5 paper

	//Cφ(1/η))^−1 = 4.0 / (1.0 - 2.0 * C1(1.0 / eta))

	//(4)^-1 = 1 / 4

	//(4 * Cφ(1 / η)) ^ −1 = 1.0 / (1.0 - 2.0 * C1(1.0 / eta))


	const double Cp_norm = 1.0 / (1.0 - 2.0 * C1(1.0 / eta));

	const double A = (1.0 - Ce) / (2.0 * Cp);
	const double min_sigma_tr = min3(sigma_tr.x, sigma_tr.y, sigma_tr.z);

	//Directional Dipole
	//papper reference 20

	double Sp_d(const vec3& x, const vec3& w, const double& r, const vec3& n, const int& j)
	{
		// evaluate the profile
		const double s_tr_r = sigma_tr[j] * r;
		const double s_tr_r_one = 1.0 + s_tr_r;
		const double x_dot_w = x.dot(w);
		const double r_sqr = r * r;

		const double t0 = Cp_norm * (1.0 / (4.0 * pi * pi)) * exp(-s_tr_r) / (r * r_sqr);
		const double t1 = r_sqr / D[j] + 3.0 * s_tr_r_one * x_dot_w;
		const double t2 = 3.0 * D[j] * s_tr_r_one * w.dot(n);
		const double t3 = (s_tr_r_one + 3.0 * D[j] * (3.0 * s_tr_r_one + s_tr_r * s_tr_r) / r_sqr * x_dot_w) * x.dot(n);

		return t0 * (Cp * t1 - Ce * (t2 - t3));
	}

	double bssrdf(const vec3& xi, const vec3& ni, const vec3& wi, const vec3& xo, const vec3& no, const vec3& wo, const int j)
	{
		// distance
		const vec3 xoxi = xo - xi;
		const double r = xoxi.length();

		// modified normal
		const vec3 ni_s = (xoxi.norm()).cross((ni.cross(xoxi).norm()));

		// directions of ray sources
		const double nnt = 1.0 / eta, ddn = -wi.dot(ni);
		const vec3 wr = (wi * -nnt - ni * (ddn * nnt + sqrt(1.0 - nnt * nnt * (1.0 - ddn * ddn)))).norm();
		const vec3 wv = wr - ni_s * (2.0 * wr.dot(ni_s));

		// distance to real sources
		const double cos_beta = -sqrt((r * r - xoxi.dot(wr) * xoxi.dot(wr)) / (r * r + de[j] * de[j]));
		double dr;
		const double mu0 = -no.dot(wr);
		if (mu0 > 0.0) {
			dr = sqrt((D[j] * mu0) * ((D[j] * mu0) - de[j] * cos_beta * 2.0) + r * r);
		}
		else {
			dr = sqrt(1.0 / (3.0 * sigma_t[j] * 3.0 * sigma_t[j]) + r * r);
		}

		// distance to virtual source
		const vec3 xoxv = xo - (xi + ni_s * (2.0 * A * de[j]));
		const double dv = xoxv.length();

		// BSSRDF
		const double result = Sp_d(xoxi, wr, dr, no, j) - Sp_d(xoxv, wv, dv, no, j);

		// clamping to zero
		return (result < 0.0) ? 0.0 : result;
	}

	bool sample(const vec3& direction_in, vec3& n, vec3& original_n, Sampling_Struct& sampling, bool& isReflect, vec3& color) const
	{

	}

	float pdf(const vec3& direction_in, const vec3& direction_out, const vec3& n, vec3& original_n) const
	{

	}

	vec3 eval(const vec3& direction_in, const vec3& direction_out, const vec3& n, vec3& original_n, vec3& color) const
	{

	}

	bool isSpecular() const
	{
		return false;
	}






};



#endif // !_SubSurface_Scattering_H_


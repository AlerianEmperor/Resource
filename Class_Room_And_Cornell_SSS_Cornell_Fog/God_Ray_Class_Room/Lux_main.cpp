//#include "Scene.h"
//#include "math.h"
//#include "Scene_Prestigious_Path_Tracing.h"
//#include "Smooth_Normal_Scene.h"
//#include "Scene.h"
//#include "Scene_Fast.h"

//#include "Scene_SBVH.h"

//#include "Scene_Zero_Condition.h"

#include "Scene.h"
//#include "Scene_Volume_2.h"

//#include "Volume.h"
#include "Tone_Mapping.h"
#include "Filter.h"
#include <omp.h>
#include <time.h>

//http://www.raytracerchallenge.com/bonus/texture-mapping.html
//87002 usemtl drawing.002 syndra
//99569 usemtl drawing.004 sohee

using namespace std;
#define Width 960
#define Height 540

#define iWidth 1 / Width
#define iHeight 1 / Height
#define aspect_ratio Width / Height
#define golden_ratio 1.61803398875


int ns = 8192;
int step = 16;
int sqrt_ns = sqrt(ns);
float isqrt_ns = 1 / sqrt_ns;

const int num_shadow = 400;
float ins = 1.0f / ns;


//cornell box
//vec3 look_from(-3.0660074, 2.724191, 2.846079);
//vec3 view_direction(0.043121, -0.046571, -0.997984);


//vec3 look_from(-3.255610, 2.725992, 3.18545);
//vec3 view_direction(0.028125, -0.061548, -0.997708);

//vec3 look_from(-3.055610, 2.725992, 3.18545);
//vec3 view_direction(0, 0, -1);

//-2.78 2.74177 -2.79818
//vec3 look_from(-2.78, 2.74177, 2.18545);
//vec3 view_direction(0.0, 0.0, -1.000);

//God Ray Box
//vec3 look_from(-3.0660074, 3.924191, 2.846079);
//vec3 view_direction(0.043121, -0.066571, -0.997984);

//Huge Box

//vec3 look_from(-89.285355, 21.785948, -1.892215);
//vec3 view_direction(0.992008, -0.0991464, 0.086916);

//vec3 look_from(-2.78, 2.44177, 4.09818);
//vec3 view_direction(0, 0, -1);

//float fov = 60.0f;


//ClassRoom
vec3 look_from(1.888095, 1.169218, 3.067723);
vec3 view_direction(-0.054527, 0.000239, -0.998512);

//Volumetric Lighting Room
//vec3 look_from(-3.890314, 3.042109, 1.319433);
//vec3 view_direction(-0.700312, -0.187293, -0.688828);

float fov = 60.0f;

const vec3 w = -view_direction;
const vec3 up = vec3(0, 1, 0);
const vec3 u = up.cross(w).norm();
const vec3 v = w.cross(u);
const float tan_theta = tanf(fov * 0.5f * pi / 180.0f);

const float asspect_tant = aspect_ratio * tan_theta;



static void new_multi_jitter(vector<vec2>& sample)
{
	int sqrt_sample = sqrt_ns;
	float sub_cell_width = 1.0f / float(ns);

	for (int i = 0; i < sqrt_sample; ++i)
	{
		for (int j = 0; j < sqrt_sample; ++j)
		{
			sample[i * sqrt_sample + j].x = i * sqrt_sample * sub_cell_width + j * sub_cell_width + randf() * sub_cell_width;
			sample[i * sqrt_sample + j].y = j * sqrt_sample * sub_cell_width + i * sub_cell_width + randf() * sub_cell_width;

		}
	}

	for (int i = 0; i < sqrt_sample; ++i)
	{
		for (int j = 0; j < sqrt_sample; ++j)
		{
			int k = j + int(randf() * (sqrt_sample - j - 1));
			swap(sample[i * sqrt_sample + j].x, sample[i * sqrt_sample + k].x);
		}
	}

	for (int i = 0; i < sqrt_sample; ++i)
	{
		for (int j = 0; j < sqrt_sample; ++j)
		{
			int k = j + int(randf() * (sqrt_sample - j - 1));
			swap(sample[j * sqrt_sample + i].x, sample[k * sqrt_sample + i].x);
		}
	}
}


static vec3 denan(const vec3& v)
{
	vec3 temp(v);
	if (!(temp.x == temp.x))
		temp.x = 0.0f;
	if (!(temp.y = temp.y))
		temp.y = 0.0f;
	if (!(temp.z == temp.z))
		temp.z = 0.0f;

	return temp;
}

static float Luminance(const vec3& v)
{
	return 0.2126f * v.x + 0.7152f * v.y + 0.0722 * v.z;
}


int ind = 9;

void de_purple(string& file_name)
{
	ifstream ifs(file_name);

	if (!ifs)
		cout << "invalid";

	char line[512];

	ifs.getline(line, 512);
	ifs.getline(line, 512);
	ifs.getline(line, 512);

	vector<vec3> color;

	int count = 0;

	
	

	while (ifs.getline(line, 512))
	{
		//++count;
		float x, y, z;

		ifs >> x >> y >> z;
		//ifs >> y;
		//ifs >> z;
		//cout << x << " " << y << " " << z << "\n";
		color.push_back(vec3(x, y, z));
	}

	//getchar();
	int size = color.size();
	cout << size;

	

	//getchar();
	int dx[4] = { 1, 0, -1, 0 };
	int dy[4] = { 0, 1,  0, -1 };
	for (int i = 0; i < size; ++i)
	{
		if (color[i].y == 0 && color[i].x > 0 && color[i].z > 0)
		{
			int x_coord = i % Width;
			int y_coord = i / Width;

			int count_none_purple = 0;
			
			vec3 sum(0.0f);

			for (int j = 0; j < 4; ++j)
			{
				int new_x = x_coord + dx[j];
				int new_y = y_coord + dy[j];

				if (new_x >= 0 && new_x < Width && new_y >= 0 && new_y < Height)
				{
					int ind = new_y * Width + new_x;

					if (color[ind].y > 0.0f)
					{
						++count_none_purple;
						sum += color[ind];
					}
				}
			}

			if (count_none_purple >= 3)
				color[i] = sum / count_none_purple;
		}
	}

	ofstream ofs(file_name + "De_Purple.ppm");

	ofs << "P3\n" << Width << " " << Height << "\n255\n";

	for (int i = 0; i < size; ++i)
		ofs << color[i].x << " " << color[i].y << " " << color[i].z << "\n";

}

void main_2()
{
	//string file_name = "Volume_10240_10240_9.ppm";

	//string file_name = "Class_Room_God_Ray_Bigger_god_ray_2_All_Diff_16_9.ppm";

	string file_name = "1_Cornell_10240_10240_9.ppm";

	de_purple(file_name);
	getchar();
}

void main()
{
	sqrt_ns = sqrt(ns);
	isqrt_ns = 1 / sqrt_ns;
	ins = 1.0f / ns;

	look_from -= view_direction * 1.0f;

	clock_t t1;
	t1 = clock();

	//Scene scn("E:\\Models\\VolumeTric_Lighting_Room\\Volumetric_light_Monkey.obj", false, true, vec3(-1.0f, 1.0f, -0.2f), 0.04, true, false, "Sky_19.hdr", 26.0f);// , look_from, false);

	//Scene scn("E:\\Models\\VolumeTric_Lighting_Room\\Volumetric_light_Monkey.obj");

	//Scene scn("E:\\Models\\ClassRoom\\ClassRoom_No_Window\\classroom_No_Window.obj", false, true, vec3(1.0f, 1.0f, 0.2f), 0.01, true, false, "Sky_19.hdr", 26.0f);// , look_from, false);
																									  //1.0f, 0.6f, 0.2f, 0.01
	//Scene scn("E:\\Models\\ClassRoom\\ClassRoom_No_Window\\classroom_No_Window.obj", false, true, vec3(1.0f, 0.5f, 0.2f), 0.02, true, false, "Sky_19.hdr", 26.0f);// , look_from, false);

	//Scene scn("E:\\Models\\Cornell_Box\\Huge_Box_chimney.obj", false, true, vec3(0.0f, 1.0f, 0.0f), 0.01, true, false, "Sky_19.hdr", 26.0f);

	//Scene scn("E:\\Models\\Cornell_Box\\Gox_Box_No_Front.obj");

	//Scene scn("E:\\Models\\Cornell_Box\\God_Box_Small_Hole.obj", false, true, vec3(0.0f, -1.0f, 0.0f), 0.01, true, false, "Sky_19.hdr", 26.0f);

	//Scene scn("E:\\Models\\Stanford_Bunny\\Cornell_Bunny\\Cornell_Bunny.obj");
	
	//Scene scn("E:\\Models\\ClassRoom\\ClassRoom_No_Window\\classroom_No_Window.obj", "Sky_19.hdr");

	Scene scn("classroom_No_Window.obj", false, true, vec3(1.0f, 0.2f, 0.2f), 0.02, true, false, "", 26.0f);// , look_from, false);

	//Scene scn("E:\\Models\\ClassRoom\\ClassRoom_No_Window\\classroom_sohee_syndra.obj", false, true, vec3(1.0f, 0.2f, 0.2f), 0.02, true, false, "Sky_19.hdr", 26.0f);

	scn.use_directional_light = true;
	scn.compute_sky_color = true;
	//Scene scn("E:\\Models\\ClassRoom\\ClassRoom_No_Window\\classroom_LightSource.obj", false, true, vec3(1.0f, 0.2f, 0.2f), 0.01, true, false, "Sky_19.hdr", 26.0f);// , look_from, false);
	

	//Scene scn("God_Ray_Room_High_Monkey.obj", false, true, vec3(-1.0f, 1.0f, -0.2f), 0.04, true, false, "Sky_19.hdr", 26.0f);
	//Scene scn("God_Ray_Room_High_Monkey.obj", false, true, vec3(-1.0f, 1.0f, -0.2f), 0.04, true, false, "Sky_19.hdr", 26.0f);

	//Scene scn("E:\\Models\\Cornell_Box\\cornell-box.obj");

	//Scene scn("E:\\Models\\Cornell_Box\\Cornell_Box_Triangulate_Fix_Short_Box.obj");

	
	//Scene scn("Cornell_Box_Triangulate_Fix_Short_Box.obj");

	//Scene scn("E:\\Models\\Cornell_Box\\Cornell-Box_Final.obj");

	//Scene scn("cornell_box.obj");

	//Scene scn("E:\\Models\\Stanford_Bunny\\Cornell_Bunny\\Cornell_Bunny.obj");

	//Scene scn("E:\\Models\\ClassRoom\\classroom\\textures\\classroom.obj");

	scn.blur_directional_light = true;
	
	scn.use_directional_light = true;
	
	//scn.directional_blur_radius = 0.04f;
	//scn.infinite_light_direction = vec3(-1.0, -0.2, 0.1);
	
	//scn.sun_color = vec3(30, 20, 10);

	//scn.sun_color = vec3(15, 10, 5);

	//https://www.color-name.com/sunlight.color#:~:text=Sunlight%20has%20the%20hex%20code,a%20brightness%20value%20of%2096%25.
	//scn.sun_color = vec3(39, 37, 25);

	//scn.sun_color = vec3(19.5, 18.5, 12.5);

	//scn.sun_color = vec3(9.75, 9.25, 6.25);

	scn.use_enviroment = true;
	scn.blur_enviroment = true;

	//
	
	//Scene scn("E:\\Models\\Cornell_Box\\cornell-box.obj");

	//Scene scn("cornell-box.obj");

	t1 = clock() - t1;
	cout << "Reading time: " << t1 / 1000.0f << "s" << "\n";

	clock_t t_build;
	t_build = clock();


	cout << sizeof(node) << "\n";// 
	cout << scn.trs.size() << "\n";
	node* n;

	
	scn.build_bvh(n, 4);

	vec3 center = n->box.c();
	cout << center.x << " " << center.y << " " << center.z << "\n";

	//box
	
	//-2.78 2.74177 -2.79818


	//Class Room
	//-0.579745 1.54874 0.675832

	scn.delete_before_render();

	scn.medium.bound = n->box;

	t_build = clock() - t_build;
	cout << "Building time: " << t_build / 1000.0f << "s" << "\n";


	int size = Width * Height;

	vector<vec3> c;
	c.resize(size);

	vector<int> sample_count;
	sample_count.resize(size);

	vector<vec3> raw_color;
	raw_color.resize(size);






	vector<bool> visible(Width * Height, true);
	
	int num_sample = 6;
	int num_shadow = 12;

	//vec3 prev_normal(view_direction);
	omp_set_num_threads(512);

	bool sort_shadow = false;

	if (sort_shadow && scn.use_area_light)
	{
		clock_t t2 = clock();
		int count = 0;
		int count_light = 0;

		
		t2 = clock() - t2;
		cout << "Visibility check: " << t2 / 1000.0f << "s" << "\n";

		
	}



	clock_t t_render = clock();

	
	//cout <<"Light size: "<< scn.Ls.Ke.size();
	//getchar();
	for (int i = 0; i < Width; ++i)
	{
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", ns, 100.0f * i / (Width - 1));
		#pragma omp parallel for schedule(dynamic)
		for (int j = 0; j < Height; ++j)
		{
			//vector<vec2> sample;
			//sample.resize(ns);

			//new_multi_jitter(sample);

			
			vec3 sum(0.0f);

			
			float sum_sqr = 0.0f;
			float sum_so_far = 0.0f;

			int pixel_index = j * Width + i;

			int num_sample_use = 0;
			float convergence_rate = visible[pixel_index] ? 0.001f : 0.0001f;
			bool converge = false;
			
			for (int s = 0; s < ns; s += step)
			{
				for (int num = 0; num < step; ++num)
				{

					float p = ((float)i + randf()) * iWidth;
					float q = ((float)j + randf()) * iHeight;

					p = (2.0f * p - 1.0f) * aspect_ratio * tan_theta;
					q = (1.0f - 2.0f * q) * tan_theta;

					Ray new_ray(look_from, (u * p + v * q - w).norm());
					
					//vec3 current_color(denan(scn.Radiance_Volume(n, new_ray)));

					//vec3 current_color(denan(scn.sibenik_tracing_artificial_light_normal_Ke_SubSurface(n, new_ray)));

					

					//vec3 current_color(denan(scn.god_ray_pbrt_test_blur_2(n, new_ray)));

					//vec3 current_color(denan(scn.god_ray_pbrt(n, new_ray)));

					//God Ray
					//vec3 current_color(denan(scn.god_ray_2(n, new_ray)));
					//vec3 current_color(denan(scn.god_ray(n, new_ray)));

					//Subsurface Scattering

					//1 failed
					//vec3 current_color(denan(scn.Subsurface_Scattering_2_brute_force_original_modify(n, new_ray)));

					//2 failed
					//vec3 current_color(denan(scn.Subsurface_Scattering_2_brute_force_original(n, new_ray)));

					//3
					//vec3 current_color(denan(scn.Subsurface_Scattering_2_brute_force(n, new_ray)));

					//4
					//vec3 current_color(denan(scn.Subsurface_Scattering_2(n, new_ray)));

					//5
					//vec3 current_color(denan(scn.Subsurface_Scattering_2_fix_Reflection(n, new_ray)));

					//6
					//vec3 current_color(denan(scn.Subsurface_Scattering_2_fix(n, new_ray)));

					//scatter at NEE
					//vec3 current_color(denan(scn.sibenik_tracing_directional_lighting_only_scatter_if_see_light(n, new_ray)));
					//vec3 current_color(denan(scn.sibenik_tracing_directional_lighting_NEE_scatter_original_2(n, new_ray)));
					//vec3 current_color(denan(scn.sibenik_tracing_directional_lighting_NEE_scatter_different_sigma_sun(n, new_ray)));
					//vec3 current_color(denan(scn.sibenik_tracing_directional_lighting_NEE_scatter_original(n, new_ray)));
					//vec3 current_color(denan(scn.sibenik_tracing_directional_lighting_NEE_scatter_2(n, new_ray)));
					//vec3 current_color(denan(scn.sibenik_tracing_directional_lighting_NEE_scatter(n, new_ray)));

					//original
					//vec3 current_color(denan(scn.sibenik_tracing_directional_lighting(n, new_ray)));

					//----box----

					//vec3 current_color(denan(scn.sibenik_tracing_artificial_light_normal_Ke_SubSurface(n, new_ray)));

					//vec3 current_color(denan(scn.VolumeTric_Rendering_Homogeneous_Attenuate_NEE_Glass_Sphere_Direct_Lighting(n, new_ray)));
					//vec3 current_color(denan(scn.VolumeTric_Rendering_Homogeneous_Attenuate_NEE_Glass_Sphere_2(n, new_ray)));
					//vec3 current_color(denan(scn.VolumeTric_Rendering_Homogeneous_Attenuate_NEE_Glass_Sphere_2_NEE(n, new_ray)));
					//vec3 current_color(denan(scn.VolumeTric_Rendering_Homogeneous_Attenuate_NEE_Glass_Sphere_2_NEE_fog(n, new_ray)));

					//vec3 current_color(denan(scn.VolumeTric_Rendering_Homogeneous_Attenuate_NEE_Glass_Sphere(n, new_ray)));
					//vec3 current_color(denan(scn.VolumeTric_Rendering_Homogeneous_Attenuate_NEE(n, new_ray)));
					//vec3 current_color(denan(scn.VolumeTric_Rendering_Homogeneous_2(n, new_ray)));
					//----box----


					//vec3 current_color(denan(scn.sibenik_tracing_blur_enviroment(n, new_ray)));

					//vec3 current_color(denan(scn.sibenik_tracing_artificial_light_normal_Ke_alpha(n, new_ray)));

					//vec3 current_color(denan(scn.sibenik_tracing_artificial_light_normal_Ke(n, new_ray)));
					
					//vec3 current_color(scn.sibenik_tracing_blur_enviroment(n, new_ray));

					

					//vec3 current_color(denan(scn.VolumeTric_Rendering_Homogeneous_All_Transmitance(n, new_ray)));

					//vec3 current_color(denan(scn.VolumeTric_Rendering_Homogeneous_Directional(n, new_ray)));


					//chay God Ray
					vec3 current_color(denan(scn.god_ray_pbrt_test_blur(n, new_ray)));
					
					//chay direction
					//vec3 current_color(denan(scn.sibenik_tracing_directional_lighting_2(n, new_ray)));
					
					float lux = Luminance(current_color);

					sum_sqr += lux * lux;
					sum_so_far += lux;

					sum += current_color;
				}
				num_sample_use += step;

				float mean_sum = sum_so_far / num_sample_use;

				float variance = (sum_sqr / num_sample_use - mean_sum * mean_sum) / (num_sample_use - 1);

				vec3 value(sum / num_sample_use);


				if (variance < convergence_rate && value.minc() > 0.2f)
				{
					c[j * Width + i] = sum / num_sample_use;
					raw_color[j * Width + i] = sum;
					sample_count[j * Width + i] = num_sample_use;
					converge = true;
					break;
				}
			}

			if (!converge)
			{
				c[j * Width + i] = sum * ins;
				raw_color[j * Width + i] = sum;
				sample_count[j * Width + i] = ns;
			}
		}
	}
	
	string image_name = "Final_ClassRoom_8192_Waifu.ppm";

	//string image_name = "ClassRoom_Waifu_2_big_fog_sun_sigma_a_0.6_sun_sigma_s_0.05_god_ray_9_6_3_8192_sky_4.ppm";

	//string image_name = "god_ray_Class_Room_Final_Test_Before_Render_fog_2_final_Real_Class_Room.ppm";

	//string image_name = "a_a_a_a_a_a_Subsurface_2.ppm";

	//string image_name = "Test_Diffuse_Orenayar_roughness_0.9_with_cos_o_weight_correct_cos_a_minus_b_formular_abs_cos_theta";

	//string image_name = "Radiance_Volume_pdf_scatter_2_box";

	//string image_name = "God_Ray_Room_High_Monkey";

	//string image_name = "Cornell_SubSurface_divide_pdf_light";

	//string image_name = "Cornell_10240_sigma_a_0_sigma_s_0.12";

	//string image_name = "Cornell_Box_Final_NEE_no_sigma_a_attenuate";

	//string image_name = "Class_Room_Blur_God_Ray_Correct_3_2_1_1.0_0.5_0.2_sun_color_god_ray_seperate_value";

	//string image_name = "Volume_Room_Blur_Shadow_0.04_no_blur_at_all";

	//string image_name = "Class_Room_Blur_Ray_Blur_set_to_True_pbrt_2_Correct_Final";//_pi_4_r^2";

	//string image_name = "God_Ray_Class_Room_pdf_light_4_pi_radius_2_3000_2000_1000";

	//string image_name = "God_Ray_Room_Volumetric_Lighting_3_2_1_pdf_light_1_Change_Light_Direction_-1.0,_1.0_-0.2_keep_weight_at_scatter";

	//string image_name = "God_Ray_1_pbrt_original_Correct_No_Blur_3000_2000_1000";

	//string image_name = "Class_Room_God_Ray_PBRT_Final_No_Blur_No_Sigma_s_pdf";
	

	//string image_name = "Cornell_God_Box_No_Front_Direct_Lighting_Scatter_attenuate_mis_weight_directional";

	//string image_name = "Class_Room_Nee_Original_2_Correct_No_Light_Source_0.036_sun_0.036_no_T_scatter_Big";

	//string image_name = "Class_Room_God_Ray_Bigger_god_ray_2_All_Diff";

	//string image_name = "Class_Room_God_Ray_Scatter_If_See_Light_Original_0.06_0.002";

	//string image_name = "Calss_Room_God_Ray_2_blur_radius_0.2_sun_high_scatter";

	//string image_name = "Class_Room_God_Ray_different_sigma_for_sun_no_scatter_at_light_direction";

	//string image_name = "Cornell_Subsurface_2_Fix_Reflection_s_20.0_a_0.0005";

	//string image_name = "Cornell_SubSurface_2_original_return vec3_0_s_64.0_a_0.0005_brute_force_original_modify";

	//string image_name = "Cornell_Box_a_0.000_s_0.12_ior_1.9";

	//string image_name = "Cornell_Box_Yellow_Light_a_0.000_s_0.12_Super_Bright";

	//string image_name = "Bunny_SSS_No_Medium";

	//string image_name = "Cornell_SSS_Correct_prev_mtl_yellow";

	//string image_name = "Cornell_SSS_fix_2_reflect_100_percent_subsurface_sigma_Kd_a_0.005_s_0.02_tall_box_Fix_Reflection_F_1F";

	//string image_name = "Cornell_Normal_a_0.005_s_0.09_Glass_Sphere_Correct_Less_Closer_2";

	//string image_name = "Cornell_Glass_Fresnel_Direction";
	//string image_name = "Cornell_Bunny";

	//string image_name = "1_Cornell_Box_Correct_Short_Box_Glass_Realistic_No_Volume_ior_1.5_Final";

	//string image_name = "1_Class_Room_NEE_2_scatter_a_0.006_s_0.036";

	//string image_name = "1_Cornell_Volume_Triangulate_a_0.006_s_0.036_64";

	//string image_name = "1_Cornell_Volume_NEE_no_move_back";

	//string image_name = "1_Class_sun_realistic_1.0_0.2_0.2_blur_0.01_Medium_2_0.036_sigma_s_Metal_Chair_Cr_Big";

	//string image_name = "1_Class_Room_Correct_Volume_960_540_move_back_1_Sun_Lighting_2_sun_direction_good_Sky_19_SanMiguel_1.0_0.4_0.2";

	//string image_name = "a_Light_Normale_Ke_Fog_No_inv_p_0.06_0.095";

	//string image_name = "Volume_10240_Standard_Box_Homogeneous_2_sigma_a_0.006_sigma_t_0.036_ScaleBy";

	//string image_name = "Volume_10240_Standard_Buox_increase_sigma_s_0.096_add_T_0.05_out";

	//string image_name = "Volume_Scatter_No_Sample_Light_No_sigma_h_sigma_a_0.006_sigma_s_0.036_Box_Glass";

	//string image_name = "Volume_Small_VPT_Fix";

	//string image_name = "Directional_Lighting";
	//string image_name = "1.sigma_vec3_standard";

	//string image_name = "a_a_a_Fix_Error_Transmitance_Get_Transmitance_Back";//"a_g_0.6_negative_new_ray";

	clock_t t_write = clock();

	string num_sample_image = to_string(ns);
	string image_ind = to_string(ind);

	
	ofstream ofs(image_name + '_' + num_sample_image + '_' + image_ind + ".ppm");
	//ofstream ofs_raw(image_name + "_raw_" + num_sample_image + '_' + image_ind + ".txt");


	ofs << "P3\n" << Width << " " << Height << "\n255\n";


	for (int i = 0; i < size; ++i)
	{
		//vec3 raw = raw_color[i];
		//ofs_raw << raw.x << " " << raw.y << " " << raw.z << " " << sample_count[i] << "\n";
		
		//vec3 color = vec3(powf(c[i].x, 1.0f / 2.2f), powf(c[i].y, 1.0f / 2.2f), powf(c[i].z, 1.0f / 2.2f));//Filmic_unchart2(c[i]);
		
		vec3 color = Filmic_unchart2(c[i]);

		color *= 255.99f;

		ofs << color.x << " " << color.y << " " << color.z << "\n";
	}

	ofs.close();


	vector<vec3>().swap(c);
	vector<vec3>().swap(raw_color);
	vector<int>().swap(sample_count);
	
	delete_bvh(&n);
	scn.delete_struct();

	t_write = clock() - t_write;


	t_render = clock() - t_render;

	float render_time = ((double)t_render) / CLOCKS_PER_SEC;
	float write_time = ((double)t_write) / CLOCKS_PER_SEC;

	

	std::cout << "\nRendering time: " << render_time << "\n";

	ofstream ofs_log(image_name + "_log.txt");

	ofs_log << "Render Time: " << render_time << "\n";
	ofs_log << "Write Time : " << write_time << "\n";

	
}

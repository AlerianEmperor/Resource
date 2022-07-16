//#include "Scene.h"
//#include "math.h"
//#include "Scene_Prestigious_Path_Tracing.h"
//#include "Smooth_Normal_Scene.h"
//#include "Scene.h"
//#include "Scene_Save_Space.h"
#include "Scene_16.h"
//#include "Scene_16_No_Compress_Vertex.h"
#include "Tone_Mapping.h"
#include "Filter.h"
#include <time.h>
#include <omp.h>

//http://www.raytracerchallenge.com/bonus/texture-mapping.html

using namespace std;
#define Width 200
#define Height 200

#define iWidth 1 / Width
#define iHeight 1 / Height
#define aspect_ratio Width / Height
#define golden_ratio 1.61803398875

int ns = 64;
const int step = 16;
int sqrt_ns = sqrt14(ns);
float isqrt_ns = 1 / sqrt_ns;

const int num_shadow = 400;
#define ins  1 / ns

//original
vec3 look_from(1899.319, 1158.376, -4950);//-4958.338);
vec3 view_direction(0.013167, 0.0, -0.998763);

//vec3 view_direction(0.013167, -0.047944, -0.998763);
//vec3 look_from(1028.566, 699.557, -7470.727);
//vec3 view_direction(-0.983, -0.175, -0.038370);

//picture
//vec3 look_from(1311.808716, 1673.683594, -8711.320313);
//vec3 view_direction(-0.998753, -0.036582, -0.033959);

//rem
//vec3 look_from(1781.678833, 1496.433716, -7602.271484);
//vec3 view_direction(0.959052, -0.051572, -0.278495);

//Test Env
//vec3 look_from(1781.678833, 1496.433716, -7602.271484);
//vec3 view_direction(0.959052, -0.051572, -0.278495);

//test den
//vec3 look_from(2866.224121, 2554.399902, -9921.352539);
//vec3 view_direction(0.464310, 0.270013, -0.843510);

//Test Cua Sau 
//vec3 look_from(1935.008667, 1177.956421, -5388.114746);
//vec3 view_direction(0.994702, 0.038399, 0.095359);

//curtain view
//vec3 look_from(3762.888184, 2010.962646, -6864.411133);
//vec3 view_direction(0.321776, -0.205122, -0.924329);

float fov = 74.0f;



const vec3 w = -view_direction;
const vec3 up = vec3(0, 1, 0);
const vec3 u = up.cross(w).norm();
const vec3 v = w.cross(u);
const float tan_theta = tanf(fov * 0.5f * pi  / 180.0f);

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

void __fastcall determine_visibility(node*& n, Scene& scn, vector<bool>& visible, vector<bool>& in_light, int num_sample, int num_shadow, int& count, int& count_light)
{
	for (int j = 0; j < Height; ++j)
	{
		#pragma omp for schedule(guided)
		for (int i = 0; i < Width; ++i)
		{
			int light_block = 0;

			float u_index = i;
			float v_index = j;

			for (int s = 0; s < num_sample; ++s)
			{
				HitRecord rec;

				float p = (u_index + randf()) * iWidth;
				float q = (v_index + randf()) * iHeight;

				p = (2.0f * p - 1.0f) * aspect_ratio * tan_theta;
				q = (1.0f - 2.0f * q) * tan_theta;
				Ray r(look_from, u * p + v * q - w);

				bool b = scn.all_hit(n, r, rec);
				int mtl = scn.trs[rec.ti].mtl;

				if (b && !scn.mats[mtl].isLight)//
				{
					//cout << "a";
					int ti = rec.ti;
					
					//int vn0 = scn.trs[ti].vn[0];
					//int vn1 = scn.trs[ti].vn[1];
					//int vn2 = scn.trs[ti].vn[2];

					vec3 rec_normal = rec.n;//vn0 < 0 ? scn.trs[ti].n : scn.intp(scn.vn[vn0], scn.vn[vn1], scn.vn[vn2], rec.u, rec.v);

					
					//vec3 rec_normal(scn.trs[ti].n);

					vec3 hit_point(r.o + r.d * rec.t + rec_normal * 0.0002f);
					for (int j = 0; j < num_shadow; ++j)
					{
						int li = scn.Ls.sample();
						vec3 light_point(scn.Ls.sample_light(li));

						vec3 light_direction(light_point - hit_point);

						float length = light_direction.length() * 0.999f;
						Ray shadow_ray(hit_point, light_direction);

						//if (scn.hit_anything_range(n, shadow_ray, length))
						++light_block;
					}
				}
			}
			if (light_block == num_sample * num_shadow)
			{
				visible[j * Width + i] = false;
				++count;
			}
			if (light_block == 0)
			{
				in_light[j * Width + i] = true;
				++count_light;
			}
		}
	}
}

void compute_color_normal(node*& n, Scene& scn, vector<vec3>& color, vector<vec3>& normal)
{
	for (int j = 0; j < Height; ++j)
	{
		#pragma omp parallel for schedule(guided)
		for (int i = 0; i < Width; ++i)
		{
			vec3 sum_color(0.0f);
			vec3 sum_norm(0.0f);

			for (int s = 0; s < 20; ++s)
			{
				float p = ((float)i + randf()) * iWidth;
				float q = ((float)j + randf()) * iHeight;

				p = (2.0f * p - 1.0f) * aspect_ratio * tan_theta;
				q = (1.0f - 2.0f * q) * tan_theta;

				Ray new_ray(look_from, (u * p + v * q - w).norm());

				vec3 col(0.0f);
				vec3 norm(0.0f);

				scn.ray_casting_color_normal(n, new_ray, col, norm);

				sum_color += col;
				sum_norm += norm;
			}

			int ind = j * Width + i;
			normal[ind] = sum_norm * 0.05f;
			color[ind] = sum_color * 0.05f;
		}
	}
}

float compute_shadow_order(node*& n, Scene& scn, const int& ind)
{
	clock_t t = clock();

	for (int j = 0; j < Height; ++j)
	{
		#pragma omp parallel for schedule(guided)
		for (int i = 0; i < Width; ++i)
		{
			for (int s = 0; s < 10; ++s)
			{
				float p = ((float)i + randf()) * iWidth;
				float q = ((float)j + randf()) * iHeight;

				p = (2.0f * p - 1.0f) * aspect_ratio * tan_theta;
				q = (1.0f - 2.0f * q) * tan_theta;

				Ray new_ray(look_from, (u * p + v * q - w).norm());

				HitRecord rec;

				if (scn.all_hit(n, new_ray, rec))
				{
					vec3 hit_point(new_ray.o + new_ray.d * rec.t + rec.n * 0.0002f);
					for (int t = 0; t < 5; ++t)
					{
						int li = scn.Ls.sample();

						vec3 light_point(scn.Ls.sample_light(li));

						vec3 light_direction(light_point - hit_point);

						float length = light_direction.length();

						Ray light_ray(hit_point, light_direction);

						scn.hit_anything_range_test(n, light_ray, length * 0.999f, ind);
					}
				}
			}
		}
	}

	t = (clock() - t) * 0.001f;

	return t;
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
/*
static vec3 clamp_vector(vec3& v, const float& low, const float& high)
{
	float x = min(max(low, v.x), high);
	float y = min(max(low, v.y), high);
	float z = min(max(low, v.z), high);

	return vec3(x, y, z);
}
*/
//Ghi chu
//su dung prev_normal lam cho nguon sang bi xam nhin ko that
//prev_normal su dung de tinh pdf light
//khong nen su dung
/*
vec3 fix_camera_scale(Scene& scn, vec3& look_from)
{
	vec3 v_min = scn.v_min;
	vec3 v_max = scn.v_max;

	vec3 box_extend = v_max - v_min;

	vec3 inv_dimension(1.0f / box_extend.x, 1.0f / box_extend.y, 1.0f / box_extend.z);

	vec3 scale = vec3(65535.0f) * inv_dimension;

	vec3 vertex_scale = look_from * scale; //(look_from - v_min) * scale;

	//vec3 look_from_scale = 65535.0f * (look_from - v_min) / (v_max - v_min);

	//return look_from_scale;

	return vertex_scale;
}
*/
void main()
{
	//cout << "How many Samples would you like: ";
	//int num_sample_render = 32;

	//cin >> num_sample_render;

	//ns = num_sample_render;
	//sqrt_ns = sqrt14(ns);
	//isqrt_ns = 1 / sqrt_ns;

	//kitchen
	//look_from.x -= 0.4f;
	//look_from.z -= 0.4f;

	//salle de bain full wall
	//final
	//look_from -= view_direction * 3000.1f;

	//look_from.y += 3.0f;
	//look_from.x += 5.0f;
	//look_from.z -= 2.0f;
	//look_from -= view_direction * 14.1f;

	//Sang
	//victorian house
	//look_from += view_direction * 3.6f;
	//look_from.x -= 2.7f;

	//coffe maker
	//look_from -= view_direction * 0.12f;
	//look_from.x += 0.020f;

	/*
	look_from -= view_direction * 0.12f;
	look_from.x += 0.069f;

	look_from += view_direction * 0.29f;

	look_from.y += 0.12f;
	*/

	//look_from += view_direction * 0.06f;

	//look_from -= view_direction * 0.05f;

	//look_from += view_direction * 0.1f;

	//Stair Case
	//look_from += view_direction * 0.12f;

	//look_from += view_direction * 0.2f;

	//look_from += view_direction * 2700.0f;

	//look_from -= view_direction * 12.0f;

	look_from += view_direction * 4000.0f;

	//look_from.x += 20;

	//look_from -= view_direction * 300.0f;

	//look_from.x += 1000;

	clock_t t1;
	t1 = clock();

	//----------------------------------BEDROOM----------------------------

	Scene scn("D:\\AutoDesk\\Model\\9701_Download Free 3D Interior Bedroom Model By Lee Nguyen\\maps\\BedRoom.obj");
	
	//scn.Read_Sheen("SheenAlbedo.png");

	//scn.mats[340].Tr = 1.0f;

	/*for (int i = 0; i < scn.mats.size(); ++i)
	{
		if (scn.mats[i].Tr > 0.0f && scn.mats[i].Tr != 0.8f)
			scn.mats[i].Tr == 0.0f;
			//cout << scn.mats[i].Matname << "\n";
	}*/
		//	if (scn.mats[i].type == Back_Ground_type)
	//		cout << scn.mats[i].Matname << "\n";

	t1 = clock() - t1;
	cout << "Reading time: " << t1 / 1000.0f << "s" << "\n";

	
	//cout<<"Normal size: " << scn.vn.size() << "\n";

	/*for (int i = 0; i < scn.trs.size(); ++i)
	{
		cout << scn.trs[i].vn[0] << " " << scn.trs[i].vn[1] << " " << scn.trs[i].vn[2] << "\n";
	}*/

	clock_t t_build;
	t_build = clock();

	cout << sizeof(node) << "\n";
	cout << scn.trs.size() << "\n";
	cout << "v  size: " << scn.v.size()  << "\n";
	cout << "vt size: " << scn.vt.size() << "\n";
	cout << "vn size: " << scn.vn.size() << "\n";
	node* n;

	
	scn.build_bvh(n, 6);

	//scn.delete_before_render();
	//scn.n = n;


	t_build = clock() - t_build;
	cout << "Building time: " << t_build / 1000.0f << "s" << "\n";

	int size = Width * Height;

	vector<vec3> c;
	c.resize(size);
	
	vector<int> sample_count;
	sample_count.resize(size);

	vector<vec3> raw_color;
	raw_color.resize(size);

	clock_t t2 = clock();

	
	//int num_sample = 6;
	//int num_shadow = 12;

	//vec3 prev_normal(view_direction);
	omp_set_num_threads(512);
	/*
	if (scn.use_area_light)
	{
		int count = 0;
		int count_light = 0;
		determine_visibility(n, scn, visible, in_light, num_sample, num_shadow, count, count_light);
		//cout << "Num pixels in shadow: " << count << "\n";
		//cout << "Num pixels totally in light: " << count_light << "\n";

		t2 = clock() - t2;
		cout << "Visibility check: " << t2 / 1000.0f << "s" << "\n";

		clock_t t_shadow = clock();

		//int first_node = compute_shadow_order(n, scn,);

		//scn.first_shadow_node = first_node;

		float t_left = compute_shadow_order(n, scn, 0);
		float t_right = compute_shadow_order(n, scn, 1);

		scn.first_shadow_node = t_left < t_right ? 1 : 0;

		scn.second_shadow_node = 1 - scn.first_shadow_node;

		t_shadow = clock() - t_shadow;

		cout << "Compute Shadow node order:" << t_shadow / 1000.0f << "s" << "\n";
	}
	*/
	
	clock_t t_render = clock();
	for (int j = 0; j < Height; ++j)
	{
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", ns, 100.0f * j / (Height - 1));
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < Width; ++i)
		{

			vector<vec2> sample;
			sample.resize(ns);

			new_multi_jitter(sample);

			int num_sample_use = 0;

			//int index = j * Width + i + rand() % Width;
			//int sample_so_far = 0;
			bool converge = false;

			vec3 sum(0.0f);

			float sum_sqr = 0.0f;
			float sum_so_far = 0.0f;

			int pixel_index = j * Width + i;

			float convergence_rate;

			//if (visible[pixel_index])
			//	convergence_rate = 0.001f;
			//else
				convergence_rate = 0.0001f;

			
			for (int s = 0; s < ns; s += step)
			{

				//float color = 0.0f;

				//float p_index = i;
				//float q_index = j;
				for (int num = 0; num < step; ++num)
				{

					float p = ((float)i + sample[s + num].x) * iWidth;
					float q = ((float)j + sample[s + num].y) * iHeight;

					p = (2.0f * p - 1.0f) * aspect_ratio * tan_theta;
					q = (1.0f - 2.0f * q) * tan_theta;

					Ray new_ray(look_from, (u * p + v * q - w).norm());

					//vec3 current_color(scn.sibenik_tracing_artificial_light_normal_Ke_alpha(n, new_ray));
					
					//vec3 current_color(scn.sibenik_tracing_artificial_light_normal_Ke_alpha_Transparent(n, new_ray));

					//vec3 current_color(scn.sibenik_tracing_artificial_light_normal_Ke_alpha_Transparent_to_BackGround(n, new_ray));

					//vec3 current_color(scn.sibenik_tracing_directional_lighting_2(n, new_ray));

					//vec3 current_color(scn.sibenik_tracing_directional_lighting_circle(n, new_ray));

					//vec3 current_color(scn.sibenik_tracing_directional_lighting_sphere(n, new_ray));

					//vec3 current_color(scn.sibenik_tracing_directional_lighting_circle_2(n, new_ray));

					//final good
					//vec3 current_color(scn.sibenik_tracing_directional_lighting_circle_blur(n, new_ray));

					//vec3 current_color(scn.sibenik_tracing_directional_lighting_SSS(n, new_ray));

					//vec3 current_color(scn.sibenik_tracing_directional_lighting_Fall_Off_Mask(n, new_ray));

					//vec3 current_color(scn.sibenik_tracing_directional_lighting_sun_mirror(n, new_ray));
					
					//vec3 current_color(scn.sibenik_tracing_directional_lighting_circle_blur_final(n, new_ray));

					//vec3 current_color(scn.sibenik_tracing_directional_lighting_circle_blur_final_2(n, new_ray));

					//vec3 current_color(scn.cone_light(n, new_ray));

					vec3 current_color(scn.sibenik_tracing_directional_lighting_circle_blur_final_version(n, new_ray));

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
				sample_count[j * Width + i] = num_sample_use;
			}

		}
	}
	
	//string image_name = "BedRoom_transparent_diffuse_curtain_1.4_diffuse_out_gammar_cloth_0.9f_0.7f.ppm";

	//string image_name = "Test_Rem_normal_no_abs_raw_revese_rec_ke_9.ppm";

	//string image_name = "a_Fix_Picture_and_Light_Diffuse_orange_vanh_vang_big_2_check_light_LT.ppm";

	//string image_name = "a_Velvet_cos_i_ipi_sun_light_transparent_diffuse_out.ppm";

	//string image_name = "a_directional_light_Plastic_Table_No_Swap_u_v_64.ppm";

	//string image_name = "Test_Env_Render_view_From_Original.ppm";

	//string image_name = "Cloth_rough_0.06_with_diffuse_big_neg_z_0.2_spec_color_2_T_Directional_Lighting_spec_cos_o.ppm";

	//string image_name = "Cloth_Lambert_Color.ppm";

	//string image_name = "Test_Curtain_Only_return_color_fix_Tr_big.ppm";

	//string image_name = "a_a_a_a_test_den_che_tam_giac_lay_lai_den_white_3_no_blue_no_6_Ke.ppm";

	//string image_name = "a_a_a_a_test_den_tam_giac_upper_lower_middle.ppm";

	//string image_name = "a_a_a_a_view_co_rem_sky_color_white_4_big_rem_Tr_1.0_back_ground_type_mul_4.ppm";

	//string image_name = "a_a_a_a_original_view_final_correct_fix_Tf_check_floor_cua_finalin_big.ppm";

	//string image_name = "a_a_a_a_sun_transparent_final_light_sun_sau_lung_www_3dsmj_comHyndrangea5.ppm";

	//string image_name = "Test_Door_ko_co_T_mul_color_background.ppm";

	//string image_name = "Test_Original_View_Light_Diffuse_As_Diffuse_bo_rem_sun_mis_Back_Ground_scale_1.ppm";
	//string image_name = "Test_pillow_Velvet_1_0_-0.5_back_ground_scale_9.5_spec_0.6_Dark_Pillow_1.6_big_fall_off_0.5.ppm";//Original_Directional_lighting_2_sun_12_8_4_1.0_0_-0.2_no_env_No_more_Tf.ppm";

	//string image_name = "cloth_curtain_rough_plastic_GGX_spec_0.9_sun_brightter_6_blur_0.27_big_book_y_0.2_rough_plastic_back.ppm";

	//string image_name = "Rem_Fall_Off_Mask_floatOut_-cos_incident_2.0.ppm";

	//string image_name = "Picture_Cone_fix_condition_dist_2500_limit_y_2420_bsdf_val_mis_div_pdf_light.ppm";

	//string image_name = "Picture_Pos_2_cos_0.48_limit_y_2620_up_200.ppm";

	//string image_name = "Picture_Light_Cone_cos_view_0.78_only_plus_light_lower_2190_forward_200_dist_2200_up_200_inv_length2.ppm";

	string image_name = "back_800_up_600_cos_view_0.6.ppm";

	ofstream ofs(image_name);

	ofs << "P3\n" << Width << " " << Height << "\n255\n";

	//string raw_name = image_name + "_" + to_string(ns) + "_raw_.txt";

	//ofstream ofs_raw(raw_name);
	//ofs_raw << "P3\n" << Width << " " << Height << "\n255\n";

	

	clock_t t_write = clock();

	for (int i = 0; i < size; ++i)
	{
		//vec3 color = Filmmic_ToneMapping(c[i]);
		//vec3 raw = raw_color[i];
		//ofs_raw << raw.x << " " << raw.y << " " << raw.z << " " << sample_count[i] << "\n";

		vec3 color = Filmic_unchart2(c[i]);
		//vec3 color = reinhard_extended_luminance(c[i], 4.0f);
		//vec3 color = c[i];

		//color = clamp_vector(color, 0.0f, 1.0f);
		color *= 255.99f;

		ofs << color.x << " " << color.y << " " << color.z << "\n";
	}
	
	
	ofs.close();
	
	vector<vec3>().swap(c);
	vector<vec3>().swap(raw_color);
	vector<int>().swap(sample_count);

	delete_bvh(&n);
	scn.delete_struct();

	t_render = clock() - t_render;

	float render_time = ((double)t_render) / CLOCKS_PER_SEC;

	std::cout << "\nRendering time: " << ((double)t_render) / CLOCKS_PER_SEC << "\n";

	t_write = clock() - t_write;

	/*float write_time = ((double)t_write) / CLOCKS_PER_SEC;
	ofs_log << "Render Time: " << render_time << "\n";
	ofs_log << "Write Time: " << write_time << "\n";

	ofstream ofs_log(image_name + "_log.txt");*/

	scn.delete_struct();

	//cout << "Delete used data\n";

	//getchar();
	
}


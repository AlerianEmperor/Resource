#include "Scene.h"
//#include "math.h"
#include "Tone_Mapping.h"
//#include "Filter.h"
#include <time.h>
#include <omp.h>

//http://www.raytracerchallenge.com/bonus/texture-mapping.html

using namespace std;
#define Width 400
#define Height 400
//#define Width 400
//#define Height 250
#define iWidth 1 / Width
#define iHeight 1 / Height
#define aspect_ratio Width / Height
#define golden_ratio 1.61803398875

#define g 1.32471795724474602596
#define a1 1.0f / g
#define a2 1.0f / (g * g)

const int ns = 64;
const int sqrt_ns = sqrt14(ns);
const float isqrt_ns = 1 / sqrt_ns;
#define ins  1 / ns

//bath room

//vec3 look_from(0.35473, 1.024759, -0.604429); //view_direction * 0.6f;
const vec3 look_from(-0.18993, 1.012759, -0.855806); const vec3 view_direction = vec3(-0.907783, -0.019999, -0.418963).norm();

//const vec3 look_from(-1.654151f, 2.002430f, -1.132994f); const vec3 view_direction(-0.496300f, -0.305059f, 0.812789f);

//const vec3 look_from(-0.642925, 0.581257, -1.244025); const vec3 view_direction(-0.591846, -0.501213, 0.631271);

//const vec3 look_from(0.008461, 1.045267, -0.645139);
//const vec3 view_direction = vec3(-0.825316, 0.059964, -0.561477).norm();

//const vec3 look_from(0.127673, 1.100172, -0.726294);
//const vec3 view_direction = vec3(-0.927660, 0.064954, -0.367732).norm();

//const vec3 look_from(0.307766, 1.144862, -0.763209);
//const vec3 view_direction = vec3(-0.943369, -0.023173, -0.330963).norm();

//fireplace room
//const vec3 look_from(5.101118f, 1.083746f, -2.756308f);
//const vec3 view_direction = vec3(-0.93355, -0.004821, 0.358416).norm();

//sibenik
//const vec3 look_from(-16.567905, -12.021462, -0.521549);const vec3 view_direction = vec3(0.984883, 0.164252, 0.055007).norm();

//BedRoom
//const vec3 look_from(-0.18078, -3.332532, 0.975756);
//const vec3 view_direction = vec3(-0.000336, 0.995084, 0.099038).norm();

//Living Room
//const vec3 look_from(2.511603, 2.086595, 6.875006);
//const vec3 view_direction = vec3(-0.426998, -0.223106, -0.876297).norm();

const vec3 w = -view_direction;
const vec3 up = vec3(0, 1, 0);
const vec3 u = up.cross(w).norm();
const vec3 v = w.cross(u);
const float tan_theta = tanf(45.0f * pi / 180.0f);//tanf(43.001194f * pi / 180.0f);
const float asspect_tant = aspect_ratio * tan_theta;
const vec3 prev_normal = -w;

struct vec2
{
	float x;
	float y;

	vec2() {}
	vec2(float x_, float y_) : x(x_), y(y_) {}
};

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

void __fastcall determine_visibility(Scene& scn, vector<bool>& visible, int num_sample, int num_shadow, int& count)
{
	for (int j = 0; j < Height; ++j)
	{
		//#pragma omp for schedule(guided)
		for (int i = 0; i < Width; ++i)
		{
			int count_hit_light = 0;

			for (int s = 0; s < num_sample; ++s)
			{
				HitRecord rec;
				float p = ((float)i + randf()) * iWidth;
				float q = ((float)j + randf()) * iHeight;

				p = (2.0f * p - 1.0f) * aspect_ratio * tan_theta;
				q = (1.0f - 2.0f * q) * tan_theta;
				Ray r(look_from, u * p + v * q - w);

				//bool b = scn.all_hit(n, r, rec);
				bool b = scn.all_hit_flat_nodes(r, rec);
				int mtl = scn.trs[rec.ti].mtl;

				if (b && !scn.mats[mtl].isLight)
				{
					int ti = rec.ti;
					vec3 rec_normal(scn.trs[ti].n);

					vec3 hit_point(r.o + r.d * rec.t + rec_normal * 0.0002f);
					for (int j = 0; j < num_shadow; ++j)
					{
						int li = scn.Ls.sample();
						vec3 light_point(scn.Ls.sample_light(li));

						vec3 light_direction(light_point - hit_point);

						float length = light_direction.length() * 0.999f;
						Ray shadow_ray(hit_point, light_direction);

						//if (scn.hit_anything_range(n, shadow_ray, length))
						if (scn.hit_anything_range_flat_node(shadow_ray, length))
							++count_hit_light;
					}
				}
			}
			if (count_hit_light == num_sample * num_shadow)
			{
				visible[j * Width + i] = false;
				++count;
			}
		}
	}
}

static vec3 denan(const vec3& v)
{
	vec3 temp(v);
	if (temp.x != temp.x)
		temp.x = 0.0f;
	if (temp.y != temp.y)
		temp.y = 0.0f;
	if (temp.z != temp.z)
		temp.z = 0.0f;

	return temp;
}

static float Luminance(const vec3& v)
{
	return 0.2126f * v.x + 0.7152f * v.y + 0.0722 * v.z;
	//return 0.299f * v.x + 0.587f * v.y + 0.114f * v.z;
}

void main()
{
	//look_from += view_direction * 0.6f;
	clock_t t1;
	t1 = clock();
	
	Scene scn("E:\\Models\\blender_to_obj\\bathroom\\blender\\textures\\contemporary_bathroom.obj");

	//Scene scn("E:\\Models\\sibenik\\sibenik.obj");
	//Scene scn("E:\\Models\\Living_Room\\living_room\\textures\\living_room.obj");
	//Scene scn("E:\\Models\\blender\\bedroom\\textures\\BedRoom.obj");
	//Scene scn("E:\\a_Sang_Ray_Tracing\\Advance_Features\\40_Sponza_Palace\\Sponza\\Models\\fireplace_room\\fireplace_room.obj");

	//vec3 view_direction((look_at - look_from).norm());

	for (int i = 0; i < scn.Ls.size; ++i)
	{
		if (view_direction.dot(scn.Ls.normal[i]) >= 0.0f)
			scn.Ls.normal[i] = -scn.Ls.normal[i];
	}

	for (auto &trs : scn.trs)
	{
		if (view_direction.dot(trs.n) >= 0.0f)
			trs.n = -trs.n;
	}

	cout << sizeof(node) << "\n";
	cout << sizeof(FlatNode) << "\n";
	cout << scn.trs.size() << "\n";
	node* n;

	scn.build_bvh(n, 6);
	scn.delete_before_render();
	//delete_tree(n);
	delete_bvh(&n);

	//scn.n = n;
	

	t1 = clock() - t1;
	cout << "Building time: " << t1 / 1000.0f << "s" << "\n";

	int size = Width * Height;

	vector<vec3> c;

	c.resize(size);
	//0.0001f 62.37s
	//0.00001f 62.5s

	//std::ofstream ofs("after_optimized_bvh_forward_hit_th_greater_than_0.0f_tl_less_than_mint_shadow_hit_th_greater_than_0.0f.ppm", std::ios::out | std::ios::binary);//6.793
	
	//no shadow ray in unconverge pixel
	//std::ofstream ofs("Final_3_3_0.001f_0.0001f_converge_1024spp.ppm");
	//std::ofstream ofs("adaptive_convergence_only_64spp.ppm");
	//std::ofstream ofs("sibenik_1600_1000_0.01.ppm");

	//std::ofstream ofs("optimize_first_axis_bathroom.ppm");
	//std::ofstream ofs("hit_color2.ppm");
	

	clock_t t2 = clock();
		
	
	cout << "light size :" << scn.Ls.size << "\n";
	
	int step = 16;

	vector<bool> visible(Width * Height, true);
	//16, 6 -> 11s, 72082 unconverge
	//16, 16 -> 26s, 68926 unconverge
	//12, 4 -> 6s, 75259 unconverge
	//12, 12 -> 15s, 70595s
	//6, 4 -> 81000 +
	//4, 4 ->
	int num_sample = 3;
	int num_shadow = 3;

	vec3 prev_normal(view_direction);
	omp_set_num_threads(128);
	
	int count = 0;
	//determine_visibility(n, scn, visible, num_sample, num_shadow, count);
	determine_visibility(scn, visible, num_sample, num_shadow, count);
	cout << "Num pixels in shadow: " << count << "\n";
	
	t2 = clock() - t2;
	cout << "Visibility check: " << t2 / 1000.0f << "s" << "\n";

	clock_t t_render = clock();
	for (int j = 0; j < Height; ++j)
	{
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", ns, 100.0f * j / (Height - 1));
		#pragma omp parallel for schedule(guided)
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

			float sum_square = 0.0f;
			float sum_so_far = 0.0f;

			int pixel_index = j * Width + i;

			for (int s = 0; s < ns; s += step)
			{
				float convergence_rate;

				if (visible[pixel_index])
					//c = denan(scn.path_tracing(n, new_ray));
					//c = denan(scn.path_tracing_for_noisy(n, new_ray));
					convergence_rate = 0.001f;
				else
					convergence_rate = 0.0001f;
				//c = denan(scn.path_tracing(n, new_ray));

				for (int num = 0; num < step; ++num)
				{
					
					float p = ((float)i + sample[s + num].x) * iWidth;
					float q = ((float)j + sample[s + num].y) * iHeight;

					p = (2.0f * p - 1.0f) * aspect_ratio * tan_theta;
					q = (1.0f - 2.0f * q) * tan_theta;

					Ray new_ray(look_from, u * p + v * q - w);
					
					vec3 c(denan(scn.path_tracing_flat_node(new_ray)));
					
					
					float lux = Luminance(c);

					sum_square += lux * lux;
					sum_so_far += lux;

					sum += c;
				}
				num_sample_use += step;

				float mean_sum = sum_so_far / num_sample_use;

				float variance = (sum_square / num_sample_use - mean_sum * mean_sum) / (num_sample_use - 1);

				vec3 value(sum / num_sample_use);

				if (variance < convergence_rate && value.minc() > 0.2f)
				{
					c[j * Width + i] = sum / num_sample_use;
					converge = true;
					break;
				}
			}
			if (!converge)
				c[j * Width + i] = sum * ins;
				//c[j * Width + i] = sum / number_of_sample;
		}
	}

	//std::ofstream ofs("conductor_light_1024.ppm");
	//std::ofstream ofs("Flat_Nodes.ppm");
	std::ofstream ofs("28_bytes.ppm");
	ofs << "P3\n" << Width << " " << Height << "\n255\n";

	for (int i = 0; i < size; ++i)
	{		
		//vec3 color = c[i];
		//vec3 color = Filmmic_ToneMapping(c[i]);	
		//vec3 color = Filmic_yblein(c[i]);
		//vec3 color = ACES_Tone_Mapping(c[i]);
		//vec3 color = ACESFilm(c[i]);

		//vec3 color = Filmic_yblein(out[i]);

		vec3 color = Filmic_yblein(c[i]);

		ofs << color.x * 255.99 << " " << color.y * 255.99 << " " << color.z * 255.99 << "\n";
	
	}

	ofs.close();//vector<SimpleTriangle>().swap(c);
	vector<vec3>().swap(c);
	//delete_tree(n);
	
	scn.delete_struct();
	t_render = clock() - t_render;
	//double time_taken = ((double)t) / CLOCKS_PER_SEC; // in seconds 

	std::cout << "\nRendering time: " << ((double)t_render) / CLOCKS_PER_SEC << "\n";

	//cout << "no hit: " << count_no_hit << "\n";

	vector<bool>().swap(visible);
	//vector<vec3>().swap(out);
	getchar();
}






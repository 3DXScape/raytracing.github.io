//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "../common/rtweekend.h"

#include "../common/camera.h"
#include "../common/color.h"
#include "aarect.h"
#include "box.h"

#include "hittable_list.h"
#include "material.h"
#include "sphere.h"
#include "triangle.h"
#define STB_IMAGE_IMPLEMENTATION
#include "../common/external/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../common/external/stb_image_write.h"
#define STBI_MSC_SECURE_CRT
#include <iostream>
#include <fstream>
#include <time.h>
#include <cstring>
#include <string>


color ray_color(
    const ray& r,
    const color& background,
    const hittable& world,
    shared_ptr<hittable> lights,
    int depth
) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

    // If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.001, infinity, rec))
        return background;

    scatter_record srec;
    color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, srec))
        return emitted;

    if (srec.is_specular) {
        return srec.attenuation
            * ray_color(srec.specular_ray, background, world, lights, depth - 1);
    }

    auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
    mixture_pdf p(light_ptr, srec.pdf_ptr);
    ray scattered = ray(rec.p, p.generate(), r.time());
    auto pdf_val = p.value(scattered.direction());

    return emitted
        + srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
        * ray_color(scattered, background, world, lights, depth - 1)
        / pdf_val;
}


hittable_list random_scene() {
    hittable_list world;

    auto ground_material = make_shared<lambertian>(color(0.5, 0.5, 0.5));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));
    auto triangle_material = make_shared<lambertian>(color(1.0, 1.0, 0.0));
    world.add(make_shared<triangle>(point3(3, 2.0, 0), point3(3, 0.0, 2.0), point3(3, 1.4, 1.4), triangle_material));
    triangle_material = make_shared<lambertian>(color(0.0, 1.0, 1.0));
    world.add(make_shared<triangle>(point3(2.5, 1.0, 0 + 2.0), point3(2.5, 0.0, 1.0 + 2.0), point3(2.5, 0.7, 0.7 + 2.0), triangle_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}


int main(int argc, char* argv[]) {

    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 4800;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 40;
    const int max_depth = 50;

#define CHANNELS 3

    uint8_t* pixels = new uint8_t[image_width * image_height * CHANNELS];    // World

    auto lights = make_shared<hittable_list>();
    lights->add(make_shared<xz_rect>(213, 343, 227, 332, 554, shared_ptr<material>()));
    lights->add(make_shared<sphere>(point3(5, 5, 5), 90, shared_ptr<material>()));

    //auto world = cornell_box();

    color background(0, 0, 0);
    auto world = random_scene();

    // Camera

    //point3 lookfrom(13, 2, 3);
    point3 lookfrom(13, 2, 10);
    point3 lookat(0,0,0);
    vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.01;// 0.1;

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);

    // Render

    // current date/time based on current system

    struct tm ltm;
    time_t now = time(0);
    localtime_s(&ltm, &now);


    // convert now to string form
    char* dt = ctime(&now);
    char yearBuffer[4];
    char monBuffer[4];
    char mdayBuffer[4];
    char hourBuffer[4];
    char minBuffer[4];
    char secBuffer[4];
    snprintf(yearBuffer, 4, "%02d", ltm.tm_year + 1900);
    snprintf(monBuffer, 4, "%02d", ltm.tm_mon);
    snprintf(mdayBuffer, 4, "%02d", ltm.tm_mday);
    snprintf(hourBuffer, 4, "%02d", ltm.tm_hour);
    snprintf(minBuffer, 4, "%02d", ltm.tm_min);
    snprintf(secBuffer, 4, "%02d", ltm.tm_sec);
    std::string fdt = std::string(yearBuffer) + std::string(monBuffer) + std::string(mdayBuffer) + std::string(hourBuffer) + std::string(minBuffer) + std::string(secBuffer);

    std::cout << "The local date and time is: " << dt + fdt<< std::endl;

    int index = 0;

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0,0,0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (j + random_double()) / (image_height-1);
                ray r = cam.get_ray(u, v);
                //pixel_color += ray_color(r, world, max_depth);
                pixel_color += ray_color(r, background, world, lights, max_depth);
            }
            write_color(pixels, index, pixel_color, samples_per_pixel);
            index += 3;
        }
    }
    std::string filename = "rendered.001." + fdt + ".png";
    stbi_write_png(filename.c_str(), image_width, image_height, 3, pixels, image_width * 3);
    delete[] pixels;
    std::cerr << filename << "\nDone.\n";
    std::cerr << "\nDone.\n";
}
/*


int index = 0;
for (int j = height - 1; j >= 0; --j)
{
    for (int i = 0; i < width; ++i)
    {
        float r = (float)i / (float)width;
        float g = (float)j / (float)height;
        float b = 0.2f;
        int ir = int(255.99 * r);
        int ig = int(255.99 * g);
        int ib = int(255.99 * b);

        pixels[index++] = ir;
        pixels[index++] = ig;
        pixels[index++] = ib;
    }
}

// if CHANNEL_NUM is 4, you can use alpha channel in png
stbi_write_png("stbpng.png", width, height, CHANNEL_NUM, pixels, width * CHANNEL_NUM);

// You have to use 3 comp for complete jpg file. If not, the image will be grayscale or nothing.
stbi_write_jpg("stbjpg3.jpg", width, height, 3, pixels, 100);
delete[] pixels;
*/

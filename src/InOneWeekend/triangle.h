#pragma once
#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "../common/rtweekend.h"
#include "hittable.h"
#include "onb.h"

class triangle : public hittable {
public:
    triangle() {}

    triangle(point3 p0, point3 p1, point3 p2, shared_ptr<material> m)
        : vertex0(p0), vertex1(p1), vertex2(p2), mat_ptr(m) {};

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;
// stolen from xzrect - make work for triangle
    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override {
        // The bounding box must have non-zero width in each dimension, so pad the Y
        // dimension a small amount.
        double minX = fmin(fmin(vertex0.x(), vertex1.x()), vertex2.x());
        double minY = fmin(fmin(vertex0.y(), vertex1.y()), vertex2.y());
        double minZ = fmin(fmin(vertex0.z(), vertex1.z()), vertex2.z());
        double maxX = fmax(fmax(vertex0.x(), vertex1.x()), vertex2.x());
        double maxY = fmax(fmax(vertex0.y(), vertex1.y()), vertex2.y());
        double maxZ = fmax(fmax(vertex0.z(), vertex1.z()), vertex2.z());
        output_box = aabb(point3(minX, minY, minZ), point3(maxX, maxY, maxZ));
        return true;
    }

    // stolen from xzrect - make work for triangle
    virtual double pdf_value(const point3& origin, const vec3& v) const override {
        hit_record rec;
        if (!this->hit(ray(origin, v), 0.001, infinity, rec))
            return 0;

        auto area = 1.0;//(x1 - x0) * (z1 - z0);
        auto distance_squared = rec.t * rec.t * v.length_squared();
        auto cosine = fabs(dot(v, rec.normal) / v.length());

        return distance_squared / (cosine * area);
    }

    // stolen from xzrect - make work for triangle
    virtual vec3 random(const point3& origin) const override {
        double lambda0 = random_double(0.0, 1.0);
        double lambda1 = random_double(0.0, 1.0);
        double lambda2 = random_double(0.0, 1.0);
        double length = std::sqrt(lambda0 * lambda0 + lambda1 * lambda1 + lambda2 * lambda2);
        lambda0 = lambda0 / length;
        lambda1 = lambda1 / length;
        lambda1 = lambda1 / length;
        vec3 v1 = vertex0 - vertex2;
        vec3 v0 = vertex1 - vertex2;
        auto random_point = vertex2 + vertex0 * lambda0 + vertex1 * lambda1 + vertex2 * lambda2;
        //auto random_point = point3(random_double(vertex0.x(), vertex1.x()), k, random_double(vertex0.z(), vertex1.z()));
        return random_point - origin;
    }


public:
    point3 vertex0;
    point3 vertex1;
    point3 vertex2;
    shared_ptr<material> mat_ptr;
};

extern bool RayIntersectsTriangle(vec3 rayOrigin,
    vec3 rayVector,
    const triangle* inTriangle,
    vec3& outIntersectionPoint,
    double& t);
/*
// stolen from sphere - make work for triangle
double triangle::pdf_value(const point3& o, const vec3& v) const 
{
    hit_record rec;
    if (!this->hit(ray(o, v), 0.001, infinity, rec))
    {
        return 0;
    }
    auto cos_theta_max = sqrt(1 - radius * radius / (center - o).length_squared());
    auto solid_angle = 2 * pi * (1 - cos_theta_max);

    return  1 / solid_angle;
}
// stolen from sphere - make work for triangle
vec3 triangle::random(const point3& o) const 
{
    vec3 direction = center - o;
    auto distance_squared = direction.length_squared();
    onb uvw;
    uvw.build_from_w(direction);
    return uvw.local(random_to_sphere(radius, distance_squared));
}
// stolen from sphere - make work for triangle
bool triangle::bounding_box(double time0, double time1, aabb& output_box) const {
    output_box = aabb(
        center - vec3(radius, radius, radius),
        center + vec3(radius, radius, radius));
    return true;
}
*/

bool triangle::hit(const ray& r, double t_min, double t_max, hit_record& rec) const 
{
    vec3 rayOrigin(r.orig);
    vec3 rayVector(r.dir);
    vec3 intersectionPoint(0.0, 0.0, 0.0);
    double t = 0.0;
    bool isHit = RayIntersectsTriangle(rayOrigin, rayVector, this, intersectionPoint, t);
    if (isHit)
    {
        // compute surface normal
        vec3 e0(vertex1 - vertex2);
        vec3 e1(vertex0 - vertex2);
        vec3 normal(cross(e0, e1));
        double normalLength = normal.length();
        // require that we are working in metres so the tolerance is 1/10 mm
        if (fabs(normalLength < 0.0001))
        {
            return false;
        }
        normal = normal / normal.length();
        rec.set_face_normal(r, normal);
        rec.mat_ptr = mat_ptr;
        rec.t = t;
        rec.p = r.at(rec.t);
        return true;
    }
    return false;
}
/*
inline double dot(const vec3 &u, const vec3 &v) 
inline vec3 cross(const vec3 &u, const vec3 &v)
*/

#endif
bool RayIntersectsTriangle(vec3 rayOrigin,
    vec3 rayVector,
    const triangle* inTriangle,
    vec3& outIntersectionPoint,
    double& t)
{
    const double EPSILON = 0.0000001;
    vec3 vertex0 = inTriangle->vertex0;
    vec3 vertex1 = inTriangle->vertex1;
    vec3 vertex2 = inTriangle->vertex2;
    vec3 edge1, edge2;

    edge1 = vertex1 - vertex0;
    edge2 = vertex2 - vertex0;
    vec3 h = cross(rayVector, edge2);//rayVector.crossProduct(edge2);
    double a = dot(edge1, h);
    if (a > -EPSILON && a < EPSILON)
    {
        return false;    // This ray is parallel to this triangle.
    }
    double f = 1.0 / a;
    vec3 s = rayOrigin - vertex0;
    double u = f * dot(s, h);
    if (u < 0.0 || u > 1.0)
    {
        return false;
    }
    vec3 q = cross(s, edge1);
    double v = f * dot(rayVector, q);
    if (v < 0.0 || u + v > 1.0)
    {
        return false;
    }
    // At this stage we can compute t to find out where the intersection point is on the line.
    t = f * dot(edge2, q);
    if (t > EPSILON) // ray intersection
    {
        outIntersectionPoint = rayOrigin + rayVector * t;
        return true;
    }
    else // This means that there is a line intersection but not a ray intersection.
    {
        return false;
    }
}
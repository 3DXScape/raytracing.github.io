#pragma once
#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "../common/rtweekend.h"
#include "hittable.h"

class triangle : public hittable {
public:
    triangle() {}

    triangle(point3 p0, point3 p1, point3 p2, shared_ptr<material> m)
        : vertex0(p0), vertex1(p1), vertex2(p2), mat_ptr(m) {};

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;

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
    const float EPSILON = 0.0000001;
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
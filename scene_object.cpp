/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	//return false;

	Vector3D normal_model = Vector3D(0, 0, 1);
	double sideLength = 0.5;

	// transform everything into model space
	Point3D ray_origin_model = worldToModel * ray.origin;
	Vector3D ray_dir_model = worldToModel * ray.dir;
	double distanceFactor = ray.dir.length() / ray_dir_model.length();
	ray_dir_model.normalize();

	// ray parallel to xy-plane
	if (ray_dir_model[2] == 0)
		return false;
	
	// find intersection in model space
	double t_model = -ray_origin_model[2] / ray_dir_model[2];

	// intersection is behind camera
	if (t_model <= 0)
		return false;

	// intersection point
	Point3D intPoint_model = (ray_origin_model + t_model * ray_dir_model);
	double t_world = t_model * distanceFactor;
	
	// check intersection
	bool is_further_intersection = (ray.intersection.none == false) 
		&& (t_world > ray.intersection.t_value);
	bool is_intersect = 
		(intPoint_model[0] > -sideLength && intPoint_model[0] < sideLength
		&& intPoint_model[1] > -sideLength && intPoint_model[1] < sideLength);
	if (!is_intersect || is_further_intersection)
		return false;

	
	// update intersection only when there is intersection
	ray.intersection.none = false;
	ray.intersection.point = modelToWorld * intPoint_model;
	ray.intersection.normal = transNorm(modelToWorld, normal_model);
	ray.intersection.normal.normalize();
	ray.intersection.t_value = t_world;

	return true;
	
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	Point3D sphere_center = Point3D(0, 0, 0);
	int radius = 1;

	// transform everything into model space
	Point3D ray_origin_model = worldToModel * ray.origin;
	Vector3D ray_dir_model = worldToModel * ray.dir;
	double distanceFactor = ray.dir.length() / ray_dir_model.length();
	ray_dir_model.normalize();
	
	
	// find the shortest distance from ray to sphere_center
	double t1 = (sphere_center - ray_origin_model).dot(ray_dir_model);
	Point3D closest_point = ray_origin_model + t1 * ray_dir_model;
	double d = (closest_point - sphere_center).length();
	double t_model = t1 - sqrt(radius - d*d); // subtract to get the closer solution
	double t_world = t_model * distanceFactor;

	// intersection is behind camera
	if (t_model <= 0)
		return false;
	
	// is the intersection is further than previous intersection
	bool is_further_intersection = (ray.intersection.none == false) && (t_world > ray.intersection.t_value);
	bool is_intersect = (d <= radius);
	if (!is_intersect || is_further_intersection)
		return false;
	
	Point3D intPoint_model = ray_origin_model + t_model * ray_dir_model;
	//std::cout << (intPoint_model - sphere_center).length();

	// update intersection
	ray.intersection.none = false;
	ray.intersection.point = modelToWorld * intPoint_model;
	ray.intersection.normal = transNorm(worldToModel, (intPoint_model - sphere_center));
	ray.intersection.normal.normalize();
	ray.intersection.t_value = t_world;
	
	return true;
}


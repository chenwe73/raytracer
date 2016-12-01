/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include <algorithm>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  
	
	double ka = 1.0;
	double kd = 1.0;
	double ks = 1.0;
	
	Vector3D normal = ray.intersection.normal;
	normal.normalize();

	// light
	Vector3D light_dir = get_position() - ray.intersection.point;
	light_dir.normalize();
	// reflection
	Vector3D reflect_dir = 2 * (light_dir.dot(normal)) * normal - light_dir;
	reflect_dir.normalize();
	// vision
	Vector3D vision_dir = -ray.dir;
	vision_dir.normalize();
	
	//std::cout << light_dir.dot(normal) << " ";
	// ambient
	Colour ambientC = ka * ray.intersection.mat->ambient;
	// diffuse
	Colour diffuseC = kd 
		* std::max(light_dir.dot(normal), 0.0) 
		* ray.intersection.mat->diffuse;
	// specular
	Colour specularC = ks 
		* pow( std::max(reflect_dir.dot(vision_dir), 0.0), ray.intersection.mat->specular_exp)
		* ray.intersection.mat->specular;
	
	ray.col = ambientC + diffuseC + specularC;
	//ray.col = ambientC + diffuseC;
	//ray.col = ray.intersection.mat->diffuse;
	//ray.col.clamp();
}


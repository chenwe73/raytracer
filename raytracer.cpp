/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}
void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
    traverseScene(node,ray,_modelToWorld,_worldToModel);
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray, const Matrix4x4& modelToWorld, const Matrix4x4& worldToModel ) {
	SceneDagNode *childPtr;

	// Applies transformation of the current node to the global
	// transformation matrices.
	Matrix4x4 myModelToWorld = modelToWorld*node->trans;
	Matrix4x4 myWorldToModel = node->invtrans*worldToModel;
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, myWorldToModel, myModelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray, myModelToWorld,myWorldToModel);
		childPtr = childPtr->next;
	}

}

void Raytracer::computeShading( Ray3D& ray ) {
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		// shade camera ray
		curLight->light->shade(ray);

		// shadow
		Vector3D light_dir = curLight->light->get_position() - ray.intersection.point;
		light_dir.normalize();
		Ray3D shadowRay(ray.intersection.point, light_dir);

		traverseScene(_root, shadowRay);
		if (!shadowRay.intersection.none)
			ray.col = ray.intersection.mat->ambient;
		
		curLight = curLight->next;
	}
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}


Colour Raytracer::shadeRay( Ray3D& ray, int depth ) {
	Colour col(0.0, 0.0, 0.0);
	traverseScene(_root, ray); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray);
		col = ray.col;

		// specular reflection
		Colour reflectCol(0.0, 0.0, 0.0);
		double specularity = ray.intersection.mat->specularity;
		if (depth < MAXDEPTH && specularity > 0)
		{
			Ray3D reflectionRay = ray.reflectionRay(-ray.dir);
			// recursive call
			reflectCol = shadeRay(reflectionRay, depth + 1);
		}
		col = (1 - specularity) * col + specularity * reflectCol;
	}
	else
	{ // background color
		col = Colour(0, 0, 0);
	}

	// You'll want to call shadeRay recursively (with a different ray, 
	// of course) here to implement reflection/refraction effects.  

	return col; 
}

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2)/tan(fov*M_PI/360.0);

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);

	// Construct a ray for each pixel.
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			// Sets up ray origin and direction in view space, 
			// image plane is at z = -1.
			Point3D origin(0, 0, 0);
			Point3D imagePlane;
			imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
			imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
			imagePlane[2] = -1;
			
			// anti-aliasing (2*2)
			Colour average(0, 0, 0);
			for (int k = -0; k <= 1; k += 2) {
				for (int h = -0; h <= 1; h += 2) {
					Point3D multiSample;
					double r = 0.5/2; //(double)rand() / (RAND_MAX) / 2; // 0 ~ 0.5
					multiSample[0] = imagePlane[0] + r * k /factor;
					r = 0.5/2; //(double)rand() / (RAND_MAX) / 2; // 0 ~ 0.5
					multiSample[1] = imagePlane[1] + r * h /factor;
					multiSample[2] = -1;

					// TODO: Convert ray to world space and call 
					// shadeRay(ray) to generate pixel colour. 	
					Ray3D ray;
					ray.origin = eye;
					ray.dir = (viewToWorld * (multiSample - origin));
					ray.dir.normalize();

					Colour col = shadeRay(ray, 0);
					average = average + col;
				}
			}
			average = 1.0 * average;

			_rbuffer[i*width+j] = int(average[0]*255);
			_gbuffer[i*width+j] = int(average[1]*255);
			_bbuffer[i*width+j] = int(average[2]*255);
			/*
			if (ray.intersection.none)
				std::cout << ".";
			else
				std::cout << "#";
			*/
		}
		//std::cout << "\n";
	}

	flushPixelBuffer(fileName);
}

void setupScene(Raytracer &raytracer)
{
	//return raytracer;
}

int main(int argc, char* argv[])
{	
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	int width = 200; 
	int height = 200; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, 0, 5);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;
	
	// Defines a material for shading.
	Material gold( Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648), 
			Colour(0.628281, 0.555802, 0.366065), 
			51.2, 0 );
	Material jade( Colour(0.0, 0.0, 0.0), Colour(0.54, 0.89, 0.63), 
			Colour(0.316228, 0.316228, 0.316228), 
			12.8, 0 );
	Material reflective(Colour(0.3, 0.3, 0.3), Colour(0.2, 0.2, 0.5),
		Colour(0.628281, 0.555802, 0.366065),
		51.2, 0.5);

	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 4, 0), 
				Colour(0.9, 0.9, 0.9) ) );

	// sphere
	SceneDagNode* sphere = raytracer.addObject(new UnitSphere(), &gold);
	raytracer.translate(sphere, Vector3D(0, -1, 0));	
	double factor1[3] = { 1.0, 1.0, 1.0 };
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);
	
	SceneDagNode* sphere2 = raytracer.addObject(new UnitSphere(), &reflective);
	raytracer.translate(sphere2, Vector3D(-1, 1, -1));
	
	double f = 4;
	double factor2[3] = { f, f, f };
	// floor
	SceneDagNode* floor = raytracer.addObject(new UnitSquare(), &jade);
	raytracer.translate(floor, Vector3D(0, -2, 0));
	raytracer.rotate(floor, 'x', 90);
	raytracer.scale(floor, Point3D(0, 0, 0), factor2);

	// backWall
	SceneDagNode* backWall = raytracer.addObject(new UnitSquare(), &jade);
	raytracer.translate(backWall, Vector3D(0, 0, -2));
	//raytracer.rotate(backWall, 'x', 90);
	raytracer.scale(backWall, Point3D(0, 0, 0), factor2);

	// leftWall
	SceneDagNode* leftWall = raytracer.addObject(new UnitSquare(), &jade);
	raytracer.translate(leftWall, Vector3D(-2, 0, 0));
	raytracer.rotate(leftWall, 'y', -90);
	raytracer.scale(leftWall, Point3D(0, 0, 0), factor2);
	
	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	
	// Render it from a different point of view.
	Point3D eye2(4, 4, 4);
	Vector3D view2(-1, -1, -1);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");
	
	return 0;
}


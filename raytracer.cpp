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

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray, 
		const Matrix4x4& modelToWorld, const Matrix4x4& worldToModel ) {
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
		Point3D light_pos = curLight->light->get_position();

		int hit = 0;
		for (int i = 0; i < SHADOW_SAMPLE; i++)
		{
			// sample light position
			Point3D light_pos_sample = light_pos 
				+ ((double)rand() / RAND_MAX - 1.0 / 2) * LIGHT_SQAURE_WIDTH * LIGHT_U
				+ ((double)rand() / RAND_MAX - 1.0 / 2) * LIGHT_SQAURE_WIDTH * LIGHT_V;

			// shadow ray
			Ray3D shadowRay;
			shadowRay.origin = ray.intersection.point 
				+ REFLECTION_OFFSET * ray.intersection.normal; // shadow acne
			shadowRay.dir = light_pos_sample - ray.intersection.point;
			double t_light = shadowRay.dir.length();
			shadowRay.dir.normalize();

			// intersect with anything in path
			traverseScene(_root, shadowRay);
			// make sure it's inbetween light and intersection point
			if (!shadowRay.intersection.none && shadowRay.intersection.t_value < t_light)
				hit++;
		}
		if (SHADOW_SAMPLE > 0)
		{
			double r = (double)hit / SHADOW_SAMPLE; // 0 ~ 1
			ray.col = (1 - r) * ray.col + (r)* ray.intersection.mat->ambient;
		}
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
	
	// Don't bother shading if the ray didn't hit anything.
	if (!ray.intersection.none) {
		computeShading(ray);
		col = ray.col;

		Colour reflectCol = Colour(0, 0, 0);
		Colour refractCol = Colour(0, 0, 0);
		double specularity = ray.intersection.mat->specularity;
		double refractivity = ray.intersection.mat->refractivity;
		double glossiness = ray.intersection.mat->glossiness;

		Vector3D normal = ray.intersection.normal;
		normal.normalize();
		Vector3D incident = ray.dir;
		incident.normalize();
		double t = ray.intersection.t_value;

		if (depth < MAXDEPTH)
		{
			if (specularity > MIN_SPECULARITY)
			{
				Vector3D reflect_dir = reflect(normal, -incident);

				int reflectSample;
				if (glossiness > MIN_GLOSSINESS)
					reflectSample = REFLECTION_SAMPLE;
				else
					reflectSample = 1;

				// sample square
				Vector3D reflect_u = reflect_dir.cross(Vector3D(1, 0, 0));
				if (reflect_u.length() == 0)
					reflect_u = reflect_dir.cross(Vector3D(0, 1, 0));
				Vector3D reflect_v = reflect_dir.cross(reflect_u);

				// glossy reflection sampling
				for (int i = 0; i < reflectSample; i++)
				{
					// random square offset
					Vector3D reflect_dir_sample = reflect_dir
						+ ((double)rand() / RAND_MAX - 1.0 / 2) * glossiness * reflect_u
						+ ((double)rand() / RAND_MAX - 1.0 / 2) * glossiness * reflect_v;

					// reflectionRay
					Ray3D reflectRay;
					reflectRay.origin =
						ray.intersection.point + REFLECTION_OFFSET * normal; // acne
					reflectRay.dir = reflect_dir_sample;
					reflectRay.dir.normalize();

					// recursive call
					reflectCol = reflectCol + shadeRay(reflectRay, depth + 1);
				}

				reflectCol = 1.0 / reflectSample * reflectCol;
				col = (1 - specularity) * col + specularity * reflectCol;
			}

			if (refractivity > MIN_REFRACTIVITY)
			{
				const double a0 = 0;
				const double a1 = 0;
				const double a2 = 0.4;
				double cosA = incident.dot(normal);
				double k0;
				double k1;
				double k2;

				Vector3D refract_dir;
				double c;
				if (cosA < 0)
				{
					refract_dir = refract(normal, incident, Nt);
					c = -cosA;
					k0 = 1;
					k1 = 1;
					k2 = 1;
				}
				else
				{
					refract_dir = refract(-normal, incident, 1.0 / Nt);
					c = refract_dir.dot(normal);
					k0 = exp(-a0*t);
					k1 = exp(-a1*t);
					k2 = exp(-a2*t);
				}

				// refractRay
				Ray3D refractRay;
				refractRay.origin =ray.intersection.point + REFLECTION_OFFSET * refract_dir; // acne
				refractRay.dir = refract_dir;
				refractRay.dir.normalize();

				refractCol = shadeRay(refractRay, depth + 1);

				double R0 = ((Nt - 1) * (Nt - 1)) / ((Nt + 1) * (Nt + 1));
				double d = 1 - c;
				double R = R0 + (1 - R0) * (d*d*d*d*d);

				col = (R * reflectCol + (1 - R) * refractCol);
				col[0] = k0 * col[0];
				col[1] = k1 * col[1];
				col[2] = k2 * col[2];
			}
		}
	}
	else
	{ // background color
		col = Colour(0, 0, 0);
	}

	if (depth == 1)
		col.clamp();
	
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
			// antialiasing
			Colour antialias_col(0, 0, 0);
			for (int k = -0; k < ANTIALIASING_SAMPLE; k++)
			{
				for (int l = -0; l < ANTIALIASING_SAMPLE; l++)
				{
					// image plane
					Point3D origin(0, 0, 0);
					Point3D imagePlane;
					double r1 = 0.5; // (double)rand() / RAND_MAX; // for stratified
					double r2 = 0.5; // (double)rand() / RAND_MAX;
					imagePlane[0] = j + (k + r1) / ANTIALIASING_SAMPLE;
					imagePlane[1] = i + (l + r2) / ANTIALIASING_SAMPLE;
					imagePlane[0] = (-double(width) / 2 + imagePlane[0]) / factor;
					imagePlane[1] = (-double(height) / 2 + imagePlane[1]) / factor;
					imagePlane[2] = -1;

					// focal plane
					Point3D focalPlane;
					focalPlane[2] = -FOCAL_PLANE_DIST;
					focalPlane[0] = imagePlane[0] * FOCAL_PLANE_DIST;
					focalPlane[1] = imagePlane[1] * FOCAL_PLANE_DIST;

					// FOV
					Colour fov_col(0, 0, 0);
					for (int m = 0; m < FOV_SAMPLE; m++)
					{
						// sample square
						Point3D origin_sample = origin 
							+ ((double)rand() / RAND_MAX - 1.0 / 2) * LENS_WIDTH * LENS_U
							+ ((double)rand() / RAND_MAX - 1.0 / 2) * LENS_WIDTH * LENS_V;

						// camera ray
						Ray3D ray;
						ray.origin = viewToWorld * origin_sample;
						ray.dir = (viewToWorld * (focalPlane - origin_sample));
						ray.dir.normalize();

						Colour col = shadeRay(ray, 1); // !
						fov_col = fov_col + col;
					}
					fov_col = 1.0 / FOV_SAMPLE * fov_col;
					antialias_col = antialias_col + fov_col;
				}
			}
			antialias_col = 1.0 / (ANTIALIASING_SAMPLE*ANTIALIASING_SAMPLE) * antialias_col;

			_rbuffer[i*width + j] = int(antialias_col[0] * 255);
			_gbuffer[i*width + j] = int(antialias_col[1] * 255);
			_bbuffer[i*width + j] = int(antialias_col[2] * 255);
		}

		double time = (double)(clock() - t) / CLOCKS_PER_SEC;
		t = clock();
		double total_time = (double)(clock() - tt) / CLOCKS_PER_SEC;
		std::cout << "row = " << i << " | " << "row time = " << time 
			<< " | " << "totoal time = " << total_time << "\n";
	}

	flushPixelBuffer(fileName);
}

void setupScene(Raytracer &raytracer)
{
	//return raytracer;
}

int main(int argc, char* argv[])
{
	tt = clock();
	std::cout << "rendering... \n";
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Defines a point light source.
	raytracer.addLightSource( new PointLight(Point3D(0, 4.9, 0), Colour(0.9, 0.9, 0.9) ) );

	// sphere
	double factor1[3] = { 2.5, 2.5, 2.5 };
	SceneDagNode* sphere = raytracer.addObject(new UnitSphere(), &glass);
	raytracer.translate(sphere, Vector3D(5, -2.5, 1));	
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);
	
	SceneDagNode* sphere2 = raytracer.addObject(new UnitSphere(), &mirror);
	raytracer.translate(sphere2, Vector3D(-3, -2.5, -2));
	raytracer.scale(sphere2, Point3D(0, 0, 0), factor1);

	//SceneDagNode* sphere3 = raytracer.addObject(new UnitSphere(), &ruby);
	//raytracer.translate(sphere3, Vector3D(2, 2, -4));
	//raytracer.scale(sphere3, Point3D(0, 0, 0), factor1);
	
	double f = 10;
	double factor2[3] = { 20, 20, 20 };
	// floor
	SceneDagNode* floor = raytracer.addObject(new UnitSquare(), &grey);
	raytracer.translate(floor, Vector3D(0, -f/2, 0));
	raytracer.rotate(floor, 'x', 90);
	raytracer.scale(floor, Point3D(0, 0, 0), factor2);
	// backWall
	SceneDagNode* backWall = raytracer.addObject(new UnitSquare(), &grey);
	raytracer.translate(backWall, Vector3D(0, 0, -f/2));
	raytracer.rotate(backWall, 'z', 0);
	raytracer.scale(backWall, Point3D(0, 0, 0), factor2);
	// leftWall
	SceneDagNode* leftWall = raytracer.addObject(new UnitSquare(), &red);
	raytracer.translate(leftWall, Vector3D(-f / 2 * GOLDEN_RATIO, 0, 0));
	raytracer.rotate(leftWall, 'y', -90);
	raytracer.scale(leftWall, Point3D(0, 0, 0), factor2);
	// rightWall
	SceneDagNode* rightWall = raytracer.addObject(new UnitSquare(), &blue);
	raytracer.translate(rightWall, Vector3D(+f / 2 * GOLDEN_RATIO, 0, 0));
	raytracer.rotate(rightWall, 'y', +90);
	raytracer.scale(rightWall, Point3D(0, 0, 0), factor2);
	// ceiling
	SceneDagNode* ceiling = raytracer.addObject(new UnitSquare(), &grey);
	raytracer.translate(ceiling, Vector3D(0, +f/2, 0));
	raytracer.rotate(ceiling, 'x', -90);
	raytracer.scale(ceiling, Point3D(0, 0, 0), factor2);
	// light plane
	double ligth_size[3] = { LIGHT_SQAURE_WIDTH, LIGHT_SQAURE_WIDTH, LIGHT_SQAURE_WIDTH };
	SceneDagNode* lightPlane = raytracer.addObject(new UnitSquare(), &light_mat);
	raytracer.translate(lightPlane, Vector3D(0, f/2 - 0.01, 0));
	raytracer.rotate(lightPlane, 'x', -90);
	raytracer.scale(lightPlane, Point3D(0, 0, 0), ligth_size);
	
	// Render the scene
	// Camera parameters.
	Point3D eye(0, 0, 15);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 50;
	char fileName[] = "view1.bmp";
	raytracer.render(width, height, eye, view, up, fov, fileName);
	
	// Render it from a different point of view.
	Point3D eye2(4.5, 4.5, 4.5);
	Vector3D view2(-1, -1, -1);
	char fileName2[] = "view2.bmp";
	//raytracer.render(width, height, eye2, view2, up, fov, fileName2);
	
	Point3D eye3(0, 0, -100);
	Vector3D view3(0, 0, 1);
	char fileName3[] = "view3.bmp";
	//raytracer.render(width, height, eye3, view3, up, fov, fileName3);
	
	return 0;
}

// reflection of source_dir at ray.intersection
// source points away from intersection
// a is the reflection sample square side length
Ray3D reflectionRay(Ray3D& ray, double a)
{
	Vector3D source_dir = -ray.dir;
	source_dir.normalize();

	// reflection direction
	Vector3D normal = ray.intersection.normal;
	Vector3D reflect_dir = 2 * (source_dir.dot(normal)) * normal - source_dir;
	reflect_dir.normalize();

	// glossy reflection
	double r1 = (double)rand() / (RAND_MAX)-1.0 / 2; // -1/2 ~+1/2
	double r2 = (double)rand() / (RAND_MAX)-1.0 / 2;
	Vector3D u = reflect_dir.cross(Vector3D(1, 0, 0));
	if (u.length() == 0)
		u = reflect_dir.cross(Vector3D(0, 1, 0));
	Vector3D v = reflect_dir.cross(u);
	Vector3D sample_dir = reflect_dir + a * r1 * u + a * r2 * v;
	sample_dir.normalize();

	// reflectRay
	Ray3D reflectionRay;
	reflectionRay.origin = 
		ray.intersection.point + REFLECTION_OFFSET * ray.intersection.normal; // acne
	reflectionRay.dir = sample_dir;
	reflectionRay.dir.normalize();

	return reflectionRay;
}

Vector3D reflect(Vector3D normal, Vector3D incident)
{
	normal.normalize();
	Vector3D reflect = 2 * (incident.dot(normal)) * normal - incident;
	reflect.normalize();
	return reflect;
}

double totalInternalReflection(double n_dot_i, double nt)
{
	const double n = 1;
	return (1 - (n*n) / (nt*nt) * (1 - n_dot_i * n_dot_i));
}

Vector3D refract(Vector3D normal, Vector3D incident, double nt)
{
	const double n = 1;
	normal.normalize();
	incident.normalize();
	double cosA = normal.dot(incident);
	Vector3D refract = (n / nt) * (incident - (cosA)* normal)
		- sqrt(1 - (n*n) / (nt*nt) * (1 - cosA * cosA)) * normal;
	refract.normalize();

	if (1 - (n*n) / (nt*nt) * (1 - cosA * cosA) < 0)
		std::cout << "!"; // todo

	return refract;
}

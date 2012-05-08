// RaytracingNormal.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <vector>
#include <cmath>
#include <glut.h>

struct geometry {
	geometry() {
		id = id_seed++;
	}
	int id;
	static int id_seed;
};
int geometry::id_seed = 0;

struct point {
	point() {}
	point(double x, double y, double z)
		: x(x), y(y), z(z)
	{
	}
	double x, y, z;
	double squared_length() const {
		return x*x+y*y+z*z;
	}
	point operator-(const point& sec) const {
		return point(x-sec.x,y-sec.y,z-sec.z);
	}
	point operator+(const point& sec) const {
		return point(x+sec.x,y+sec.y,z+sec.z);
	}
	double operator*(const point& sec) const {
		return x*sec.x+y*sec.y+z*sec.z;
	}
	point operator*(const double scalar) const {
		return point(scalar*x,scalar*y,scalar*z);
	}
	void normalize() {
		*this = *this*(1/sqrt(squared_length()));
 	}
};
/*point operator*(const double scalar, const point& vector) const {
	return point(scalar*vector.x,scalar*vector.y,scalar*vector.z);
}*/
struct sphere : geometry {
	point center;
	double radius;
};
//struct direction : point {};
typedef point direction;
//#define direction point

struct ray {
	point start;
	direction dir;
};
struct plane : geometry {
	plane(const point& p1, const point& p2, const point& p3) {
		// p1*N + d = 0
		// p2*N + d = 0
		// p3*N + d = 0
		// (p1+p2)/2*N + d = 0
		
	}
	direction N;
	double d;
};
struct intersection_result {
	intersection_result(){}
	intersection_result(bool intersected, ray result=ray(), int intersected_object_id = 0)
		: intersected(intersected)
		, result(result)
		, intersected_object_id(intersected_object_id)
	{
	}
	ray result;
	bool intersected;
	int intersected_object_id;
	double length;
};
struct not_intersected : intersection_result {
	not_intersected() : intersection_result(false) {}
};
intersection_result intersect(const ray& r, const plane& p) {
	//intersection_result ret;
	//	a*
	return not_intersected();
}
intersection_result intersect(ray& r, const sphere& s) {
	r.dir.normalize();
	intersection_result ret;
	ret.intersected_object_id = s.id;
	// ray: p0 = r.start+r.dir*t
	// sphere: (p0.x - s.center.x)^2 + (p0.y - s.center.y)^2 + (p0.z - s.center.z)^2 = s.radius^2
	// ->
	// 
	// ((r.start+r.dir*t).x-s.center.x)^2+... = s.radius^2
	// (r.start.x+r.dir.x*t-s.center.x)^2+... = s.radius^2
	// (x0+vx*t-xc)^2+... = r^2
	// (x0^2+x0*vx*t-x0*xc +vx*t*x0+vx^2*t^2-vx*t*xc -xc*x0-xc*vx*t+xc^2)+... = r^2
	// t^2*(vx^2)+t*(x0*vx+vx*x0-vx*xc-xc*vx)+(x0^2-x0*xc-xc*x0+xc^2)+... -r^2 = 0
	// t^2*(vx^2)+t*(2*vx*(x0-xc))+(x0^2-2*x0*xc+xc^2)+... -r^2 = 0
	double A = 1;//r.dir.squared_length();
	point from_sphere_to_eye = r.start-s.center;
	double B = 2*(from_sphere_to_eye*r.dir);
	double C = from_sphere_to_eye.squared_length()-s.radius*s.radius;
	// At^2+Bt+C=0
	double d = B*B-4*A*C;
	if(d>0) {
		d = sqrt(d);
		double t1 = (d-B)/(2*A);
		double t2 = (-d-B)/(2*A);
		double t;
		if(t1>0) {
			if(t2 > 0)
				if(t1<t2)
					t = t1;
				else
					t = t2;
			else
				t = t1;
		} else {
			if(t2 > 0)
				t = t2;
			else
				return not_intersected();
		}
		ret.intersected = true;
		ret.result.start = point(r.start+r.dir*t);
		ret.result.dir = (ret.result.start-s.center)*(1/s.radius);
		ret.result.dir.normalize();
		ret.length = t;
	} else {
		ret = not_intersected();
	}
	return ret;
}
struct color {
	color(){}
	color(double r, double g, double b, double a) : r(r), g(g), b(b), a(a) {}
	double r,g,b,a;
};
struct material {
	// TODO
	color col;
	double reflect;
};

/*struct GObject {
	GObject() {
	}
	//geometry* g;
	sphere s;
	material m;
};*/
//std::vector<GObject> GObjects;
#define SNOWMAN_MODE 1
#ifdef SNOWMAN_MODE
	const int OBJECT_NUMBER = 4;
#else
	const int OBJECT_NUMBER = 100;
#endif
sphere spheres[OBJECT_NUMBER];
material materials[OBJECT_NUMBER];
color background_color(0.02,0.04,0.085,1);
//const int raytrace_steps = 2; -> change by M key
int raytrace_steps = 2;
//const bool DOF_ENABLED = false; -> change by D key
bool DOF_ENABLED = false;

ray camera;

void DrawCallback() {
	static double t = 0;
	t += 0.1;
	glClearColor(0,1,0,0);
	glClear(GL_COLOR_BUFFER_BIT);
	ray eye;

	float* buffer = new float[1366*768*3];
	intersection_result* irBuff = new intersection_result[1366*768];
	memset(irBuff,0,1366*768*sizeof(intersection_result));
	//for(double x = -0.5; x<0.5; x+=(1.0/1366.0)*10) {
	//	for(double y = -0.5; y<0.5; y+=(1.0/768.0)*10) {
	int screen_width = 512/1, screen_height = 512/1;
	for(int i = 0; i<screen_width; i++) {
		for(int j = 0; j<screen_height; j++) {
			eye = camera;
			//double x = (i-136/2)/136.0;
			//double y = (i-76/2)/76.0;
			eye.dir.x += 0;
			eye.dir.y += (i-screen_width/2)*512/screen_width;
			eye.dir.z += (j-screen_height/2)*512/screen_height;
			eye.start = eye.start + eye.dir*0.005;
			eye.dir.normalize();
			color col = background_color;
			//for each(sphere& s in spheres) {
			for(int oi = 0; oi<OBJECT_NUMBER; ++oi) {
				sphere& s = spheres[oi];
				intersection_result res = intersect(eye,s);
				if(res.intersected) {
					if(irBuff[i+j*1366].intersected && irBuff[i+j*1366].length>0 && irBuff[i+j*1366].length<res.length) {
						irBuff[i+j*1366] = res;
						continue;
					}
					res.result.dir.normalize();
					double cosRef = -(res.result.dir*eye.dir);
					col.r = cosRef*materials[res.intersected_object_id].col.r;
					col.g = cosRef*materials[res.intersected_object_id].col.g;
					col.b = cosRef*materials[res.intersected_object_id].col.b;
					col.a = cosRef*materials[res.intersected_object_id].col.a;
					
					if(raytrace_steps>1) {
						if(res.length<500) {
							for(int oi2 = 0; oi2<OBJECT_NUMBER; ++oi2) {
								sphere& s2 = spheres[oi2];
								intersection_result res2 = intersect(res.result,s2);
								if(res2.intersected) {
									res2.result.dir.normalize();
									double cosRef2 = -(res.result.dir*eye.dir);
									col.r += cosRef2*cosRef*materials[res2.intersected_object_id].col.r*materials[res.intersected_object_id].reflect;
									col.g += cosRef2*cosRef*materials[res2.intersected_object_id].col.g*materials[res.intersected_object_id].reflect;
									col.b += cosRef2*cosRef*materials[res2.intersected_object_id].col.b*materials[res.intersected_object_id].reflect;
									col.a += cosRef2*cosRef*materials[res2.intersected_object_id].col.a*materials[res.intersected_object_id].reflect;
								} else {
									col.r += cosRef*background_color.r*materials[res.intersected_object_id].reflect;
									col.g += cosRef*background_color.g*materials[res.intersected_object_id].reflect;
									col.b += cosRef*background_color.b*materials[res.intersected_object_id].reflect;
									col.a += cosRef*background_color.a*materials[res.intersected_object_id].reflect;
								}
							}
						}
					}
					irBuff[i+j*1366] = res;
					//break;
				}
			}
			buffer[i*3+j*1366*3] = col.r;
			buffer[i*3+j*1366*3+1] = col.g;
			buffer[i*3+j*1366*3+2] = col.b;
		}
	}
	// DOF
	if(DOF_ENABLED) {
		float* tmp = new float[1366*768*3];
		memcpy(tmp,buffer,1366*768*3*sizeof(float));
		double fullscreencolor_lum = 0;
		for(int i = 0; i<screen_width; i++) {
			for(int j = 0; j<screen_height; j++) {
				if(irBuff[i+j*1366].intersected) {
					int r = irBuff[i+j*1366].length*0.05/512*screen_width;
					int color_num = 1;
					double lum = buffer[i*3+j*1366*3]*0.26+buffer[i*3+j*1366*3+1]*0.58+buffer[i*3+j*1366*3+2]*0.16;
					for(int x = i-r; x<i+r; ++x) {
						if(x<0 || x>=1365)
							continue;
						for(int y = j-r; y<j+r; ++y) {
							double from_center_r_sqr = (j-y)*(j-y)+(i-x)*(i-x);
							if(y<0 || y>=768 || r*r<from_center_r_sqr)
								continue;
							++color_num;
							buffer[i*3+j*1366*3] += tmp[x*3+y*1366*3];
							buffer[i*3+j*1366*3+1] += tmp[x*3+y*1366*3+1];
							buffer[i*3+j*1366*3+2] += tmp[x*3+y*1366*3+2];
							/*//HDR Bloom
							if(lum*lum>from_center_r_sqr) {
								double mul_k = (r*r-from_center_r_sqr)/(r*r);
								buffer[x*3+y*1366*3] += tmp[i*3+j*1366*3]*mul_k;
								buffer[x*3+y*1366*3] /= 1+mul_k;
								buffer[x*3+y*1366*3+1] += tmp[i*3+j*1366*3+1]*mul_k;
								buffer[x*3+y*1366*3+1] /= 1+mul_k;
								buffer[x*3+y*1366*3+2] += tmp[i*3+j*1366*3+2]*mul_k;
								buffer[x*3+y*1366*3+2] /= 1+mul_k;
							}*/
						}
					}
					buffer[i*3+j*1366*3] /= color_num;
					buffer[i*3+j*1366*3+1] /= color_num;
					buffer[i*3+j*1366*3+2] /= color_num;
				}
				{
					fullscreencolor_lum += buffer[i*3+j*1366*3]*0.26+buffer[i*3+j*1366*3+1]*0.58+buffer[i*3+j*1366*3+2]*0.16;
				}
			}
		}
		double midcolor = fullscreencolor_lum/(screen_width*screen_height);
		for(int i = 0; i<screen_width; i++) {
			for(int j = 0; j<screen_height; j++) {
				buffer[i*3+j*1366*3] = buffer[i*3+j*1366*3]*0.2/midcolor;//buffer[i*3+j*1366*3]*(0.1 + 0.01*0.9/midcolor);
				buffer[i*3+j*1366*3+1] = buffer[i*3+j*1366*3+1]*0.2/midcolor;//buffer[i*3+j*1366*3+1]*(0.1 + 0.01*0.9/midcolor);// * 0.5/midcolor;
				buffer[i*3+j*1366*3+2] = buffer[i*3+j*1366*3+2]*0.2/midcolor;//buffer[i*3+j*1366*3+2]*(0.1 + 0.01*0.9/midcolor);// * 0.5/midcolor;
			}
		}
		delete [] tmp;
	}
	// draw
	glDrawPixels(1366,768,GL_RGB,GL_FLOAT,buffer);
	delete [] buffer;
	delete [] irBuff;
	glutSwapBuffers();
}
void IdleCallback() {
}
struct {
	direction moving_dir;
	bool moving;
} player;
void TimerCallback(int i) {
	glutTimerFunc(16,TimerCallback,i);
	DrawCallback();
#ifdef SNOWMAN_MODE
	static direction velocity(0,0,0);
	velocity = velocity + point(sin(rand()%1000/10000.0-0.05),sin(rand()%1000/10000.0-0.05),0);
	spheres[3].center = spheres[3].center + velocity;
#else
	for(int i = 0; i<OBJECT_NUMBER; ++i) {
		spheres[i].center = spheres[i].center + point(rand()%1000/100000.0,rand()%1000/100000.0,rand()%1000/100000.0);
	}
#endif
	background_color.r *= 1.01;
	background_color.g *= 1.01;
	background_color.b *= 1.01;
	if(player.moving) {
		player.moving_dir.normalize();
		camera.start = camera.start + player.moving_dir;
	}
}
void MouseMotionCallback(int x, int y) {
	static int mlx = x;
	static int mly = y;
	int dx = x - mlx;
	int dy = y - mly;
	camera.start = camera.start + point(0,-dx,dy);
	mlx = x;
	mly = y;
}
void KeyboardCallback(unsigned char key,int,int) {
	if(key=='w') {
		player.moving = true;
		player.moving_dir = player.moving_dir + camera.dir;
	} else if(key=='s') {
		player.moving = true;
		player.moving_dir = player.moving_dir - camera.dir;
	} else if('m'==key) {
		if(1==raytrace_steps)
			raytrace_steps = 2;
		else if(2==raytrace_steps)
			raytrace_steps = 1;
	} else if('d'==key) {
		DOF_ENABLED = !DOF_ENABLED;
	}
}
void KeyboardUpCallback(unsigned char key,int,int) {
	if(key=='w') {
		player.moving = false;
		player.moving_dir = point(0,0,0);
	} else if(key=='s') {
		player.moving = false;
		player.moving_dir = point(0,0,0);
	}
}

int _tmain(int argc, char* argv[])
{
	srand(9786342);
	// UnitTest:
	{
		//---------------------------------------------------------+
		point p(rand()%10002,rand()%3235,rand()%rand()%1356);
		p.normalize();
		double l = p.squared_length();
		if(l<0.99999 || l>1.000001)
			perror("Unit test not passed!!! Vector normalization fails!");
		ray r0;
		r0.dir = point(rand()%100,rand()%100,rand()%100);
		r0.start = point(rand()%100,rand()%100,rand()%100);
		r0.dir.normalize();
		r0.start.normalize();
		l = r0.start.squared_length()+r0.dir.squared_length();
		if(l<1.99999 || l>2.000001)
			perror("Unit test not passed!!! Vector normalization fails!");

		//---------------------------------------------------------+		
		ray r;
		r.dir = point(rand()%100,rand()%100,rand()%100);
		r.start = point(rand()%100,rand()%100,rand()%100);
		sphere s;
		s.center = r.start+r.dir*(rand()%100);
		s.radius = rand()%100+2;
		intersection_result res = intersect(r,s);
		if(res.intersected) {
			double l =(res.result.start-s.center).squared_length()-s.radius*s.radius;
			if(l>0.001 || l<-0.001)
				perror("Unit test not passed!!! Ray-sphere intersection fails!");
		}
		
		//---------------------------------------------------------+
	}
	camera.dir = point(200,0,0);
	camera.start = point(0,0,0);
	//
	sphere s;// = go.s;
	material m;// = go.m;
#ifdef SNOWMAN_MODE
	spheres[0].id = 0;
	spheres[0].center = point(200,0,0);
	spheres[0].radius = 20;
	materials[0].col = color(1,0,0,1);
	materials[0].reflect = 0.3;

	spheres[1].id = 1;
	spheres[1].center = point(200,0,15);
	spheres[1].radius = 15;
	materials[1].col = color(0,1,0,1);
	materials[1].reflect = 0.3;

	spheres[2].id = 2;
	spheres[2].center = point(200,0,30);
	spheres[2].radius = 10;
	materials[2].col = color(0,0,1,1);
	materials[2].reflect = 0.3;
	
	spheres[3].id = 3;
	spheres[3].center = point(200,0,30);
	spheres[3].radius = 25;
	materials[3].col = color(0.3,0.6,0.4,1);
	materials[3].reflect = 0.3;
#else
	for(int i = 0; i<OBJECT_NUMBER; ++i) {
		s.radius = rand()%1000/2000.0+rand()%1000/2000.0+0.15;
		s.center = point(rand()%10000/1000.0-5,rand()%10000/1000.0-5,rand()%10000/1000.0-5);
		m.col.r = rand()%500/1000.0+0.05;
		m.col.g = rand()%500/1000.0+0.05;
		m.col.b = rand()%500/1000.0+0.05;
		m.col.a = rand()%500/1000.0+0.05;
		m.reflect = rand()%1000/1000.0;
		s.id = i;
		spheres[i] = s;
		materials[i] = m;
	}
#endif
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition(10,10);
	glutInitWindowSize(512,512);
	glutCreateWindow("Raytracing App [Ivajkin Timofej, 2011 (c)]");
	glutIdleFunc(IdleCallback);
	glutDisplayFunc(DrawCallback);
	glutTimerFunc(100,TimerCallback,0);
	glutMotionFunc(MouseMotionCallback);
	glutKeyboardFunc(KeyboardCallback);
	glutKeyboardUpFunc(KeyboardUpCallback);
	glutMainLoop();
	return 0;
}


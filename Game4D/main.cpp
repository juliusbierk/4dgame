// Include standard headers
#include <stdio.h>
#include <stdlib.h>

// Include GLEW
// #define GLEW_STATIC
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>
GLFWwindow* window;

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
using namespace glm;

#include <shader.hpp>
#include <shader.cpp>

#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <chrono>
#include <algorithm>
#include <string>
#include <fstream>

using namespace std;

double sqr(double x) {
	return x*x;
}

class xyzcolPoint {
public:
	double x[6];
	xyzcolPoint();
	xyzcolPoint(float x_, float y, float z, float c1, float c2, float c3);

	double operator[] (size_t i) { return x[i]; }
	const double operator[] (size_t i) const { return x[i]; }
	bool operator <(const xyzcolPoint &p) const {
		return x[0] < p.x[0] || (x[0] == p.x[0] && x[1] < p.x[1]);
	}
};

xyzcolPoint::xyzcolPoint() {
	x[0] = 0; x[1] = 0; x[2] = 0; x[3] = 0; x[4] = 0; x[5] = 0;
}

xyzcolPoint::xyzcolPoint(float x_, float y, float z, float c1, float c2, float c3) {
	x[0] = x_;
	x[1] = y;
	x[2] = z;
	x[3] = c1;
	x[4] = c2;
	x[5] = c3;
}


double crossdiff2D(const xyzcolPoint &O, const xyzcolPoint &A, const xyzcolPoint &B) {
	return (A.x[0] - O.x[0]) * (B.x[1] - O.x[1]) - (A.x[1] - O.x[1]) * (B.x[0] - O.x[0]);
}

vector<xyzcolPoint> convex_hull(vector<xyzcolPoint> P)
{
	int n = P.size(), k = 0;
	vector<xyzcolPoint> H(2 * n);

	// Sort points
	sort(P.begin(), P.end());

	// Build lower hull
	for (int i = 0; i < n; ++i) {
		while (k >= 2 && crossdiff2D(H[k - 2], H[k - 1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	// Build upper hull
	for (int i = n - 2, t = k + 1; i >= 0; i--) {
		while (k >= t && crossdiff2D(H[k - 2], H[k - 1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	H.resize(k - 1);
	return H;
}

class Point {
public:
	double x[4];
	Point();
	Point(double x_, double y, double z, double w);

	double operator[] (size_t i) { return x[i]; }
	const double operator[] (size_t i) const { return x[i]; }

	double norm();
	double norm2();
};

Point::Point() {
	x[0] = 0; x[1] = 0; x[2] = 0; x[3] = 0;
}

Point::Point(double x_, double y, double z, double w) {
	x[0] = x_;
	x[1] = y;
	x[2] = z;
	x[3] = w;
}

double Point::norm() {
	return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
}

double Point::norm2() {
	return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
}


Point operator*(float x, const Point& p) {
	return Point(x*p[0], x*p[1], x*p[2], x*p[3]);
}

Point operator/(const Point& p, float x) {
	return Point(p[0] / x, p[1] / x, p[2] / x, p[3] / x);
}

double operator*(const Point &p1, const Point& p2) { // Dot product
	return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2] + p1[3] * p2[3];
}

Point operator+(const Point &p1, const Point& p2) { // Sum
	return Point(p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2], p1[3] + p2[3]);
}

Point operator-(const Point &p1, const Point& p2) { // Difference
	return Point(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2], p1[3] - p2[3]);
}

std::ostream &operator<<(std::ostream &os, Point const &p) {
	return os << "(" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << ")";
}

class Cube {
public:
	Point c0;
	Point n[3];
	Point unit_ortho;
	Cube();
	Cube(Point c0_, Point n1_, Point n2_, Point n3_);
	Point col;

	vector<Point> intersection(Point ortho, Point x0);
	void update_ortho();
};

void Cube::update_ortho() {
	unit_ortho = Point(-n[1][3] * n[0][2] * n[2][1] + n[1][2] * n[0][3] * n[2][1] + n[1][3] * n[0][1] * n[2][2] - n[1][1] * n[0][3] * n[2][2] - n[1][2] * n[0][1] * n[2][3] + n[1][1] * n[0][2] * n[2][3],
		n[1][3] * n[0][2] * n[2][0] - n[1][2] * n[0][3] * n[2][0] - n[1][3] * n[0][0] * n[2][2] + n[1][0] * n[0][3] * n[2][2] + n[1][2] * n[0][0] * n[2][3] - n[1][0] * n[0][2] * n[2][3],
		-n[1][3] * n[0][1] * n[2][0] + n[1][1] * n[0][3] * n[2][0] + n[1][3] * n[0][0] * n[2][1] - n[1][0] * n[0][3] * n[2][1] - n[1][1] * n[0][0] * n[2][3] + n[1][0] * n[0][1] * n[2][3],
		n[1][2] * n[0][1] * n[2][0] - n[1][1] * n[0][2] * n[2][0] - n[1][2] * n[0][0] * n[2][1] + n[1][0] * n[0][2] * n[2][1] + n[1][1] * n[0][0] * n[2][2] - n[1][0] * n[0][1] * n[2][2]);
	unit_ortho = unit_ortho / unit_ortho.norm();
}

Cube::Cube() {
	;
}

Cube::Cube(Point c0_, Point n1_, Point n2_, Point n3_) {
	c0 = c0_;
	n[0] = n1_;
	n[1] = n2_;
	n[2] = n3_;
	update_ortho();
}

vector<Point> Cube::intersection(Point ortho, Point x0) {
	vector<Point> corners;

	int choice = -1;
	for (int i = 0; i<3; i++) {
		if (abs(ortho*n[i])>0.001) {
			choice = i;
			break;
		}
	}
	if (choice == -1) {
		return corners;
	}

	double norm = ortho*n[choice];
	double a = (ortho*x0 - ortho*c0) / norm;
	int si = (choice + 1) % 3;
	int ui = (choice + 2) % 3;
	double b = -(ortho*n[si]) / norm;
	double c = -(ortho*n[ui]) / norm;

	double sol;
	vector<double> t;
	vector<double> s;
	vector<double> u;

	if (abs(c)>0.0001) { // s=+-1, t=+-1
		for (int i = -1; i<2; i = i + 2) {
			for (int j = -1; j<2; j = j + 2) {
				sol = (i - a - j*b) / c;
				if (abs(sol) <= 1) {
					t.push_back(i); s.push_back(j); u.push_back(sol);
				}
			}
		}
	}
	if (abs(b)>0.0001) { // u=+-1, t=+-1
		for (int i = -1; i<2; i = i + 2) {
			for (int j = -1; j<2; j = j + 2) {
				sol = (i - a - j*c) / b;
				if (abs(sol) <= 1) {
					t.push_back(i); s.push_back(sol); u.push_back(j);
				}
			}
		}
	}
	for (int i = -1; i<2; i = i + 2) {
		for (int j = -1; j<2; j = j + 2) {
			if (abs(a + i*b + j*c) <= 1) {
				s.push_back(i); u.push_back(j); t.push_back(a + i*b + j*c);
			}
		}
	}

	// Find 4D coordinates of corners
	for (int i = 0; i<s.size(); i++) {
		corners.push_back(c0 + t[i] * n[choice] + s[i] * n[si] + u[i] * n[ui]);
	}

	return corners;
}

class Map {
public:
	vector<Cube> cubes;
};

class Camera {
public:
	//Point front = Point(-1, -1, -1, -1); // should be normed
	//Point up = Point(1, -1, 1, -1); // should be normed
	//Point side = Point(-1, 1, 1, -1); // should be normed
	Point front = Point(0, 0, 1, 0);
	Point up = Point(0, 1, 0, 0);
	Point side = Point(1, 0, 0, 0);

	Point pos = -5*up;
	double focal_dist = 1.0;
	int rot4d = 0;
	double rot4d_speed = 0.0;

	Point orthogonal();

	void rotate(double dtheta, int axis);
	void rotate(double dtheta, Point axis);
	vector<Point> rotateIn3Plane(double dtheta, Point axis3D);
	void rotate4D(double dtheta);
	void update4d(double dt);
};

void Camera::update4d(double dt) {
	double sensitivity = 0.0005;
	double gamma = 0.02;
	rot4d_speed += sensitivity*rot4d;
	rot4d = 0;
	rot4d_speed = rot4d_speed / (1 + gamma*dt);
	if (abs(rot4d_speed*dt) > 0.000001) {
		rotate4D(rot4d_speed*dt);
	}
}

Point Camera::orthogonal() {
	return Point(-side[3] * front[2] * up[1] + side[2] * front[3] * up[1] + side[3] * front[1] * up[2] - side[1] * front[3] * up[2] - side[2] * front[1] * up[3] + side[1] * front[2] * up[3],
		side[3] * front[2] * up[0] - side[2] * front[3] * up[0] - side[3] * front[0] * up[2] + side[0] * front[3] * up[2] + side[2] * front[0] * up[3] - side[0] * front[2] * up[3],
		-side[3] * front[1] * up[0] + side[1] * front[3] * up[0] + side[3] * front[0] * up[1] - side[0] * front[3] * up[1] - side[1] * front[0] * up[3] + side[0] * front[1] * up[3],
		side[2] * front[1] * up[0] - side[1] * front[2] * up[0] - side[2] * front[0] * up[1] + side[0] * front[2] * up[1] + side[1] * front[0] * up[2] - side[0] * front[1] * up[2]);
}

void Camera::rotate(double dtheta, int axis) {
	//axis = 0 is rotatation around side, axis=1 is around up
	Point new_front;
	Point new_side;
	Point new_up;
	// axis in basis {front, side, up}
	if (axis == 0) {
		new_front = cos(dtheta) * front + sin(dtheta) * up;
		new_side = side;
		new_up = -sin(dtheta) * front + cos(dtheta) * up;
	}
	else {
		new_front = cos(dtheta) * front - sin(dtheta) * side;
		new_side = sin(dtheta) * front + cos(dtheta) * side;
		new_up = up;
	}

	front = new_front / new_front.norm();
	side = new_side / new_side.norm();
	up = new_up / new_up.norm();
}

void Camera::rotate(double dtheta, Point axis) {
	//Assume that axis lies inside span(front, side, up).
	//Want to rotate such that axis is fixed and orthogonal is fixed.

	//axis3d is axis in basis {front, side, up}
	Point axis3d = Point(axis*front, axis*side, axis*up, 0);
	vector<Point> rotation_matrix = rotateIn3Plane(dtheta, axis3d);

	Point new_front = rotation_matrix[0][0] * front + rotation_matrix[0][1] * side + rotation_matrix[0][2] * up;
	Point new_side = rotation_matrix[1][0] * front + rotation_matrix[1][1] * side + rotation_matrix[1][2] * up;
	Point new_up = rotation_matrix[2][0] * front + rotation_matrix[2][1] * side + rotation_matrix[2][2] * up;

	front = new_front / new_front.norm();
	side = new_side / new_side.norm();
	up = new_up / new_up.norm();
}

vector<Point> Camera::rotateIn3Plane(double dtheta, Point axis3D) {
	Point cross;
	//Axis 3d is (?,?,?,0)
	//front is (1,0,0), side is (0,1,0), up is (0,0,1)

	Point front3d = Point(1, 0, 0, 0);
	cross.x[0] = 0;
	cross.x[1] = axis3D[2];
	cross.x[2] = -axis3D[1];
	front3d = cos(dtheta)*front3d + sin(dtheta)*cross + (axis3D*front3d)*(1 - cos(dtheta))*axis3D;

	Point side3d = Point(0, 1, 0, 0);
	cross.x[0] = -axis3D[2];
	cross.x[1] = 0;
	cross.x[2] = axis3D[0];
	side3d = cos(dtheta)*side3d + sin(dtheta)*cross + (axis3D*side3d)*(1 - cos(dtheta))*axis3D;

	Point up3d = Point(0, 0, 1, 0);
	cross.x[0] = axis3D[1];
	cross.x[1] = -axis3D[0];
	cross.x[2] = 0;
	up3d = cos(dtheta)*up3d + sin(dtheta)*cross + (axis3D*up3d)*(1 - cos(dtheta))*axis3D;

	vector<Point> rot_mat = { front3d, side3d, up3d };
	return rot_mat;
}

void Camera::rotate4D(double dtheta) {
	//Front and up stay fixed.
	side = cos(dtheta)*side + sin(dtheta)*orthogonal();
}

class Scene {
public:
	Map map;
	Camera cam;
	vector< vector<xyzcolPoint> > render();
};

vector< vector<xyzcolPoint> > Scene::render() {
	Point cam_ortho = cam.orthogonal();
	vector<xyzcolPoint>  draw;
	vector< vector<xyzcolPoint> > draws;
	Point x;
	double t, s, u;
	double div = 1 / (sqr(cam.side*cam.up) - 1);
	double c1, c2;
	double dist;
	for (int i = 0; i<map.cubes.size(); i++) {
		vector<Point> v = map.cubes[i].intersection(cam_ortho, cam.pos);

		draw.clear();
		if (v.size() >= 3) {
			for (int j = 0; j<v.size(); j++) {
				t = cam.front*(v[j] - cam.pos);
				dist = ((t > 0) - (t < 0)) * (v[j] - cam.pos).norm(); // This is wrong: should be calculated as a z-coordinate, orthogonal to viewing plane.
				if (abs(t)<0.0001) {
					continue;
				}
				t = cam.focal_dist / t;

				x = cam.pos + t * (v[j] - cam.pos);

				c1 = ((t>0) - (t<0)) * cam.side*(x - cam.pos - cam.focal_dist*cam.front);
				c2 = ((t>0) - (t<0)) * cam.up*(x - cam.pos - cam.focal_dist*cam.front);
				s = div * (c2 * cam.side*cam.up - c1);
				u = -div * (c1 * cam.side*cam.up - c2);
				draw.push_back(xyzcolPoint(s, u, dist, map.cubes[i].col[0], map.cubes[i].col[1], map.cubes[i].col[2]));
			}
		}
		if (draw.size() >= 3) {
			draw = convex_hull(draw);
			draws.push_back(draw);
		}
	}
	return draws;
}

void load_map(string filename, Scene &scene) {
	double number;
	int state = 0;
	ifstream mapfile(filename);
	Cube cube;
	while (mapfile >> number) {
		if (state < 4) {
			cube.c0.x[state] = number;
		}
		else if (state < 8) {
			cube.n[0].x[state-4] = number;
		}
		else if (state < 12) {
			cube.n[1].x[state - 8] = number;
		}
		else if (state < 16) {
			cube.n[2].x[state - 12] = number;
		}
		else if (state < 19) {
			cube.col.x[state - 16] = number;
		}
		state++;
		if (state == 19) {
			cube.update_ortho();
			scene.map.cubes.push_back(cube);
			state = 0;
		}
	}
}

void generate_cube(Scene &scene) {
	// Add hollow tesseract to scene
	vector<Point> cube_points;
	cube_points.push_back(Point(10, 0, 0, 0));
	cube_points.push_back(Point(0, 10, 0, 0));
	cube_points.push_back(Point(0, 0, 10, 0));
	cube_points.push_back(Point(0, 0, 0, 10));
	Cube cube;

	srand(1);
	for (int i = 0; i<4; i++) {
		cube = Cube(cube_points[i], cube_points[(i + 1) % 4],
			cube_points[(i + 2) % 4], cube_points[(i + 3) % 4]);
		cube.col = Point((100 + rand() % 100) / 255.f, (100 + rand() % 100) / 255.f, (100 + rand() % 100) / 255.f, 0);
		scene.map.cubes.push_back(cube);

		cube = Cube(-1 * cube_points[i], cube_points[(i + 1) % 4],
			cube_points[(i + 2) % 4], cube_points[(i + 3) % 4]);
		cube.col = Point((100 + rand() % 255) / 255.f, (100 + rand() % 100) / 255.f, (100 + rand() % 100) / 255.f, 0);
		scene.map.cubes.push_back(cube);
	}
}

Scene scene;
bool use_mouse = true;
bool fly = true;

void scollfun(GLFWwindow *window, double x, double y) {
	scene.cam.rot4d += y;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
		fly = !fly;
	}

	if (key == GLFW_KEY_E && action == GLFW_PRESS) {
		glfwSetCursorPos(window, 500, 500);
		use_mouse = !use_mouse;
	}
}

void move_pos(Point dx) {
	Point x;
	Point gravity(0, 1, 0, 0);
	int n;
	double s, t, u, f;
	double player_size = 0.1;
	double player_move_size = 1.5;
	double player_height = 5.0;
	bool  res;
	if (!fly) {
		for (int i = 0; i < scene.map.cubes.size(); i++) {
			x = scene.cam.pos - scene.map.cubes[i].c0;
			s = x*scene.map.cubes[i].n[0] / scene.map.cubes[i].n[0].norm2();
			t = x*scene.map.cubes[i].n[1] / scene.map.cubes[i].n[1].norm2();
			u = x*scene.map.cubes[i].n[2] / scene.map.cubes[i].n[2].norm2();
			f = x*scene.map.cubes[i].unit_ortho;
			res = false;
			if (abs(scene.map.cubes[i].unit_ortho * gravity) > 0.1) {
				if (abs(s) < 1 + player_size && abs(t) < 1 + player_size && abs(u) < 1 + player_size && abs(f) < player_height) {
					res = true;
				}
			}
			else {
				if (abs(s) < 1 + player_size && abs(t) < 1 + player_size && abs(u) < 1 + player_size && abs(f) < player_move_size) {
					res = true;
				}
			}
			if (res) {
				f = (f*scene.map.cubes[i].unit_ortho)*dx;
				if (f < 0) {
					dx = dx - (dx*scene.map.cubes[i].unit_ortho)*scene.map.cubes[i].unit_ortho;
				}
			}
		}
	}

	scene.cam.pos = scene.cam.pos + dx;
}


int main(void)
{
	// Initialise GLFW
	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(1024, 768, "Game4D", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible.");
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("TransformVertexShader.vertexshader", "ColorFragmentShader.fragmentshader");

	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");
	glm::mat4 MVP = glm::mat4(1.0f); // glm::ortho(-1.f, 1.f, -1.f, 1.f, 0.01f, 1000.f); 
	/*for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << MVP[i][j] << " ";
		}
		cout << endl;
	}*/
	MVP[0][0] = 1.5; // temp fix
	MVP[1][1] = 1.5; // temp fix for z coord currently being dist
	float zfar = 1000.f;
	float znear = 0.001f;
	MVP[3][2] = -(zfar + znear) / (zfar - znear); // This is needed to have a clip-plane.
	MVP[2][2] /= (zfar - znear);

	// generate_cube(scene);
	load_map("map3d.txt", scene);
	
	GLuint vertexbuffer;
	GLuint colorbuffer;
	glGenBuffers(1, &vertexbuffer);
	glGenBuffers(1, &colorbuffer);
	int triangles;
	int c;

	double xpos, ypos;
	double speed = 1.5 / 60.;
	double dx, dy;
	double sensitivity = 0.015;
	auto t1 = chrono::high_resolution_clock::now();
	auto t2 = chrono::high_resolution_clock::now();
	chrono::duration<double, milli> dt_ms;
	glfwSetCursorPos(window, 500, 500);
	int t = 0;
	double dt = 0;

	glfwSetScrollCallback(window, scollfun);

	vector<GLfloat> g_vertex_buffer_data;
	vector<GLfloat> g_color_buffer_data;
	Point move;
	Point gravity = Point(0, 0.01, 0, 0);
	vector< vector<xyzcolPoint> > draw;
	do {
		// cout << scene.cam.front << endl;
		t2 = chrono::high_resolution_clock::now();
		dt_ms = t2 - t1;
		dt = dt_ms.count();
		t1 = t2;

		// Rotation
		if (use_mouse) {
			glfwGetCursorPos(window, &xpos, &ypos);
			glfwSetCursorPos(window, 500, 500);

			dx = xpos - 500;
			dy = ypos - 500;
			scene.cam.rotate(-sensitivity*dx, 1);
			scene.cam.rotate(sensitivity*dy, 0);
			scene.cam.update4d(dt);
		}

		// Enable/Disable look
		glfwSetKeyCallback(window, key_callback);

		// Movement
		move = Point();
		// Move forward
		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
			move = move + scene.cam.front;
		}
		// Move backward
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
			move = move - scene.cam.front;
		}
		// Strafe right
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
			move = move - scene.cam.side;
		}
		// Strafe left
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
			move = move + scene.cam.side;
		}
		if (move.norm() > 0) {
			move_pos(speed*dt*move / move.norm());
		}
		if (!fly) {
			move_pos(dt*gravity);
		}
		

		// Rendering
		draw = scene.render();

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		triangles = 0;
		g_vertex_buffer_data.clear();
		g_color_buffer_data.clear();
		if (draw.size() > 0) {
			for (int i = 0; i < draw.size(); i++) {
				c = draw[i].size();
				for (int j = 0; j < c - 2; j++) {
					for (int k = 0; k < 3; k++) {
						g_vertex_buffer_data.push_back(draw[i][0][k]);
						g_color_buffer_data.push_back(draw[i][0][k + 3]);
					}
					for (int k = 0; k < 3; k++) {
						g_vertex_buffer_data.push_back(draw[i][j + 1][k]);
						g_color_buffer_data.push_back(draw[i][j + 1][k + 3]);
					}
					for (int k = 0; k < 3; k++) {
						g_vertex_buffer_data.push_back(draw[i][j + 2][k]);
						g_color_buffer_data.push_back(draw[i][j + 2][k + 3]);
					}
					triangles += 1;
				}
			}

			glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
			glBufferData(GL_ARRAY_BUFFER, g_vertex_buffer_data.size()*sizeof(g_vertex_buffer_data[0]), &g_vertex_buffer_data[0], GL_STREAM_DRAW);
			glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
			glBufferData(GL_ARRAY_BUFFER, g_color_buffer_data.size()*sizeof(g_color_buffer_data[0]), &g_color_buffer_data[0], GL_STREAM_DRAW);

			// Use our shader
			glUseProgram(programID);

			// Send our transformation to the currently bound shader, 
			// in the "MVP" uniform
			glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

			// 1rst attribute buffer : vertices
			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
			glVertexAttribPointer(
				0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
				3,                  // size
				GL_FLOAT,           // type
				GL_FALSE,           // normalized?
				0,                  // stride
				(void*)0            // array buffer offset
				);

			// 2nd attribute buffer : colors
			glEnableVertexAttribArray(1);
			glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
			glVertexAttribPointer(
				1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
				3,                                // size
				GL_FLOAT,                         // type
				GL_FALSE,                         // normalized?
				0,                                // stride
				(void*)0                          // array buffer offset
				);

			// Draw the triangle !
			glDrawArrays(GL_TRIANGLES, 0, 3 * triangles);

			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);

		}

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	// Cleanup VBO and shader
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteBuffers(1, &colorbuffer);
	glDeleteProgram(programID);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}


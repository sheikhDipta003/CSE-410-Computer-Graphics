#include <windows.h>
#include <GL/glut.h>
#include <cmath>
#include <iostream>
#include <vector>
#define PI (2*acos(0.0))

using namespace std;

// left, right, up, down, tilt anticlockwise, tilt clockwise
GLfloat cam_rot_incr[6] = {5, -5, 5, -5, -5, 5};    //increment of camera rotation angle (in degrees) per key-press

struct triangle{
    float A = 2;    //maximum arm length
    float a = A;    //current arm length
    float da = 0.1;     //change in arm length per key-press
} Triangle;

struct sphere{
    float R = Triangle.A / sqrt(3);     //maximum possible radius of the sphere
    float r = 0;    //current radius
    float dr = 0.06;   //change in radius per key-press
} Sphere;

struct point
{
	float x, y, z;
};

struct camera{
    struct point position;
    struct point up;
    struct point right;
    struct point look;      // camera position, up direction, right direction, left direction
} Camera;

struct point binaryOperation(struct point p, struct point q, char type){
    switch(type){
    case 'c':       //cross product
        return {p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x};
    case 's':       //summation
        return {p.x + q.x, p.y + q.y, p.z + q.z};
    default:
        return p;
    }
}

struct point unaryOperation(struct point p, char type, float a = 1){
    switch(type){
    case 'n':       //negate
        return {-p.x, -p.y, -p.z};
    case 's':       //scalar multiplication
        if (a != 0) return {p.x * a, p.y * a, p.z * a};
    case 'N':       //normalization
        {
            float len = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
            return {p.x / len, p.y / len, p.z / len};
        }
    default:
        return p;
    }
}

struct point rotate_rodrigues_formula(struct point to_rotate, struct point axis, float angle, bool isOrthogonal=true)
{
    //the angle is expected to be provided in degrees
    axis = unaryOperation(axis, 'N');
    struct point resultCross = binaryOperation(axis, to_rotate, 'c');
    struct point p = unaryOperation( to_rotate, 's', cos(angle * PI / 180) );
    struct point r = unaryOperation(resultCross, 's', sin(angle * PI / 180));

    // if to_rotate and axis are orthogonal, their dot product is 0
	if(isOrthogonal)    return {p.x + r.x, p.y + r.y, p.z + r.z};

    float resultDot = to_rotate.x * axis.x + to_rotate.y * axis.y + to_rotate.z * axis.z;
    struct point q = unaryOperation(axis, 's', (1 - cos(angle * PI / 180)) * resultDot);

    return {p.x + q.x + r.x, p.y + q.y + r.y, p.z + q.z + r.z};
}

void init() {
	Camera.position = {0, 0, 5};
	Camera.up = {0, 1, 0};
	Camera.right = {1, 0, 0};
	Camera.look = binaryOperation(Camera.up, Camera.right, 'c');
	glClearColor(0, 0, 0, 0);
}

void drawAxes()
{
    glPushMatrix();

    glBegin(GL_LINES);
    glColor3f(1, 0, 0);     //x : red
    glVertex3f( 100,0,0);
    glVertex3f(-100,0,0);

    glColor3f(0, 1, 0);     //y : green
    glVertex3f(0,-100,0);
    glVertex3f(0, 100,0);

    glColor3f(0, 0, 1);     //z : blue
    glVertex3f(0,0, 100);
    glVertex3f(0,0,-100);
    glEnd();

    glPopMatrix();
}

void drawTriangle()
{
	glBegin(GL_TRIANGLES);
    glVertex3f(1,0,0);
    glVertex3f(0,1,0);
    glVertex3f(0,0,1);
	glEnd();
}

void drawCylHelp(float height, float radius) {
    float offset = 70.5287794 * PI / 180;

    float x = 0, y = 0, theta = 0, del_theta = 0.1;

    glBegin(GL_QUAD_STRIP);
    theta = -offset / 2;
    while( theta <= offset/2 ) {
        x = radius * cos(theta);
        y = radius * sin(theta);
        glVertex3f(x, y , -height/2);
        glVertex3f(x, y , height/2);
        theta += del_theta;
    }
    glEnd();
}

void drawSphereSegment(float radius)
{
    int N = 30;
    vector<vector<point>> spherePoints(N+1, vector<point>(N+1));

    for (int i = 0; i <= N; i++)
    {
        float theta = (float)i / N * 2 * PI;
        for (int j = 0; j <= N; j++)
        {
            float phi = (float)j / N * PI;
            spherePoints.at(i).at(j) = {radius * sin(phi) * cos(theta), radius * sin(phi) * sin(theta), radius * cos(phi)};
        }
    }

    for (int i = 0; i < N; i++)
    {
        glBegin(GL_QUAD_STRIP);
        for (int j = 0; j <= N; j++)
        {
            glVertex3f(spherePoints[i][j].x, spherePoints[i][j].y, spherePoints[i][j].z);
            glVertex3f(spherePoints[i + 1][j].x, spherePoints[i + 1][j].y, spherePoints[i + 1][j].z);
        }
        glEnd();
    }
}

void applyTransformationSphere(){
    glTranslatef(0,0,Triangle.a);
    drawSphereSegment(Sphere.r);
}

void applyTransformationCyl(){
    glTranslatef(Triangle.a / sqrt(2),0,0);
    drawCylHelp(Triangle.a * sqrt(2), Sphere.r);
}

void drawSphereSegFromVertices(){
    //two red
    glColor3f(0, 0, 1);

    glPushMatrix();
    glRotatef(90,1,0,0);
    applyTransformationSphere();
    glPopMatrix();

    glPushMatrix();
    glRotatef(270,1,0,0);
    applyTransformationSphere();
    glPopMatrix();

    //two green
    glColor3f(0, 1, 0);

    glPushMatrix();
    glRotatef(90,0,1,0);
    glTranslatef(0,0,Triangle.a);
    drawSphereSegment(Sphere.r);
    glPopMatrix();

    glPushMatrix();
    glRotatef(270,0,1,0);
    applyTransformationSphere();
    glPopMatrix();

    // two blue
    glColor3f(1, 0, 0);

    glPushMatrix();
    applyTransformationSphere();
    glPopMatrix();

    glPushMatrix();
    glRotatef(180,0,1,0);
    applyTransformationSphere();
    glPopMatrix();
}

void drawCylinderAlongEdges(){
    glColor3f(1, 1, 0);

    //top 4 edges
    glPushMatrix();
    glRotatef(225,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90,0,0,1);
    glRotatef(225,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    glPushMatrix();
    glRotatef(315,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90,0,0,1);
    glRotatef(315,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    //middle 4 edges
    glPushMatrix();
    glRotatef(90,1,0,0);
    glRotatef(45,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90,1,0,0);
    glRotatef(135,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90,1,0,0);
    glRotatef(225,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90,1,0,0);
    glRotatef(315,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    //bottom 4 edges
    glPushMatrix();
    glRotatef(45,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90,0,0,1);
    glRotatef(45,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    glPushMatrix();
    glRotatef(135,0,1,0);
    applyTransformationCyl();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90,0,0,1);
    glRotatef(135,0,1,0);
    applyTransformationCyl();
    glPopMatrix();
}

void drawOctahedron(){
    int i = 0;
    while(i++ < 4){
        glPushMatrix();
        glColor4ub(0x11, 0xED, 0xEA, 0xFF);
        glRotatef(90*i,0,1,0);
        glTranslatef((Triangle.A - Triangle.a) / 3,
                     (Triangle.A - Triangle.a) / 3,
                     (Triangle.A - Triangle.a) / 3);
        glScaled(Triangle.a,Triangle.a,Triangle.a);
        drawTriangle();
        glPopMatrix();
    }

    glPushMatrix();
    glColor4ub(0xF4, 0x5A, 0xF9, 0xFF);
    glRotatef(180,1,0,1);
    glTranslatef((Triangle.A - Triangle.a) / 3,
                 (Triangle.A - Triangle.a) / 3,
                 (Triangle.A - Triangle.a) / 3);
    glScaled(Triangle.a,Triangle.a,Triangle.a);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glColor4ub(0xF4, 0x5A, 0xF9, 0xFF);
    glRotatef(90,0,1,0);
    glRotatef(180,1,0,1);
    glTranslatef((Triangle.A - Triangle.a) / 3,
                 (Triangle.A - Triangle.a) / 3,
                 (Triangle.A - Triangle.a) / 3);
    glScaled(Triangle.a,Triangle.a,Triangle.a);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glColor4ub(0xF4, 0x5A, 0xF9, 0xFF);
    glRotatef(180,0,1,0);
    glRotatef(180,1,0,1);
    glTranslatef((Triangle.A - Triangle.a) / 3,
                 (Triangle.A - Triangle.a) / 3,
                 (Triangle.A - Triangle.a) / 3);
    glScaled(Triangle.a,Triangle.a,Triangle.a);
    drawTriangle();
    glPopMatrix();

    glPushMatrix();
    glColor4ub(0xF4, 0x5A, 0xF9, 0xFF);
    glRotatef(270,0,1,0);
    glRotatef(180,1,0,1);
    glTranslatef((Triangle.A - Triangle.a) / 3,
                 (Triangle.A - Triangle.a) / 3,
                 (Triangle.A - Triangle.a) / 3);
    glScaled(Triangle.a,Triangle.a,Triangle.a);
    drawTriangle();
    glPopMatrix();
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(Camera.position.x , Camera.position.y, Camera.position.z,
              Camera.position.x + Camera.look.x, Camera.position.y + Camera.look.y, Camera.position.z + Camera.look.z,
              Camera.up.x, Camera.up.y, Camera.up.z);

    drawAxes();
    drawOctahedron();
    drawCylinderAlongEdges();
    drawSphereSegFromVertices();

    glutSwapBuffers();
}

void reshape(int width, int height) {
    if (height == 0) height = 1;
    float aspect = (float)width / height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(60, aspect, 1, 100);
    glMatrixMode(GL_MODELVIEW);
}

void specialKeyListener(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_DOWN:
		Camera.position = binaryOperation(Camera.position, unaryOperation(Camera.look, 'n'), 's');
		break;
	case GLUT_KEY_UP:
		Camera.position = binaryOperation(Camera.position, Camera.look, 's');
		break;
	case GLUT_KEY_RIGHT:
		Camera.position = binaryOperation(Camera.position, Camera.right, 's');
		break;
	case GLUT_KEY_LEFT:
		Camera.position = binaryOperation(Camera.position, unaryOperation(Camera.right, 'n'), 's');
		break;
	case GLUT_KEY_PAGE_UP:
		Camera.position = binaryOperation(Camera.position, Camera.up, 's');
		break;
	case GLUT_KEY_PAGE_DOWN:
		Camera.position = binaryOperation(Camera.position, unaryOperation(Camera.up, 'n'), 's');
		break;
	default:
		break;
	}
	glutPostRedisplay();
}

void keyboardListener(unsigned char key, int x, int y) {
    switch (key) {
    case '1':
		Camera.look = rotate_rodrigues_formula(Camera.look, Camera.up, cam_rot_incr[0]);
		Camera.right = binaryOperation(Camera.look, Camera.up, 'c');
		break;
	case '2':
		Camera.look = rotate_rodrigues_formula(Camera.look, Camera.up, cam_rot_incr[1]);
		Camera.right = binaryOperation(Camera.look, Camera.up, 'c');
		break;
	case '3':
		Camera.look = rotate_rodrigues_formula(Camera.look, Camera.right, cam_rot_incr[2]);
		Camera.up = binaryOperation(Camera.right, Camera.look, 'c');
		break;
	case '4':
		Camera.look = rotate_rodrigues_formula(Camera.look, Camera.right, cam_rot_incr[3]);
		Camera.up = binaryOperation(Camera.right, Camera.look, 'c');
		break;
	case '5':
		Camera.right = rotate_rodrigues_formula(Camera.right, Camera.look, cam_rot_incr[4]);
		Camera.up = binaryOperation(Camera.right, Camera.look, 'c');
		break;
	case '6':
		Camera.right = rotate_rodrigues_formula(Camera.right, Camera.look, cam_rot_incr[5]);
		Camera.up = binaryOperation(Camera.right, Camera.look, 'c');
		break;
    case ',':
        Triangle.a = (Triangle.a < 0) ? 0 : (Triangle.a - Triangle.da);
        Sphere.r = (Triangle.a < 0) ? Sphere.R : (Sphere.r + Sphere.dr);
        break;
    case '.':
        Triangle.a = (Triangle.a > Triangle.A) ? Triangle.A : (Triangle.a + Triangle.da);
        Sphere.r = (Triangle.a > Triangle.A) ? 0 : (Sphere.r - Sphere.dr);
        break;
    }
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
	glutInitWindowSize(800, 600);
	glutInitWindowPosition(100, 100);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
	glutCreateWindow("Magic Cube");

	init();
	glutReshapeFunc(reshape);
	glEnable(GL_DEPTH_TEST);

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);

    glutMainLoop();
    return 0;
}

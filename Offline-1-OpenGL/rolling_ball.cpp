#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <windows.h>
#include <GL/glut.h>
#define PI (2*acos(0.0))
#include <cmath>
#include <iostream>
#include <queue>
#include <vector>
#include <chrono>

using namespace std;

GLfloat r = 0.8;        //radius of the ball
GLfloat Theta = 0;      //rotation angle around z-axis in radians   (ball)
GLfloat speed = 0.2;    //linear speed of the ball over the xy plane
const int N = 4;        //number of walls
bool to_pos_x = true;   //'true' if the arrow points towards the positive x-axis
bool eventBasedSimul = false;
bool timeBasedSimul = false;
int CC = 0;             // total collision count

// Wall positions, insert the coordinates of the 4 corners of a 2d wall that is perpendicular to xy plane
GLfloat wallPosition[][3] = {
    {6,6,2},{6,-6,2},{6,-6,0},{6,6,0},  //right
    {6,6,2},{-6,6,2},{-6,6,0},{6,6,0},      //front
    {-6,6,2},{-6,-6,2},{-6,-6,0},{-6,6,0},  //left
    {-6,-6,2},{6,-6,2},{6,-6,0},{-6,-6,0}   //back
};

// left, right, up, down, tilt anticlockwise, tilt clockwise
GLfloat cam_rot_incr[6] = {5, -5, 5, -5, -5, 5};    //increment of camera rotation angle per key-press
GLfloat del_Theta = 0.1;    //change of the angle the position vector of the ball's center makes with +x axis

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

struct point ballCenter = {0, 0};

struct event {
    long long time_to_collide;
    int cc;     //number of collisions expected upto 'time_to_collide'
    string wall_to_collide;
};

struct CompareEvent {
    bool operator()(const event& a, const event& b) const {
        return a.time_to_collide > b.time_to_collide;
    }
};

std::priority_queue<event, std::vector<event>, CompareEvent> eventQueue;

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

struct point unaryOperation(struct point p, char type, float a = 1.0){
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

struct point rotate_rodrigues_formula(struct point to_rotate, struct point axis, double angle, bool isOrthogonal=true)
{
    //the angle is expected to be provided in degrees
    axis = unaryOperation(axis, 'N');
    struct point resultCross = binaryOperation(axis, to_rotate, 'c');
    struct point p = unaryOperation( to_rotate, 's', cos(angle * PI / 180.0) );
    struct point r = unaryOperation(resultCross, 's', sin(angle * PI / 180.0));

    // if to_rotate and axis are orthogonal, their dot product is 0
	if(isOrthogonal)    return {p.x + r.x, p.y + r.y, p.z + r.z};

    float resultDot = to_rotate.x * axis.x + to_rotate.y * axis.y + to_rotate.z * axis.z;
    struct point q = unaryOperation(axis, 's', (1 - cos(angle * PI / 180.0)) * resultDot);

    return {p.x + q.x + r.x, p.y + q.y + r.y, p.z + q.z + r.z};
}

struct point reflect(struct point a, struct point n) {
    // 'a' -> incident vector, 'n' -> normal vector
    // 'n' must be a unit vector
    float a_dot_n = a.x * n.x + a.y * n.y + a.z * n.z;
    return {a.x - 2 * a_dot_n * n.x, a.y - 2 * a_dot_n * n.y, a.z - 2 * a_dot_n * n.z};
}

void drawAxes()
{
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);{
        glVertex3f( 100,0,0);
        glVertex3f(-100,0,0);

        glVertex3f(0,-100,0);
        glVertex3f(0, 100,0);

        glVertex3f(0,0, 100);
        glVertex3f(0,0,-100);
    }glEnd();
}

void drawArrow() {
    glPushMatrix();
    glTranslatef(ballCenter.x, ballCenter.y, r);
    glRotatef(Theta * 180.0 / PI, 0, 0, 1);

    struct point A, B, C;   //arrowhead : starting at the tip, anticlockwise

    if (to_pos_x) {
        A = { 2 * r, 0, 0 };
        B = { 1.3 * r, 0.2 * r, 0 };
        C = { 1.3 * r, -0.2 * r, 0 };
    }
    else {
        A = { -2 * r, 0, 0 };
        B = { -1.3 * r, 0.2 * r, 0 };
        C = { -1.3 * r, -0.2 * r, 0 };
    }

    glBegin(GL_LINES);
    glColor3f(0, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(A.x, A.y, A.z);
    glEnd();

    glBegin(GL_TRIANGLES);
    glVertex3f(A.x, A.y, A.z);
    glVertex3f(B.x, B.y, B.z);
    glVertex3f(C.x, C.y, C.z);
    glEnd();

    glPopMatrix();
}

void drawFloor() {
    int N = 10;     //number of squares
    GLfloat a = 2.0;

    for (int i = -N / 2; i < N / 2; i++) {
        for (int j = -N / 2; j < N / 2; j++) {
            if ((i + j) % 2 == 0) {
                glColor3f(0, 0, 0);
            } else {
                glColor3f(1, 1, 1);
            }

            glBegin(GL_QUADS);
            glVertex2f(i * a, j * a);
            glVertex2f((i + 1) * a, j * a);
            glVertex2f((i + 1) * a, (j + 1) * a);
            glVertex2f(i * a, (j + 1) * a);
            glEnd();
        }
    }
}

void drawBall() {
    int S = 30, T = 30;     // S -> sectors, T -> stacks
    GLfloat del_theta = 2 * PI / S;
    GLfloat del_phi = PI / T;

    glPushMatrix();
    glTranslatef(ballCenter.x, ballCenter.y, 0.0);

    for (int i = 0; i < T; i++) {
        GLfloat phi = i * del_phi - PI / 2;

        for (int j = 0; j <= S; j++) {
            GLfloat theta = j * del_theta;

            struct point bottom_left = {
                    r * cos(phi) * cos(theta),
                    r * cos(phi) * sin(theta),
                    r * sin(phi)
            };
            struct point top_left = {
                    r * cos(phi + del_phi) * cos(theta),
                    r * cos(phi + del_phi) * sin(theta),
                    r * sin(phi + del_phi)
            };
            struct point top_right = {
                    r * cos(phi + del_phi) * cos(theta + del_theta),
                    r * cos(phi + del_phi) * sin(theta + del_theta),
                    r * sin(phi + del_phi)
            };
            struct point bottom_right = {
                    r * cos(phi) * cos(theta + del_theta),
                    r * cos(phi) * sin(theta + del_theta),
                    r * sin(phi)
            };

            if ((i+j) % 2 == 0) {
                glColor3f(1, 0, 0);
            } else {
                glColor3f(0, 1, 0);
            }

            glBegin(GL_QUADS);
            glVertex3f(bottom_left.x, bottom_left.y, bottom_left.z + r);
            glVertex3f(top_left.x, top_left.y, top_left.z + r);
            glVertex3f(top_right.x, top_right.y, top_right.z + r);
            glVertex3f(bottom_right.x, bottom_right.y, bottom_right.z + r);
            glEnd();
        }
    }
    glPopMatrix();
}

void drawWalls(){
    glBegin(GL_QUADS);
    glPushMatrix();
    glColor3f(1, 0, 0);

    for(int i = 0; i < N; i++){
        glVertex3f(wallPosition[4*i][0],wallPosition[4*i][1],wallPosition[4*i][2]);
        glVertex3f(wallPosition[4*i+1][0],wallPosition[4*i+1][1],wallPosition[4*i+1][2]);
        glVertex3f(wallPosition[4*i+2][0],wallPosition[4*i+2][1],wallPosition[4*i+2][2]);
        glVertex3f(wallPosition[4*i+3][0],wallPosition[4*i+3][1],wallPosition[4*i+3][2]);
    }

    glPopMatrix();
    glEnd();
}

// (Q, R) -> endpoints of a line, P -> the point from which the perpendicular distance to this Q---R line is returned
float findDistToLine(struct point P, struct point Q, struct point R) {
    float A[4] = {P.x - Q.x, P.y - Q.y, R.x - Q.x, R.y - Q.y};
    float len_QR = A[2] * A[2] + A[3] * A[3];
    return fabs(A[0] * A[3] - A[1] * A[2]) / sqrt(len_QR);
}

void collisionPredict(){
    float t1 = cos(Theta) == 0.0 ? -1 : (wallPosition[0][1] - r - ballCenter.x) / (speed * cos(Theta));    //right
    float t2 = sin(Theta) == 0.0 ? -1 : (wallPosition[4][1] - r - ballCenter.y) / (speed * sin(Theta));    //front
    float t3 = cos(Theta) == 0.0 ? -1 : (wallPosition[9][1] + r - ballCenter.x) / (speed * cos(Theta));   //left
    float t4 = sin(Theta) == 0.0 ? -1 : (wallPosition[13][1] + r - ballCenter.y) / (speed * sin(Theta));   //back

    if(t1 >= 0) {
        auto currentTime = std::chrono::high_resolution_clock::now();
        long long elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime.time_since_epoch()).count();
        eventQueue.push({elapsedTime + (long long)t1, CC, "right"});
    }
    if(t2 >= 0) {
        auto currentTime = std::chrono::high_resolution_clock::now();
        long long elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime.time_since_epoch()).count();
        eventQueue.push({elapsedTime + (long long)t2, CC, "front"});
    }
    if(t3 >= 0) {
        auto currentTime = std::chrono::high_resolution_clock::now();
        long long elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime.time_since_epoch()).count();
        eventQueue.push({elapsedTime + (long long)t3, CC, "left"});
    }
    if(t4 >= 0) {
        auto currentTime = std::chrono::high_resolution_clock::now();
        long long elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime.time_since_epoch()).count();
        eventQueue.push({elapsedTime + (long long)t4, CC, "back"});
    }
}

void collisionHandleHelp(struct point wallStart, struct point wallEnd){
    struct point ballDirection = {cos(Theta), sin(Theta)};
    bool isCollision = findDistToLine({ ballCenter.x, ballCenter.y }, wallStart, wallEnd) <= r;
    if (isCollision == true) {
        CC++;
        struct point wallNormal = unaryOperation(binaryOperation({ wallEnd.x - wallStart.x, wallEnd.y - wallStart.y, 0.0 }, { 0.0, 0.0, 1.0 }, 'c'), 'N');
        struct point reflected = reflect(ballDirection, wallNormal);
        Theta = atan2(reflected.y, reflected.x);
        ballCenter.x += reflected.x * speed;
        ballCenter.y += reflected.y * speed;
        collisionPredict();
    }
}

void collisionHandle(int i){
    struct point wallStart = { wallPosition[4*i][0], wallPosition[4*i][1] };
    struct point wallEnd = { wallPosition[4*i+1][0], wallPosition[4*i+1][1] };
    collisionHandleHelp(wallStart, wallEnd);
}

void timer(int value){
    if(eventBasedSimul){
        if(!to_pos_x) {
            to_pos_x = true;
        }

        ballCenter.x += speed * cos(Theta);
        ballCenter.y += speed * sin(Theta);

        if(!eventQueue.empty()){
            event impendingEvent = eventQueue.top();
            eventQueue.pop();

            auto currentTime = std::chrono::high_resolution_clock::now();
            long long elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime.time_since_epoch()).count();
            if(impendingEvent.time_to_collide - elapsedTime > 0){
                if(impendingEvent.wall_to_collide == "right" && impendingEvent.cc == CC)
                    glutTimerFunc(impendingEvent.time_to_collide - elapsedTime, collisionHandle, 0);
                else if(impendingEvent.wall_to_collide == "front" && impendingEvent.cc == CC)
                    glutTimerFunc(impendingEvent.time_to_collide - elapsedTime, collisionHandle, 1);
                else if(impendingEvent.wall_to_collide == "left" && impendingEvent.cc == CC)
                    glutTimerFunc(impendingEvent.time_to_collide - elapsedTime, collisionHandle, 2);
                else if(impendingEvent.wall_to_collide == "back" && impendingEvent.cc == CC)
                    glutTimerFunc(impendingEvent.time_to_collide - elapsedTime, collisionHandle, 3);
            }
        }
    }
    glutTimerFunc(16, timer, 0);
}

void timer2(int value){
    if(timeBasedSimul){
        if(!to_pos_x) {
            to_pos_x = true;
        }

        struct point ballDirection = {cos(Theta), sin(Theta)};

        ballCenter.x += speed * cos(Theta);
        ballCenter.y += speed * sin(Theta);

        for (int i = 0; i < N; i++){
            struct point wallStart = { wallPosition[4 * i][0], wallPosition[4 * i][1] };
            struct point wallEnd = { wallPosition[4 * i + 1][0], wallPosition[4 * i + 1][1] };
            bool isCollision = findDistToLine({ ballCenter.x, ballCenter.y }, wallStart, wallEnd) <= r;

            if (isCollision == true) {
                struct point wallNormal = unaryOperation(binaryOperation({ wallEnd.x - wallStart.x, wallEnd.y - wallStart.y, 0.0 }, { 0.0, 0.0, 1.0 }, 'c'), 'N');
                struct point reflected = reflect(ballDirection, wallNormal);

                Theta = atan2(reflected.y, reflected.x);
                ballCenter.x += reflected.x * speed;
                ballCenter.y += reflected.y * speed;

                break;
            }
        }
    }
    glutTimerFunc(16, timer2, 0);
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
    case 'j':
        Theta += del_Theta;
        break;
    case 'l':
        Theta -= del_Theta;
        break;
    case 'i':
        {
            if(!to_pos_x) {
                to_pos_x = true;
                Theta += PI;
            }
            struct point ballDirection = {cos(Theta), sin(Theta)};

            ballCenter.x += speed * cos(Theta);
            ballCenter.y += speed * sin(Theta);

            for (int i = 0; i < N; i++){
                struct point wallStart = { wallPosition[4 * i][0], wallPosition[4 * i][1] };
                struct point wallEnd = { wallPosition[4 * i + 1][0], wallPosition[4 * i + 1][1] };
                bool isCollision = findDistToLine({ ballCenter.x, ballCenter.y }, wallStart, wallEnd) <= r;

                if (isCollision == true) {
                    struct point wallNormal = unaryOperation(binaryOperation({ wallEnd.x - wallStart.x, wallEnd.y - wallStart.y, 0.0 }, { 0.0, 0.0, 1.0 }, 'c'), 'N');
                    struct point reflected = reflect(ballDirection, wallNormal);

                    Theta = atan2(reflected.y, reflected.x);
                    ballCenter.x += reflected.x * speed;
                    ballCenter.y += reflected.y * speed;

                    break;
                }
            }
        }
        break;
    case 'k':
        {
            if(to_pos_x) {
                to_pos_x = false;
                Theta += PI;
            }
            struct point ballDirection = {cos(Theta), sin(Theta)};

            ballCenter.x += speed * cos(Theta);
            ballCenter.y += speed * sin(Theta);

            for (int i = 0; i < N; i++){
                struct point wallStart = { wallPosition[4 * i][0], wallPosition[4 * i][1] };
                struct point wallEnd = { wallPosition[4 * i + 1][0], wallPosition[4 * i + 1][1] };
                bool isCollision = findDistToLine({ ballCenter.x, ballCenter.y }, wallStart, wallEnd) <= r;

                if (isCollision == true) {
                    struct point wallNormal = unaryOperation(binaryOperation({ wallEnd.x - wallStart.x, wallEnd.y - wallStart.y, 0.0 }, { 0.0, 0.0, 1.0 }, 'c'), 'N');
                    struct point reflected = reflect(ballDirection, wallNormal);

                    Theta = atan2(reflected.y, reflected.x);
                    ballCenter.x += reflected.x * speed;
                    ballCenter.y += reflected.y * speed;

                    break;
                }
            }
        }
        break;
    case ' ':
        //switch between manual control and event-driven simulation mode
        eventBasedSimul = !eventBasedSimul;
        break;
    case ',':
        //switch between manual control and event-driven simulation mode
        timeBasedSimul = !timeBasedSimul;
        break;
    default:
        break;
	}
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
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0, 0, 0, 0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(Camera.position.x, Camera.position.y, Camera.position.z, Camera.position.x + Camera.look.x, Camera.position.y + Camera.look.y, Camera.position.z + Camera.look.z, Camera.up.x, Camera.up.y, Camera.up.z);

	drawAxes();
	drawFloor();
    drawWalls();
    drawArrow();
	drawBall();

	glutSwapBuffers();
}

void request_for_repaint() {
	glutPostRedisplay();
}

void init() {
	Camera.position = {0, 0, 10};
	Camera.up = {0, 1, 0};
	Camera.right = {1, 0, 0};
	Camera.look = binaryOperation(Camera.up, Camera.right, 'c');
	glClearColor(0, 0, 0, 0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, 1, 1, 100);
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

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitWindowSize(800, 600);
	glutInitWindowPosition(100, 100);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
	glutCreateWindow("Rotating ball");

	init();
	glutReshapeFunc(reshape);
	glEnable(GL_DEPTH_TEST);

	glutDisplayFunc(display);
	glutIdleFunc(request_for_repaint);
	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutTimerFunc(16, timer, 0);
	glutTimerFunc(16, timer2, 0);

	glutMainLoop();
	return 0;
}

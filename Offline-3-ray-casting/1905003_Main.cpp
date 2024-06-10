#include<bits/stdc++.h>
#include <GL/glut.h>
using namespace std;
#include "1905003_Header.h"
#include <windows.h>

int recursion_level;
vector <Object*> objects;
vector <PointLight*> pointLights;
vector <SpotLight*> spotLights;

GLfloat cam_rot_incr[6] = {0.1, -0.1, 0.1, -0.1, -0.1, 0.1};
GLfloat cam_translate_scale = 4;

struct camera{
    Point position;
    Point up;
    Point right;
    Point look;
    double fovY;
} Camera;

struct window {
    double windowWidth;
    double windowHeight;
} Window;

struct image {
    bitmap_image image_bmp;
    int imageHeight;
    int imageWidth;
    int imgCount;
} Image;

void capture()
{
	for(int i=0;i<Image.imageWidth;i++)
		for(int j=0;j<Image.imageHeight;j++)
			Image.image_bmp.set_pixel(i, j, 0, 0, 0);

	double planeDistance = Window.windowHeight / ( 2.0 * tan( (pi / 180.0) * Camera.fovY / 2 ));

	Point topLeft = Camera.position.sum( Camera.look.scale(planeDistance).sum( Camera.up.scale(Window.windowHeight / 2.0).sum(Camera.right.scale(-Window.windowWidth / 2.0)) ));

	double du = Window.windowWidth / (Image.imageWidth * 1.0);
	double dv = Window.windowHeight / (Image.imageHeight * 1.0);

	topLeft = topLeft.sum( Camera.right.scale(du / 2.0).sum(Camera.up.scale(-dv / 2.0)) );

	int nearest = -1;
	double t, t_min;

	for(int i=0;i<Image.imageWidth;i++)
	{
		for(int j=0;j<Image.imageHeight;j++)
		{
			Point currPixel = topLeft.sum( Camera.right.scale(du * i).sum(Camera.up.scale(-dv * j)) );

			Ray ray_camera(Camera.position, currPixel.sum(Camera.position.scale(-1)));
			Color color;

			t_min = -1;
			nearest = -1;
			for(int i = 0; i < objects.size(); i++)
			{
			    if(objects[i])
				t = objects[i]->intersect(ray_camera,color, 0);
				if(t>0 && (nearest == -1 || t < t_min) ){
                    t_min = t;
                    nearest = i;
				}
			}

			if(nearest != -1)
			{
				color = Color(0,0,0);
				double t = objects[nearest]->intersect(ray_camera, color, 1);

				color.r() = (color.r() < 0) ? 0 : ((color.r() > 1) ? 1: color.r());
				color.g() = (color.g() < 0) ? 0 : ((color.g() > 1) ? 1: color.g());
				color.b() = (color.b() < 0) ? 0 : ((color.b() > 1) ? 1: color.b());

				Image.image_bmp.set_pixel(i, j, 255*color.r(), 255*color.g(), 255*color.b());
			}
		}
	}

	Image.image_bmp.save_image("Output_1"+to_string(Image.imgCount)+".bmp");
	Image.imgCount++;
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
		case '0':
			capture();
			break;
		case '1':
            Camera.look = Camera.look.rotate_rodrigues(Camera.up, cam_rot_incr[0]);
            Camera.right = Camera.look.cross(Camera.up);
            break;
        case '2':
            Camera.look = Camera.look.rotate_rodrigues(Camera.up, cam_rot_incr[1]);
            Camera.right = Camera.look.cross(Camera.up);
            break;
        case '3':
            Camera.look = Camera.look.rotate_rodrigues(Camera.right, cam_rot_incr[2]);
            Camera.up = Camera.right.cross(Camera.look);
            break;
        case '4':
            Camera.look = Camera.look.rotate_rodrigues(Camera.right, cam_rot_incr[3]);
            Camera.up = Camera.right.cross(Camera.look);
            break;
        case '5':
            Camera.right = Camera.right.rotate_rodrigues(Camera.look, cam_rot_incr[4]);
            Camera.up = Camera.right.cross(Camera.look);
            break;
        case '6':
            Camera.right = Camera.right.rotate_rodrigues(Camera.look, cam_rot_incr[5]);
            Camera.up = Camera.right.cross(Camera.look);
            break;
		default:
			break;
	}
}

void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:
            Camera.position = Camera.position.sum(Camera.look.scale(-cam_translate_scale));
            break;
        case GLUT_KEY_UP:
            Camera.position = Camera.position.sum(Camera.look.scale(cam_translate_scale));
            break;
        case GLUT_KEY_RIGHT:
            Camera.position = Camera.position.sum(Camera.right.scale(cam_translate_scale));
            break;
        case GLUT_KEY_LEFT:
            Camera.position = Camera.position.sum(Camera.right.scale(-cam_translate_scale));
            break;
        case GLUT_KEY_PAGE_UP:
            Camera.position = Camera.position.sum(Camera.up.scale(cam_translate_scale));
            break;
        case GLUT_KEY_PAGE_DOWN:
            Camera.position = Camera.position.sum(Camera.up.scale(-cam_translate_scale));
            break;
		default:
			break;
	}
}

void loadData()
{
	ifstream inFile("scene.txt");
    Object *temp;
	int N, plCount, slCount;
    string type;

	inFile >> recursion_level >> Image.imageHeight >> N;
	Image.imageWidth = Image.imageHeight;

	for(int i=0;i<N;i++)
	{
		inFile >> type;
//		temp->setType(type);

		if(type == "sphere"){
			temp = new Sphere();
			inFile >> temp->reference_point.x() >> temp->reference_point.y() >> temp->reference_point.z() >> temp->length >> temp->color.r() >> temp->color.g() >> temp->color.b();

            for(int i = 0; i < 4; i++){
                inFile >> temp->coefficients[i];
            }

            inFile >> temp->shine;
		}
		else if(type == "triangle"){
			temp = new Triangle();
			inFile >> ((Triangle*) temp)->a.x() >> ((Triangle*) temp)->a.y() >> ((Triangle*) temp)->a.z()
                >> ((Triangle*) temp)->b.x() >> ((Triangle*) temp)->b.y() >> ((Triangle*) temp)->b.z()
                >> ((Triangle*) temp)->c.x() >> ((Triangle*) temp)->c.y() >> ((Triangle*) temp)->c.z()
                >> ((Triangle*) temp)->color.r() >> ((Triangle*) temp)->color.g() >> ((Triangle*) temp)->color.b();

            for(int i = 0; i < 4; i++){
                inFile >> ((Triangle*) temp)->coefficients[i];
            }

            inFile >> ((Triangle*) temp)->shine;
		}
		else if(type == "general"){
			temp = new General();
            inFile >> ((General*) temp)->A >> ((General*) temp)->B >> ((General*) temp)->C >> ((General*) temp)->D
                >> ((General*) temp)->E >> ((General*) temp)->F >> ((General*) temp)->G
                >> ((General*) temp)->H >> ((General*) temp)->I >> ((General*) temp)->J
                >> ((General*) temp)->reference_point.x() >> ((General*) temp)->reference_point.y() >> ((General*) temp)->reference_point.z()
                >> ((General*) temp)->length >> ((General*) temp)->width >> ((General*) temp)->height
                >> ((General*) temp)->color.r() >> ((General*) temp)->color.g() >> ((General*) temp)->color.b();

            for(int i = 0; i < 4; i++) {
                inFile >> ((General*) temp)->coefficients[i];
            }

            inFile >> ((General*) temp)->shine;
		}
		objects.push_back(temp);
	}

	inFile >> plCount;

	for(int i=0;i<plCount;i++){
		PointLight *pl = new PointLight();
        inFile >> pl->light_pos.x() >> pl->light_pos.y() >> pl->light_pos.z() >> pl->color.r() >> pl->color.g() >> pl->color.b();
		pointLights.push_back(pl);
	}

	inFile >> slCount;

	for(int i=0;i<slCount;i++){
		SpotLight *sl = new SpotLight();
        inFile >> sl->point_light.light_pos.x() >> sl->point_light.light_pos.y() >> sl->point_light.light_pos.z() >> sl->point_light.color.r() >> sl->point_light.color.g() >> sl->point_light.color.b()
            >> sl->light_dir.x() >> sl->light_dir.y() >> sl->light_dir.z() >> sl->cutoff_angle;
		spotLights.push_back(sl);
	}

	temp = new Floor(1000, 20);
	temp->coefficients[0] = 0.4;
	temp->coefficients[1] = 0.2;
	temp->coefficients[2] = 0.2;
	temp->coefficients[3] = 0.2;
	objects.push_back(temp);
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(Camera.position.x() , Camera.position.y(), Camera.position.z(),
              Camera.position.x() + Camera.look.x(), Camera.position.y() + Camera.look.y(), Camera.position.z() + Camera.look.z(),
              Camera.up.x(), Camera.up.y(), Camera.up.z());


	glMatrixMode(GL_MODELVIEW);

    for (int i=0; i<objects.size(); i++){
		objects[i]->draw();
	}

	for (int i=0; i<pointLights.size(); i++){
		pointLights[i]->draw();
	}

	for(int i=0;i<spotLights.size();i++){
		spotLights[i]->draw();
	}

	glutSwapBuffers();
}

void redisplay(){
	glutPostRedisplay();
}

void init(){
	loadData();
	Image.image_bmp = bitmap_image(Image.imageWidth, Image.imageHeight);

	Camera.position = {120, 70, 25};
	Camera.up = {0, 0, 1};
	Camera.right = {-1 / sqrt(2), 1 / sqrt(2), 0};
	Camera.look = Camera.up.cross(Camera.right);
	Image.imgCount = 1;
	Camera.fovY = 80.0;
	Window.windowHeight = Window.windowWidth = 500;
	glClearColor(0, 0, 0, 0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(80, 1, 1, 1000.0);
}

void clearMem(){
    objects.clear();
	pointLights.clear();
	spotLights.clear();
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

	glutCreateWindow("Assignment 3: Ray Tracing - 1905003");

	init();

	glEnable(GL_DEPTH_TEST);

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutIdleFunc(redisplay);

	glutMainLoop();

	clearMem();

	return 0;
}

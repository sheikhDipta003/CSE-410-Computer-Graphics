#ifndef HEADERS_H_INCLUDED
#define HEADERS_H_INCLUDED

#include <vector>
#include <string>
#include "bitmap_image.hpp"
#include <GL/glut.h>
#define pi (2*acos(0.0))
#define epsilon 1e-5
#define epsilon2 0.001

using namespace std;

class Point{
private:
    double x1, y1, z1, w1;

public:
    Point(){
        this->x1 = 0.0;
        this->y1 = 0.0;
        this->z1 = 0.0;
        this->w1 = 1.0;
    }

    Point(double x1, double y1, double z1){
        this->x1 = x1;
        this->y1 = y1;
        this->z1 = z1;
        this->w1 = 1.0;
    }

    Point(double x1, double y1, double z1, double w1){
        this->x1 = x1;
        this->y1 = y1;
        this->z1 = z1;
        this->w1 = w1;
    }

    Point(const Point &other) : x1(other.x1), y1(other.y1), z1(other.z1), w1(other.w1) {}

    double& x(){ return this->x1; }

    double& y(){ return this->y1; }

    double& z(){ return this->z1; }

    double& w(){ return this->w1; }

    vector<double> getCoords(){
        return {this->x1, this->y1, this->z1, this->w1};
    }

    void scale_w_to_one(){
        this->x1 /= this->w1;
        this->y1 /= this->w1;
        this->z1 /= this->w1;
        this->w1 = 1.0;
    }

    Point sum(Point other){
        return Point(this->x1 + other.x(),
                     this->y1 + other.y(),
                     this->z1 + other.z());
    }

    Point scale(double s){
        return Point(x1*s, y1*s, z1*s);
    }

    double dot(Point other){
        return (this->x1 * other.x() + this->y1 * other.y() + this->z1 * other.z());
    }

    Point cross(Point other){
        return Point(this->y1 * other.z() - this->z1 * other.y(),
                    this->z1 * other.x() - this->x1 * other.z(),
                    this->x1 * other.y() - this->y1 * other.x());
    }

    void norm(){
        double len = sqrt(x1*x1 + y1*y1 + z1*z1);
        this->x1 /= len;
        this->y1 /= len;
        this->z1 /= len;
    }

    double length(){ return sqrt(x1*x1 + y1*y1 + z1*z1); }

    Point rotate_rodrigues(Point axis, double theta_rad){
        axis.norm();
        Point res1 = this->scale( cos(theta_rad) );
        Point res2 = axis.scale( axis.dot(*this) * (1-cos(theta_rad)) );
        Point res3 = axis.cross(*this).scale( sin(theta_rad) );
        return res1.sum(res2).sum(res3);
    }
};

class Color{
private:
    double red, green, blue;

public:
    Color(){
        red = green = blue = 0.0;
    }

    Color(double _red, double _green, double _blue){
        this->red = _red;
        this->green = _green;
        this->blue = _blue;
    }

    double& r(){ return this->red; }
    double& g(){ return this->green; }
    double& b(){ return this->blue; }

    Color& operator=(const Color& other) {
        if (this != &other) {
            red = other.red;
            green = other.green;
            blue = other.blue;
        }
        return *this;
    }

    Color operator*(double scalar) const {
        return Color(red * scalar, green * scalar, blue * scalar);
    }

    Color operator^(const Color& other) const {
        return Color(red * other.red, green * other.green, blue * other.blue);
    }

    Color& operator+=(const Color& other) {
        red += other.red;
        green += other.green;
        blue += other.blue;
        return *this;
    }
};

class PointLight{
public:
    Point light_pos;
    Color color;

    void draw()
    {
        glBegin(GL_POINTS);
        glColor3f(color.r(), color.g(), color.b());
        glVertex3f(light_pos.x(), light_pos.y(), light_pos.z());
        glEnd();
    }
};

class SpotLight{
public:
    PointLight point_light;
    Point light_dir;
    double cutoff_angle;

    void draw()
    {
        glBegin(GL_POINTS);
        glColor3f(point_light.color.r(), point_light.color.g(), point_light.color.b());
        glVertex3f(point_light.light_pos.x(), point_light.light_pos.y(), point_light.light_pos.z());
        glEnd();
    }
};

struct Ray{
public:
    Point start;
    Point dir;

    Ray(Point start, Point dir){
        this->start = start;
        dir.norm();
        this->dir = dir;
    }
};

class Object;
extern vector <Object*> objects;
extern vector <PointLight*> pointLights;
extern vector <SpotLight*> spotLights;
extern int recursion_level;

class Object {
public:
		Point reference_point;
		double height, width, length;
		Color color;
		double coefficients[4];
		int shine;

		Object(){
		}

		void setColor(Color color){
            this->color = color;
        }

		void setShine(int shine){
            this->shine = shine;
        }


		void setCoefficients(double coefficients[4]){
		    for(int i = 0; i < 4; i++){
                this->coefficients[i] = coefficients[i];
		    }
        }

        virtual Color getColorAt(Point point){
            return Color(this->color.r(), this->color.g(), this->color.b());
        }

        void set_color_for_diffuse_specular_reflection
        (Point intersectionPoint, Color intersectionPointColor, Ray r, PointLight pl, Point light_pos, Point light_dir, Color& color, double beta_rad, char type)
        {
            Ray ray_1 = Ray(light_pos, light_dir);
            Ray ray_n = getNormal(intersectionPoint, ray_1);
            Ray ray_r = Ray(intersectionPoint, ray_1.dir.sum( ray_n.dir.scale(2*ray_1.dir.dot(ray_n.dir)).scale(-1) ));

            double t1 = intersectionPoint.sum(light_pos.scale(-1)).length();
            if(t1 < epsilon) return;

            bool is_ray_1_obscured = false;

            for(int i = 0; i < objects.size(); i++){
                double t2 = objects[i]->intersect(ray_1, color, 0);
                is_ray_1_obscured = (t2 > 0) && (t1 - t2 > epsilon);
                if(is_ray_1_obscured) break;
            }

            if(!is_ray_1_obscured){
                double phongValue = max( 0.0, -r.dir.dot(ray_r.dir) );
                double lambertValue = max( 0.0, -ray_1.dir.dot(ray_n.dir) );

                if( type == 's' ){
                    pl.color = pl.color * pow( cos(beta_rad), epsilon2);
                }

                Color temp(intersectionPointColor * coefficients[1] * lambertValue);
                Color temp2(intersectionPointColor * coefficients[2] * pow(phongValue, shine));

                color += pl.color ^ temp;
                color += pl.color ^ temp2;
            }
        }

        void set_reflected_color_from_ambient(Color intersectionPointColor, Color& color){
            color = intersectionPointColor * coefficients[0];
        }

        void set_reflected_color_from_pointlight_source(Point intersectionPoint, Color intersectionPointColor, Ray r, Color& color){
            for(int i = 0; i < pointLights.size(); i++){
                Point light_dir = intersectionPoint.sum( pointLights[i]->light_pos.scale(-1) );
                light_dir.norm();
                set_color_for_diffuse_specular_reflection(intersectionPoint, intersectionPointColor, r, *(pointLights[i]), pointLights[i]->light_pos, light_dir, color, 0, 'p');
            }
        }

        void set_reflected_color_from_spotlight_source(Point intersectionPoint, Color intersectionPointColor, Ray r, Color& color){
            for(int i = 0; i < spotLights.size(); i++){
                Point light_dir = intersectionPoint.sum( spotLights[i]->point_light.light_pos.scale(-1) );
                light_dir.norm();

                double beta_rad = acos( light_dir.dot(spotLights[i]->light_dir) / (light_dir.length() * spotLights[i]->light_dir.length()));

                if(fabs(beta_rad * 180.0 / pi) < spotLights[i]->cutoff_angle){
                    set_color_for_diffuse_specular_reflection(intersectionPoint, intersectionPointColor, r, spotLights[i]->point_light, spotLights[i]->point_light.light_pos, light_dir, color, beta_rad, 's');
                }
            }
        }

        void set_color_for_recursive_reflection(Point intersectionPoint, Ray r, Color& color, int level){
            if(level < recursion_level){
                Ray ray_n = getNormal(intersectionPoint,r);
                Ray ray_r = Ray(intersectionPoint, r.dir.sum( ray_n.dir.scale(2*r.dir.dot(ray_n.dir)).scale(-1) ));
                ray_r.start = ray_r.start.sum( ray_r.dir.scale(epsilon) );

                int i_min = -1;
                double t = -1, t_min = INT_MAX;

                for(int i = 0; i < objects.size(); i++)
                {
                    t = objects[i]->intersect(ray_r, color, 0);
                    if(t > 0 && t < t_min){
                        t_min = t;
                        i_min = i;
                    }
                }

                if(i_min != -1)
                {
                    Color color_reflected(1, 1, 1);
                    double t = objects[i_min]->intersect(ray_r, color_reflected, level + 1);
                    color += color_reflected * coefficients[3];
                }
            }
        }

        virtual void draw(){
            return;
        }
		virtual double intersect(Ray r, Color &color, int level){
            return -1.0;
		};
        virtual Ray getNormal(Point point, Ray incidentRay){
            return Ray(Point(0, 0, 0), Point(1, 1, sqrt(2)));
        }
};

class General : public Object{
public:
    double A,B,C,D,E,F,G,H,I,J;

    General(){}

    Ray getNormal(Point point, Ray ray_in) override
    {
        return Ray(point, Point(2*A*point.x() + D*point.y() + E*point.z() + G,
               2*B*point.y() + D*point.x() + F*point.z() + H,
               2*C*point.z() + E*point.x() + F*point.y() + I));
    }

    double intersect(Ray r, Color &color, int level) override{

        Point R0 = r.start;
        Point Rd = r.dir;

        double a = A*Rd.x()*Rd.x() + B*Rd.y()*Rd.y() + C*Rd.z()*Rd.z() + D*Rd.x()*Rd.y() + E*Rd.x()*Rd.z() + F*Rd.y()*Rd.z();
        double b = 2*A*R0.x()*Rd.x() + 2*B*R0.y()*Rd.y() + 2*C*R0.z()*Rd.z() + D*(R0.x()*Rd.y() + Rd.x()*R0.y()) + E*(R0.x()*Rd.z() + Rd.x()*R0.z()) + F*(R0.y()*Rd.z() + Rd.y()*R0.z()) + G*Rd.x() + H*Rd.y() + I*Rd.z();
        double c = A*R0.x()*R0.x() + B*R0.y()*R0.y() + C*R0.z()*R0.z() + D*R0.x()*R0.y() + E*R0.x()*R0.z() + F*R0.y()*R0.z() + G*R0.x() + H*R0.y() + I*R0.z() + J;

        double D = pow(b, 2) - 4*a*c, t = -1;
        if(D < 0) t = -1;
        if(fabs(a) < epsilon) {
            t = -c/b;
        }
        else{
            double root_1 = (-b - sqrt(D))/(2*a);
            double root_2 = (-b + sqrt(D))/(2*a);

            if(root_1 < 0 && root_2 < 0) t = -1;

            if(root_2 < root_1){
                double t4 = root_1;
                root_1 = root_2;
                root_2 = t4;
            }

            if(root_1 > 0) {
                Point intersectionPoint = r.start.sum( r.dir.scale(root_1) );

                if((fabs(length) <= epsilon || (intersectionPoint.x() >= reference_point.x() && intersectionPoint.x() <= reference_point.x() + length)) &&
                    (fabs(width) <= epsilon || (intersectionPoint.y() >= reference_point.y() && intersectionPoint.y() <= reference_point.y() + width)) &&
                    (fabs(height) <= epsilon || (intersectionPoint.z() >= reference_point.z() && intersectionPoint.z() <= reference_point.z() + height)))
                {
                    t = root_1;
                }
            }

            if(root_2 > 0) {
                Point intersectionPoint = r.start.sum( r.dir.scale(root_2) );

                if((fabs(length) <= epsilon ||(intersectionPoint.x() >= reference_point.x() && intersectionPoint.x() <= reference_point.x() + length)) &&
                     (fabs(width) <= epsilon || (intersectionPoint.y() >= reference_point.y() && intersectionPoint.y() <= reference_point.y() + width)) &&
                     (fabs(height) <= epsilon || (intersectionPoint.z() >= reference_point.z() && intersectionPoint.z() <= reference_point.z() + height)))
                {
                    t = root_2;
                }
            }
            else{
                t = -1;
            }
        }

        if(t < 0) return -1;
        if(level == 0) return t;

        Point intersectionPoint = r.start.sum( r.dir.scale(t) );
        Color intersectionPointColor = getColorAt(intersectionPoint);

        set_reflected_color_from_ambient(intersectionPointColor, color);
        set_reflected_color_from_pointlight_source(intersectionPoint, intersectionPointColor, r, color);
        set_reflected_color_from_spotlight_source(intersectionPoint, intersectionPointColor, r, color);

        set_color_for_recursive_reflection(intersectionPoint, r, color, level);

        return t;
    }
};

class Triangle: public Object
{
public:
    Point a,b,c;

    Triangle(){
    }

    Triangle(Point a, Point b, Point c)
    {
        this->a = a;
        this->b = b;
        this->c = c;
    }

    Ray getNormal(Point point, Ray ray_in) override
    {
        Point triangle_normal = b.sum( a.scale(-1) ).cross(c.sum( a.scale(-1) ));
        triangle_normal.norm();

        return (ray_in.dir.dot(triangle_normal) < 0)? Ray(point, triangle_normal.scale(-1)) : Ray(point, triangle_normal);
    }

    void draw() override{
        glColor3f(color.r(), color.g(), color.b());
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x(), a.y(), a.z());
            glVertex3f(b.x(), b.y(), b.z());
            glVertex3f(c.x(), c.y(), c.z());
        }
        glEnd();
    }

    double intersect(Ray r, Color &color, int level) override{

        double Beta[3][3] = {
            {a.x() - r.start.x(), a.x() - c.x(), r.dir.x()},
            {a.y() - r.start.y(), a.y() - c.y(), r.dir.y()},
            {a.z() - r.start.z(), a.z() - c.z(), r.dir.z()}
        };
        double Gamma[3][3] = {
            {a.x() - b.x(), a.x() - r.start.x(), r.dir.x()},
            {a.y() - b.y(), a.y() - r.start.y(), r.dir.y()},
            {a.z() - b.z(), a.z() - r.start.z(), r.dir.z()}
        };
        double T[3][3] = {
            {a.x() - b.x(), a.x() - c.x(), a.x() - r.start.x()},
            {a.y() - b.y(), a.y() - c.y(), a.y() - r.start.y()},
            {a.z() - b.z(), a.z() - c.z(), a.z() - r.start.z()}
        };
        double A[3][3] {
            {a.x() - b.x(), a.x() - c.x(), r.dir.x()},
            {a.y() - b.y(), a.y() - c.y(), r.dir.y()},
            {a.z() - b.z(), a.z() - c.z(), r.dir.z()}
        };

        double Beta_det = ( Beta[0][0] * (Beta[1][1] * Beta[2][2] - Beta[1][2] * Beta[2][1]) ) -
                        ( Beta[0][1] * (Beta[1][0] * Beta[2][2] - Beta[1][2] * Beta[2][0]) ) +
                        ( Beta[0][2] * (Beta[1][0] * Beta[2][1] - Beta[1][1] * Beta[2][0]) );

        double Gamma_det = ( Gamma[0][0] * (Gamma[1][1] * Gamma[2][2] - Gamma[1][2] * Gamma[2][1]) ) -
                        ( Gamma[0][1] * (Gamma[1][0] * Gamma[2][2] - Gamma[1][2] * Gamma[2][0]) ) +
                        ( Gamma[0][2] * (Gamma[1][0] * Gamma[2][1] - Gamma[1][1] * Gamma[2][0]) );

        double T_Det = ( T[0][0] * (T[1][1] * T[2][2] - T[1][2] * T[2][1]) ) -
                        ( T[0][1] * (T[1][0] * T[2][2] - T[1][2] * T[2][0]) ) +
                        ( T[0][2] * (T[1][0] * T[2][1] - T[1][1] * T[2][0]) );

        double A_Det = ( A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) ) -
                        ( A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) ) +
                        ( A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]) );

        double beta = Beta_det / A_Det;
        double gamma = Gamma_det / A_Det;
        double t = T_Det / A_Det;

        t = (beta + gamma < 1 && beta > 0 && gamma > 0 && t > 0) ? t : -1;

        if(t < 0) return -1;
        if(level == 0) return t;

        Point intersectionPoint = r.start.sum( r.dir.scale(t) );
        Color intersectionPointColor = getColorAt(intersectionPoint);

        set_reflected_color_from_ambient(intersectionPointColor, color);
        set_reflected_color_from_pointlight_source(intersectionPoint, intersectionPointColor, r, color);
        set_reflected_color_from_spotlight_source(intersectionPoint, intersectionPointColor, r, color);

        set_color_for_recursive_reflection(intersectionPoint, r, color, level);

        return t;
    }
};

class Sphere : public Object
{
public:

    Sphere(){
    }

    Sphere(Point center, double radius){
        reference_point = center;
        length = radius;
    }

    Ray getNormal(Point point, Ray incidentRay) override{
        return Ray(point, point.sum(reference_point.scale(-1)));
    }

    void draw()override{
        int N = 30;
        vector<vector<Point>> spherePoints(N+1, vector<Point>(N+1));

        for (int i = 0; i <= N; i++)
        {
            float theta = (float)i / N * 2 * pi;
            for (int j = 0; j <= N; j++)
            {
                float phi = (float)j / N * pi;
                spherePoints.at(i).at(j) = {length * sin(phi) * cos(theta), length * sin(phi) * sin(theta), length * cos(phi)};
            }
        }

        glPushMatrix();
        glTranslatef(reference_point.x(), reference_point.y(), reference_point.z());
        glColor3f(color.r(), color.g(), color.b());
        for (int i = 0; i < N; i++)
        {
            glBegin(GL_QUAD_STRIP);
            for (int j = 0; j <= N; j++)
            {
                glVertex3f(spherePoints[i][j].x(), spherePoints[i][j].y(), spherePoints[i][j].z());
                glVertex3f(spherePoints[i + 1][j].x(), spherePoints[i + 1][j].y(), spherePoints[i + 1][j].z());
            }
            glEnd();
        }
        glPopMatrix();
    }

    double intersect(Ray r, Color &color, int level) override{

        r.start = r.start.sum( reference_point.scale(-1) );

        double a = 1.0;
        double b = 2 * r.dir.dot(r.start);
        double c = r.start.dot(r.start) - pow(length, 2);

        double D = pow(b, 2) - 4 * a * c;
        double t = -1;

        if(D < 0){
            t = -1;
        }
        else{
            double root_1 = (-b - sqrt(D)) / (2 * a);
            double root_2 = (-b + sqrt(D)) / (2 * a);

            if(root_2 < root_1){
                double t4 = root_1;
                root_1 = root_2;
                root_2 = t4;
            }

            t = (root_1 > 0) ? root_1 : ( (root_2 > 0) ? root_2 : -1);
        }

        if(t < 0) return -1;
        if(level == 0) return t;

        Point intersectionPoint = r.start.sum( r.dir.scale(t) );
        Color intersectionPointColor = getColorAt(intersectionPoint);

        set_reflected_color_from_ambient(intersectionPointColor, color);
        set_reflected_color_from_pointlight_source(intersectionPoint, intersectionPointColor, r, color);
        set_reflected_color_from_spotlight_source(intersectionPoint, intersectionPointColor, r, color);

        set_color_for_recursive_reflection(intersectionPoint, r, color, level);

        return t;
    }
};

class Floor : public Object
{
public:

    int N;

    Floor(){
        N = 1;
    }

    Floor(int floorWidth, int tileWidth){
        N = floorWidth / tileWidth;
        reference_point = Point(-floorWidth / 2, -floorWidth / 2, 0);
        length = tileWidth;
    }

    Color getColorAt(Point point) override{
        int i = (point.x() - reference_point.x()) / length;
		int j = (point.y() - reference_point.y()) / length;

        return (i < 0 || i >= N || j < 0 || j >= N) ? Color(0, 0, 0) : ( (((i + j) % 2) != 0) ? Color(0, 0, 0): Color(1, 1, 1) );
    }

    Ray getNormal(Point point, Ray ray_in) override{
        return (ray_in.dir.z() > 0) ? Ray(point, Point(0, 0, 1)) : Ray(point, Point(0, 0, -1));
    }

    void draw() override{
        glPushMatrix();
        glTranslatef(reference_point.x(), reference_point.y(), reference_point.z());

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if ((i + j) % 2 == 0) {
                    glColor3f(0, 0, 0);
                } else {
                    glColor3f(1, 1, 1);
                }

                glBegin(GL_QUADS);
                glVertex2f(i * length, j * length);
                glVertex2f((i + 1) * length, j * length);
                glVertex2f((i + 1) * length, (j + 1) * length);
                glVertex2f(i * length, (j + 1) * length);
                glEnd();
            }
        }

        glPopMatrix();
    }

    double intersect(Ray r, Color &color, int level) override{
        Point floor_normal = Point(0, 0, 1);
        double t = -1;

        if (round(floor_normal.dot(r.dir) * 100) == 0){
			t = -1;
        }

        t = - floor_normal.dot(r.start) / floor_normal.dot(r.dir);

        Point intersectionPoint = r.start.sum( r.dir.scale(t) );

        if(intersectionPoint.x() <= reference_point.x() ||
           intersectionPoint.x() >= -reference_point.x() ||
           intersectionPoint.y() <= reference_point.y() ||
           intersectionPoint.y() >= -reference_point.y())
        {
            t = -1;
        }

        if(t < 0) return -1;
        if(level == 0) return t;

        Color intersectionPointColor = getColorAt(intersectionPoint);

        set_reflected_color_from_ambient(intersectionPointColor, color);
        set_reflected_color_from_pointlight_source(intersectionPoint, intersectionPointColor, r, color);
        set_reflected_color_from_spotlight_source(intersectionPoint, intersectionPointColor, r, color);

        set_color_for_recursive_reflection(intersectionPoint, r, color, level);

        return t;
    }
};

#endif

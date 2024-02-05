#include <vector>
#include <stack>
#include <string>
#include <algorithm>
#include <iomanip>
#include "bitmap_image.hpp"
#define PI 2.0*acos(0.0)

using namespace std;

static unsigned long long int g_seed = 1;
inline int _random(){
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>16)&0x7FFF;
}

class Color{
private:
    int red, green, blue;

public:
    Color(){
        red = green = blue = 0;
    }

    int r(){ return this->red; }
    int g(){ return this->green; }
    int b(){ return this->blue; }

    void setRGB(){
        red = _random();
        green = _random();
        blue = _random();
    }
};

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

    double dot(Point &other){
        return (this->x1 * other.x() + this->y1 * other.y() + this->z1 * other.z());
    }

    Point cross(Point &other){
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

    Point rotate_rodrigues(Point axis, double theta_rad){
        axis.norm();
        Point res1 = this->scale( cos(theta_rad) );
        Point res2 = axis.scale( axis.dot(*this) * (1-cos(theta_rad)) );
        Point res3 = axis.cross(*this).scale( sin(theta_rad) );
        return res1.sum(res2).sum(res3);
    }
};

class Triangle{
private:
    vector<Point> vertices;
    Color color;

public:
    Triangle(Point a, Point b, Point c){
        vertices = {a, b, c};
    }

    Triangle(){
        vertices = vector<Point>(3, Point());
    }

    Point& a(){ return this->vertices[0]; }

    Point& b(){ return this->vertices[1]; }

    Point& c(){ return this->vertices[2]; }

    void setVertices(vector<Point> verts){
        vertices.assign(verts.begin(), verts.end());
    }

    Color& getColor(){ return this->color; }

    void _swap(int i, int j){
        Point temp = vertices[i];
        vertices[i] = vertices[j];
        vertices[j] = temp;
    }
};

class Matrix{
private:
    vector<vector<double>> matrix4x4;

public:
    Matrix(){
        matrix4x4 = vector<vector<double>>(4, vector<double>(4, 0.0));
    }

    Matrix(vector<vector<double>> other){
        matrix4x4 = vector<vector<double>>(other.begin(), other.end());
    }

    vector<vector<double>>& getMatrix(){ return matrix4x4; }

    void setMatrix(vector<vector<double>> other){
        matrix4x4.assign(other.begin(), other.end());
    }

    Matrix multiplyMatrix(Matrix other){
        Matrix ans;
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                for(int k=0; k<4; k++){
                    ans.getMatrix().at(i).at(j) += matrix4x4[i][k] * other.getMatrix()[k][j];
                }
            }
        }
        return ans;
    }

    Point multiplyPoint(Point other){
        vector<double> vecPoint(4, 0.0);
        vector<double> point = other.getCoords();

        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                vecPoint.at(i) += matrix4x4[i][j] * point[j];
            }
        }
        Point resultPoint(vecPoint[0], vecPoint[1], vecPoint[2],vecPoint[3]);
        resultPoint.scale_w_to_one();
        return resultPoint;
    }

    Triangle multiplyTriangle(Triangle t){
        return Triangle( multiplyPoint(t.a()),
                        multiplyPoint(t.b()),
                        multiplyPoint(t.c()) );
    }

    void _translate(Point translate_amount){
        setMatrix({
            {1, 0, 0, translate_amount.x()},
            {0, 1, 0, translate_amount.y()},
            {0, 0, 1, translate_amount.z()},
            {0, 0, 0, 1}
        });
    }

    void _scale(Point scaleFactor){
        setMatrix({
            {scaleFactor.x(), 0, 0, 0},
            {0, scaleFactor.y(), 0, 0},
            {0, 0, scaleFactor.z(), 0},
            {0, 0, 0, 1}
        });
    }

    void _rotate(Point axis,double theta_rad){
        axis.norm();
        setMatrix({
            {1,0,0,0},
            {0,1,0,0},
            {0,0,1,0},
            {0,0,0,1}
        });

        Point i1(1.0, 0.0, 0.0);
        Point j1(0.0, 1.0, 0.0);
        Point k1(0.0, 0.0, 1.0);

        Point r1 = i1.rotate_rodrigues(axis, theta_rad);
        Point r2 = j1.rotate_rodrigues(axis, theta_rad);
        Point r3 = k1.rotate_rodrigues(axis, theta_rad);

        setMatrix({
            {r1.x(), r2.x(), r3.x(), 0},
            {r1.y(), r2.y(), r3.y(), 0},
            {r1.z(), r2.z(), r3.z(), 0},
            {0, 0, 0, 1}
        });
    }

    void printMatrix(string matName){
        cout << "\nPrinting the matrix \'" << matName << "\': \n";
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                cout << matrix4x4[i][j] << " ";
            }
            cout << endl;
        }
    }
};

struct camera{
    Point eye;
    Point look;
    Point up;
    Point r_hat;
    Point l_hat;
    Point u_hat;
} Camera;

struct perspective {
    double fovX;
    double fovY;
    double aspect;
    double near;
    double far;
    double t;
    double r;
} Pers;

class ZBuffer{
    ifstream fin_3, fin_c;
    ofstream fout;
    int Screen_Width, Screen_Height;
    double X_Min, X_Max, Y_Min, Y_Max, Z_Min, Z_Max;
    double dx, dy;
    double Top_Y, Left_X;
    vector<Triangle> triangles;
    vector<vector<double>> z_buffer;
    bitmap_image image;
public:
    ZBuffer(){
        X_Min = Y_Min = Z_Min = -1;
        X_Max = Y_Max = Z_Max = 1;
        Screen_Width = Screen_Height = 500;
        dx = (X_Max - X_Min) / Screen_Width;
        dy = (Y_Max - Y_Min) / Screen_Height;
        Top_Y = Y_Max - dy/2.0;
        Left_X = X_Min + dx/2.0;
    }

    void readData(){
        fin_c.open("config.txt");
        fin_c >> Screen_Width >> Screen_Height;

        fin_3.open("stage3.txt");
        Triangle t;
        while( fin_3 >> t.a().x() >> t.a().y() >> t.a().z() >> t.b().x() >> t.b().y() >> t.b().z() >> t.c().x() >> t.c().y() >> t.c().z() )
        {
            if( t.b().y() < t.a().y() ) t._swap(0,1);
            if( t.c().y() < t.a().y() ) t._swap(0,2);
            if( t.c().y() < t.b().y() ) t._swap(1,2);
            t.getColor().setRGB();
            triangles.push_back(t);
        }
    }

    void init_z_buffer(){
        z_buffer = vector<vector<double>>(Screen_Height, vector<double>(Screen_Width, Z_Max));
        for(int i=0; i<Screen_Height; i++){
            for(int j=0; j<Screen_Width; j++){
                z_buffer[i][j] = Z_Max;
            }
        }
    }

    void init_frame_buffer(){
        image.setwidth_height(Screen_Width, Screen_Height);
        image.set_all_channels(0, 0, 0);
    }

    vector<double> find_intersect_a_b(Triangle t, double y){
        double x1 = 0, z1 = 0;
        double x2 = -1, z2 = -1;

        if( t.b().y() != t.a().y() && t.a().y() != t.c().y() ){
            x1 = t.a().x() + (t.b().x() - t.a().x()) * (y - t.a().y()) / (t.b().y() - t.a().y());
            x2 = t.a().x() + (t.c().x() - t.a().x()) * (y - t.a().y()) / (t.c().y() - t.a().y());

            z1 = t.a().z() + (t.b().z() - t.a().z()) * (y - t.a().y()) / (t.b().y() - t.a().y());
            z2 = t.a().z() + (t.c().z() - t.a().z()) * (y - t.a().y()) / (t.c().y() - t.a().y());

            if( x1 > x2 ) {
                swap(x1, x2);
                swap(z1, z2);
            }
        }

        return {x1, x2, z1, z2};
    }

    vector<double> find_intersect_b_c(Triangle t, double y){
        double x1 = 0, z1 = 0;
        double x2 = -1, z2 = -1;

        if( t.b().y() != t.c().y() && t.a().y() != t.c().y() ){
            x1 = t.c().x() + (t.b().x() - t.c().x()) * (y - t.c().y()) / (t.b().y() - t.c().y());
            x2 = t.c().x() + (t.c().x() - t.a().x()) * (y - t.c().y()) / (t.c().y() - t.a().y());

            z1 = t.c().z() + (t.b().z() - t.c().z()) * (y - t.c().y()) / (t.b().y() - t.c().y());
            z2 = t.c().z() + (t.c().z() - t.a().z()) * (y - t.c().y()) / (t.c().y() - t.a().y());

            if( x1 > x2 ) {
                swap(x1, x2);
                swap(z1, z2);
            }
        }

        return {x1, x2, z1, z2};
    }

    void calc_z_buffer_update_pixels(Triangle t, double y, double x1, double x2, double z1, double z2){
        for(double x = x1; x <= x2; x += dx){
            if( x2 == x1 ) continue;

            int row = (Top_Y - y) / dy;
            int col = (x - Left_X) / dx;

            double z = z1 + (z2 - z1) * (x - x1) / (x2 - x1);
            if( z < z_buffer[row][col] && z > Z_Min ){
                z_buffer[row][col] = z;
                image.set_pixel(col, row, t.getColor().r(), t.getColor().g() , t.getColor().b());
            }
        }
    }

    void apply_scan_conversion_a_b(Triangle t, double y1, double y2){
        for(double y = y1; y <= y2; y += dy){
            vector<double> result = find_intersect_a_b(t, y);
            double x1 = result[0];
            double x2 = result[1];
            double z1 = result[2];
            double z2 = result[3];

            x1 = max(x1, X_Min);
            x2 = min(x2, X_Max);

            calc_z_buffer_update_pixels(t, y, x1, x2, z1, z2);
        }
    }

    void apply_scan_conversion_b_c(Triangle t, double y1, double y2){
        for(double y = y1; y <= y2; y += dy){
            vector<double> result = find_intersect_b_c(t, y);
            double x1 = result[0];
            double x2 = result[1];
            double z1 = result[2];
            double z2 = result[3];

            x1 = max(x1, X_Min);
            x2 = min(x2, X_Max);

            calc_z_buffer_update_pixels(t, y, x1, x2, z1, z2);
        }
    }

    void apply_clipping_scan_conversion_algo(){
        for (Triangle t : triangles){
            double y1 = max(t.a().y(), Y_Min);
            double y2 = min(t.b().y(), Y_Max);

            apply_scan_conversion_a_b(t, y1, y2);

            y1 = max(t.b().y(), Y_Min);
            y2 = min(t.c().y(), Y_Max);

            apply_scan_conversion_b_c(t, y1, y2);
        }
    }

    void save_image_z_buffer_val(){
        image.save_image("out.bmp");

        for(int i = 0; i < Screen_Height; i++){
            for(int j = 0; j < Screen_Width; j++){
                if (z_buffer[i][j] < Z_Max) {
                    fout << setprecision(6) << fixed << z_buffer[i][j] << "\t";
                }
            }
            fout<<endl;
        }
    }

    void free_mem_image_z_buffer(){
        fin_3.close();
        fin_c.close();
        fout.close();
        z_buffer.clear();
    }
};

int main(){
    ifstream fin("..\\Offline2\\IOs\\5\\scene.txt");
    ofstream fout("stage1.txt");

    fin >> Camera.eye.x() >> Camera.eye.y() >> Camera.eye.z();
    fin >> Camera.look.x() >> Camera.look.y() >> Camera.look.z();
    fin >> Camera.up.x() >> Camera.up.y() >> Camera.up.z();
    fin >> Pers.fovY >> Pers.aspect >> Pers.near >> Pers.far;

    /////////////////////////////////////////// stage-1 ///////////////////////////////////////////
    stack<Matrix> S;
    S.push(Matrix({
      {1, 0, 0, 0},
      {0, 1, 0, 0},
      {0, 0, 1, 0},
      {0, 0, 0, 1}
    }));

    string cmd;
    while(fin >> cmd){
        if( cmd == "triangle" ){
            Triangle t;
            fin >> t.a().x() >> t.a().y() >> t.a().z() >> t.b().x() >> t.b().y() >> t.b().z() >> t.c().x() >> t.c().y() >> t.c().z();

            t = S.top().multiplyTriangle(t);
            fout << fixed << setprecision(7) << t.a().x() << " " << t.a().y() << " " << t.a().z() << endl;
            fout << fixed << setprecision(7) << t.b().x() << " " << t.b().y() << " " << t.b().z() << endl;
            fout << fixed << setprecision(7) << t.c().x() << " " << t.c().y() << " " << t.c().z() << endl << endl;
        }
        else if( cmd == "translate" ){
            Point point;
            fin >> point.x() >> point.y() >> point.z();

            Matrix mat;
            mat._translate(point);
            S.top() = S.top().multiplyMatrix(mat);
        }
        else if( cmd == "scale" ){
            Point point;
            fin >> point.x() >> point.y() >> point.z();

            Matrix mat;
            mat._scale(point);
            S.top() = S.top().multiplyMatrix(mat);
        }
        else if( cmd == "rotate" ){
            double theta_deg;
            Point point;
            fin >> theta_deg >> point.x() >> point.y() >> point.z();

            Matrix mat;
            mat._rotate(point, theta_deg * PI / 180.0);
            S.top() = S.top().multiplyMatrix(mat);
        }
        else if( cmd == "push" ){
            S.push(S.top());
        }
        else if( cmd == "pop" ){
            S.pop();
        }
        else if( cmd == "end" ){
            break;
        }
    }

    fin.close();
    fout.close();


    /////////////////////////////////////////// stage-2 ///////////////////////////////////////////
    Camera.l_hat = Camera.look.sum( Camera.eye.scale(-1) );
    Camera.l_hat.norm();
    Camera.r_hat = Camera.l_hat.cross(Camera.up);
    Camera.r_hat.norm();
    Camera.u_hat = Camera.r_hat.cross(Camera.l_hat);
    Camera.u_hat.norm();

    Matrix T;
    T._translate( Camera.eye.scale(-1) );
    Matrix R({
        {Camera.r_hat.x(), Camera.r_hat.y(), Camera.r_hat.z(), 0},
        {Camera.u_hat.x(), Camera.u_hat.y(), Camera.u_hat.z(), 0},
        {-Camera.l_hat.x(), -Camera.l_hat.y(), -Camera.l_hat.z(), 0},
        {0, 0, 0, 1}
    });
    Matrix V = R.multiplyMatrix(T);

    fin.open("stage1.txt");
    fout.open("stage2.txt");
    Triangle t;
    while( fin >> t.a().x() >> t.a().y() >> t.a().z() >> t.b().x() >> t.b().y() >> t.b().z() >> t.c().x() >> t.c().y() >> t.c().z() )
    {
        t = V.multiplyTriangle(t);

        fout << fixed << setprecision(7) << t.a().x() << " " << t.a().y() << " " << t.a().z() << endl;
        fout << fixed << setprecision(7) << t.b().x() << " " << t.b().y() << " " << t.b().z() << endl;
        fout << fixed << setprecision(7) << t.c().x() << " " << t.c().y() << " " << t.c().z() << endl << endl;
    }

    fin.close();
    fout.close();


    /////////////////////////////////////////// stage-3 ///////////////////////////////////////////
    Pers.fovX = Pers.fovY * Pers.aspect;
    Pers.t = Pers.near * tan(Pers.fovY * PI / 360.0);
    Pers.r = Pers.near * tan(Pers.fovX * PI / 360.0);

    Matrix P({
        {Pers.near / Pers.r, 0, 0, 0},
        {0, Pers.near / Pers.t, 0, 0},
        {0, 0, -(Pers.far + Pers.near) / (Pers.far - Pers.near), -(2.0 * Pers.far * Pers.near) / (Pers.far - Pers.near)},
        {0, 0, -1, 0}
    });

    fin.open("stage2.txt");
    fout.open("stage3.txt");
    while( fin >> t.a().x() >> t.a().y() >> t.a().z() >> t.b().x() >> t.b().y() >> t.b().z() >> t.c().x() >> t.c().y() >> t.c().z() )
    {
        t = P.multiplyTriangle(t);

        fout << fixed << setprecision(7) << t.a().x() << " " << t.a().y() << " " << t.a().z() << endl;
        fout << fixed << setprecision(7) << t.b().x() << " " << t.b().y() << " " << t.b().z() << endl;
        fout << fixed << setprecision(7) << t.c().x() << " " << t.c().y() << " " << t.c().z() << endl << endl;
    }

    fin.close();
    fout.close();


    /////////////////////////////////////////// stage-4 ///////////////////////////////////////////
    ZBuffer obj;
    obj.readData();
    obj.init_z_buffer();
    obj.init_frame_buffer();
    obj.apply_clipping_scan_conversion_algo();
    obj.save_image_z_buffer_val();
    obj.free_mem_image_z_buffer();

    return 0;
}

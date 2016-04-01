
//
// raytrace.cpp
//


#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <float.h>              // for FLT_MAX
#include <math.h>


using namespace std;


// Structs and enums //////////////////////////////////////
struct Sphere;
struct Light_source;
struct Ray
{
    vec4 origin;
    vec4 dir;
};

struct Sphere {
    string name;
    vec4 pos;                   // position: vec4, last row is 1
    vec3 scale;                 // scale: scl_x, scl_y, scl_z
    vec4 color;                 // color: r, g, b, alpha
    vec4 coeff;                 // coefficients: Ka, Kd, Ks, Kr
    float n;                    // exponent
};

struct Light_source {
    string name;                // light source name
    vec4 pos;                   // point light source pos
    vec4 color;                 // light source color
};

enum INTERSECT_TYPES
{
    INTERSECT_1 = 0,
    INTERSECT_2 = 1,
    INTERSECT_3 = 2,
    INTERSECT_4 = 3,
    INTERSECT_5 = 4,
    INTERSECT_NONE = 5
};


// Global Variables //////////////////////////////////////////

vector<vec4> g_colors;
vector<Light_source*> lights;       // light sources(<=5)
vector<Sphere*> spheres;            // spheres(<=5)
vector<mat4> Msphere;               // original matrices for spheres
vector<mat4> Msphere_inversed;      // inversed matrices after scaling

float g_left;                       //--------------
float g_right;                      //              |
float g_top;                        //              |
float g_bottom;                     //              |----> global vars to store graphics info
float g_near;                       //              |
int g_width;                        //              |
int g_height;                       //--------------
vec4 Background;                    // background color
vec4 Ambient;                       // ambient color
string output;                      // output file name

// helper functions ///////////////////////////////////////////
inline int n_spheres(void) {
    return (int)spheres.size();
}

inline int n_lights(void) {
    return (int)lights.size();
}

// given a sphere and a point on that sphere, return normal n' = M^(-t)n
vec4 normal(const vec4& p, int i)
{
    vec4 normal = p - spheres[i]->pos;
    mat4 temp = transpose(Msphere_inversed[i]);
    normal = temp * Msphere_inversed[i] * normal;
    normal.w = 0;
    return normal;
}

// clamp color that has component larger than 1
void clamp(vec4& color)
{
    color.w = 1.0f;
    if (color.x > 1.0f)
        color.x = 1.0f;
    if (color.y > 1.0f)
        color.y = 1.0f;
    if (color.z > 1.0f)
        color.z = 1.0f;
}

// calculate the inverse Trans matrix of spheres, should only call once
void calculate_inverse(void) {
    int nSphere = n_spheres();
    // for each trans matrix, calculate the inverse matrix and store
    for (int i = 0; i < nSphere; i++) {
        mat4 temp;
        InvertMatrix(Msphere[i], temp);
        Msphere_inversed.push_back(temp);
    }
}

// ////////////////////////////////////////////////////////////


// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3, const string& s4 = "1")
{
    stringstream ss(s1 + " " + s2 + " " + s3 + " " + s4);
    vec4 result;
    ss >> result.x >> result.y >> result.z >> result.w;
    return result;
}

vec3 toVec3(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec3 result;
    ss >> result.x >> result.y >> result.z;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    if (vs[0] == "\r" || vs[0][0] == '\0') {    // skip new line
        return;
    } else if (vs[0] == "NEAR") {
        g_near = toFloat(vs[1]);
    } else if (vs[0] == "LEFT") {
        g_left = toFloat(vs[1]);
    } else if (vs[0] == "RIGHT") {
        g_right = toFloat(vs[1]);
    } else if (vs[0] == "BOTTOM") {
        g_bottom = toFloat(vs[1]);
    } else if (vs[0] == "TOP") {
        g_top = toFloat(vs[1]);
    } else if (vs[0] == "RES") {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
        
    } else if (vs[0] == "SPHERE") {
        Sphere* temp = new Sphere;
        temp->name  = vs[1];
        temp->pos = toVec4(vs[2], vs[3], vs[4]);
        temp->scale = toVec3(vs[5], vs[6], vs[7]);
        temp->color = toVec3(vs[8], vs[9], vs[10]);
        temp->coeff = toVec4(vs[11], vs[12], vs[13], vs[14]);
        temp->n = toFloat(vs[15]);
        spheres.push_back(temp);
        // also store the transformation matrix
        mat4 Mtemp = Translate(temp->pos) * Scale(temp->scale);
        Msphere.push_back(Mtemp);
        
    } else if (vs[0] == "LIGHT") {
        Light_source* temp = new Light_source;
        temp->name = vs[1];
        temp->pos = toVec4(vs[2], vs[3], vs[4]);
        temp->color = toVec4(vs[5], vs[6], vs[7]);
        lights.push_back(temp);
        
    } else if (vs[0] == "BACK") {
        Background = toVec4(vs[1], vs[2], vs[3]);
    } else if (vs[0] == "AMBIENT") {
        Ambient = toVec4(vs[1], vs[2], vs[3]);
    } else if (vs[0] == "OUTPUT") {
        output = vs[1];
    } else {
        cout << "error!" << endl;
        cout << int(vs[0][0]) << endl;
        exit(1);
    }
    
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
    
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine


/*  
 intersect() should takes a ray as an input, caculates the intersect
 point, update type to one of the INTERSECT_TYPES. The caculation is based on
 the coordinates of spheres(to find closest intersection point),
 OR light sources (for shadow rays).
*/
vec4 intersect(const Ray& ray, int& type, const float t_min, const float t_max, bool& IF_HOLLOW) {

    // need for closest intersection point
    float t_cur_min = 999;
    IF_HOLLOW = false;
    
    int ret_type = INTERSECT_NONE;
    // ray = S + Ct, and S is just the origin
    vec4 S = ray.origin;
    vec4 C = ray.dir;
    
    int nSphere = n_spheres();
    
    for (int i = 0; i < nSphere; i++) {
        
        vec4 S_t = Msphere_inversed[i] * S;
        vec4 C_t = Msphere_inversed[i] * C;

        // C^2*t^2 + 2(S*tC) + S^2 - 1 = 0
        float a = dot(C_t, C_t);
        float b = dot(S_t, C_t);
        float c = dot(S_t, S_t) - 1 - 1;    // extra 1 for w coordinate
        float det = b*b - a*c;
    
        if (det < 0)                        // no solution, the ray doesn't intersect the sphere
            continue;
        
        // NOTE t1 > t2
        float t1 = (-b + sqrt(det))/a;
        float t2 = (-b - sqrt(det))/a;

        if (t_min == 1.0f) {
            if (0 < t2 && t2 < 1 && t1 > 1)
                IF_HOLLOW = true;
        }
        
        if (det == 0) { // t1 == t2
            if (t1 <= t_min || t1 >= t_max)   // hollow sphere
                continue;
            else if (t1 < t_cur_min) {      // check if this intersection is the closest
                t_cur_min = t1;
                ret_type = i;
            }
            else continue;
        } else {
            // t1 != t2, and t2 < t1
            if (t2 > t_min && t2 < t_max) {
                if (t2 < t_cur_min) {
                    t_cur_min = t2;
                    ret_type = i;
                }
                // if the smaller valid solution t2 is not the smallest so far,
                // there's no need to check the larger t1.
                continue;
            }
            // if get here, t2 must be invalid, simply check t1 and update if
            // it is the smallest t so far.
            else if (t1 > t_min && t1 < t_max) {
                if (t1 < t_cur_min) {
                    t_cur_min = t1;
                    ret_type = i;
                }
            }
            else continue;
        }
        
    }
    
    // use the UNTRANSFORMED equation to calculate intersection
    // if no intersection, we don't care about return value
    vec4 closest_intersect_point = S + t_cur_min*C;
    
    type = ret_type;
#ifdef DEBUG1
    cout << "type: " << ret_type << endl;
#endif
    
    return closest_intersect_point;
}

/*  
 shadowRay() should takes the the closest intersect point as an input,
 and calculate if lights to that point is blocked or not. If not, it gives
 a color for the pixel (possibly a mixture of multiple lights).
*/
vec4 shadowRay(const vec4& p_intersect, int i_intersect, vec4 normal) {
    
    // time range
    const float t_max = 1.000f;
    const float t_min = 0.0001f;
    
    vec4 N = normal;
    
    vec4 Oc = spheres[i_intersect]->color;              // object color
    float Ks = spheres[i_intersect]->coeff.z;           // specular reflectance coefficient
    float Kd = spheres[i_intersect]->coeff.y;           // diffuse reflectance coefficient
    vec4 ret_color = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    Ray ray;
    ray.origin = p_intersect;
    
    vec4 diffuse, specular;
    
    int nlights = n_lights();
    for (int i = 0; i < nlights; i++)
    {
        ray.dir = lights[i]->pos - p_intersect;
        int TYPE = 0;
        bool IF_HOLLOW;
        intersect(ray, TYPE, t_min, t_max, IF_HOLLOW);
        
        ray.dir = normalize(ray.dir);
    
        N = normalize(N);
                
        if (TYPE == INTERSECT_NONE) {
            vec4 Ip = lights[i]->color;
            vec4 L = ray.dir;           // direction to lights[i]
            vec4 t_diffuse, t_specular;
            if (dot(N, L) < 0) {
                continue;
            }
            else {
                t_diffuse = Kd * Ip * dot(N, L) * Oc;          // diffuse color
                diffuse += t_diffuse;
            
                vec4 R = 2 * N * dot(N, L) - L;
                vec4 V = normalize(vec4(0.0f, 0.0f, 0.0f, 1.0f) - p_intersect);
                t_specular = Ks * Ip * pow(dot(R, V), spheres[i_intersect]->n);

                specular += t_specular;
            }
        }
    }
    
    ret_color = diffuse + specular;

    return ret_color;
}

// -------------------------------------------------------------------
// recursive ray tracing
vec4 trace(const Ray& ray, int recur_depth)
{
    // depth should be no more than 3
    if (recur_depth == 4)
        return Background;
    
    int INDEX;                          // which ball the ray is intersecting
    bool IF_HOLLOW;
    vec4 local_color, ambient_color, reflect_color, pixel_color;
    vec4 intersect_point;
    
    
    if (recur_depth == 1)   // first hitpoint:  t > 1
        intersect_point = intersect(ray, INDEX, 1.0f, FLT_MAX, IF_HOLLOW);
    else                    // else(reflection): t > 0
        intersect_point = intersect(ray, INDEX, 0.0001f, FLT_MAX, IF_HOLLOW);
    
    if (INDEX == INTERSECT_NONE)
        return Background;
    
    vec4 N = normal(intersect_point, INDEX);
    
    if(IF_HOLLOW == true) {
        N.x = -N.x; N.y = -N.y; N.z = -N.z; N.w = 0.0f;
    }
    
    //////////////// Ambient Color /////////////////
    float Ka = spheres[INDEX]->coeff.x; // Ka for Ambient Surface Reflectance Coefficient
    vec4 Oa = spheres[INDEX]->color;    // Oa for Object color
    vec4 Ia = Ambient;                  // Ia for Ambient color Intensity
    ambient_color = Ka * Ia * Oa;
    pixel_color += ambient_color;
    
    
    //////////////// Local Color(diffuse + specular) ////////////////
    local_color = shadowRay(intersect_point, INDEX, N);
    pixel_color += local_color;
    

    /////////////////// Reflect Color /////////////////////
    Ray reflected_ray;
    reflected_ray.origin = intersect_point;
    N = normalize(N);
    vec4 L = normalize(ray.dir);
    L.x = -L.x; L.y = -L.y; L.z = -L.z; L.w = 0;
    reflected_ray.dir = 2 * N * dot(N, L) - L;
    float Kr = spheres[INDEX]->coeff.w;
    
    reflect_color = trace(reflected_ray, recur_depth+1);
    
    // weird, vec comparison doesn't seem to work here
    if( reflect_color.x != Background.x ||
        reflect_color.y != Background.y ||
        reflect_color.z != Background.z ||
        reflect_color.w != Background.w )
        
        pixel_color += Kr * reflect_color;
    
    clamp(pixel_color);

    return pixel_color;
}

/*
 return direction to (ix, iy), UNNORMALIZED
 */
vec4 getDir(int ix, int iy)
{
    float x_len, y_len, x_ratio, y_ratio;
    x_len = g_right - g_left;
    y_len = g_top - g_bottom;
    x_ratio = (float)ix / g_width;
    y_ratio = (float)iy / g_height;
    
    float x = g_left + x_ratio * x_len;         // linear interpolate x coordinate
    float y = g_bottom + y_ratio * y_len;       // linear interpolate y coordinate
    
    vec4 dir = vec4(x, y, -g_near, 0.0f);
    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);  // eyes located at origin
    ray.dir = getDir(ix, iy);                   // shadow ray from eye to pixel
    
    vec4 color = trace(ray, 1);
 
    setColor(ix, iy, color);
}

void render()
{
#ifndef DEBUG1
     for (int iy = 0; iy < g_height; iy++)
         for (int ix = 0; ix < g_width; ix++)
             renderPixel(ix, iy);
#endif
#ifdef DEBUG1
    renderPixel(260, 321);
#endif
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, const char* fname, unsigned char* pixels)
{
    FILE *fp;
    const int maxVal=255;
    
    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);
    
    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }
    
    fclose(fp);
}

void saveFile()
{
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
    
    char *filename = new char[output.length() + 1];
    strcpy(filename, output.c_str());
    savePPM(g_width, g_height, filename, buf);
    delete[] filename;
}


// -------------------------------------------------------------------
// Main
//
int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    
    calculate_inverse();        // inverse transform matrix only once
    
    render();
    
    saveFile();
    
    return 0;
}



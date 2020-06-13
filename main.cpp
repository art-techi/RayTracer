/*Student Name: Margaret(Meg) Coleman
 * Student ID: 
 * CSC305 Assignment 1
 * Last Updated: May 22, 2017
 *
 * Note: Most of this code is based from the books "Ray Tracing In One Weekend"
 *       and "Ray Tracing The Next Week" by Peter Shirley. For this assignment,
 *       I walked through the code in these two books and included the components
 *       were relevant to the requirements then applied the components I needed
 *       to build the scenes I wanted to present.
 *
 * Ray tracer: Creates different scenes and images by using vectors and rays
 *             Ouput is sent to a file called output.png in your current file
 *             directory.
 *
 * Use: uncomment (and comment out other) world variables in main to see effects.
 *      make sure to change the camera parameters to appropriate values or the
 *      scene will not appear.
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <float.h>
using namespace std;
class material;

/*=============================**VEC3**===================================*/
/*Class to hold a vector object. vec3 holds x,y,z components of a vector.
 * Operator functions that can be used on vectors declared here.
 * vec3 class and vector operations are as defined by Shirley in the Ray Tracing books
 */
class vec3{
public:
    vec3() {}
    vec3(float f0, float f1, float f2) { f[0]=f0; f[1]=f1; f[2]=f2; }
    //Note: inline places a copy of the code where the function is called
    //at compile time
    //vector can hold x,y,z coordinates or r,g,b values
    inline float x() const {return f[0];}
    inline float y() const {return f[1];}
    inline float z() const {return f[2];}
    inline float r() const {return f[0];}
    inline float g() const {return f[1];}
    inline float b() const {return f[2];}

    inline const vec3& operator+() const {return *this;}
    inline vec3 operator-() const {return vec3(-f[0], -f[1], -f[2]);}
    inline float operator[](int i) const {return f[i];}
    inline float& operator[](int i) {return f[i]; };

    //operator function set up
    inline vec3& operator+=(const vec3 &v2);
    inline vec3& operator-=(const vec3 &v2);
    inline vec3& operator*=(const vec3 &v2);
    inline vec3& operator/=(const vec3 &v2);
    inline vec3& operator*=(const float t);
    inline vec3& operator/=(const float t);

    //length, square length, unit vector set up
    inline float length() const {
        return sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]); }
    inline float squared_length() const {
        return f[0]*f[0]+f[1]*f[1]+f[2]*f[2]; }
    inline void make_unit_vector();

    float f[3];
};

/*=========================**VEC OPERATIONS**===========================*/
/*Set of vector operations that can be performed on vec3 objects built here.
 *You could also use the Eigen class as demonstrated in labs.
 */

//handle file io
inline std::istream& operator>>(std::istream &is, vec3 &t){
        is>> t.f[0]>>t.f[1]>>t.f[2];
        return is;
}
inline std::ostream& operator<<(std::ostream &os, const vec3 &t){
        os<<t.f[0]<<" "<<t.f[1]<<" "<<t.f[2];
     return os;
}
//create unit vector
inline void vec3::make_unit_vector(){
        float k=1.0/sqrt(f[0]*f[0]+f[1]*f[1]+f[2]*f[2]);
        f[0] *= k; f[1] *= k; f[2] *= k;
}
//vector addition
inline vec3 operator+(const vec3 &v1, const vec3 &v2){
        return vec3(v1.f[0] +v2.f[0], v1.f[1]+v2.f[1], v1.f[2]+v2.f[2]);
}
//vector subtraction
inline vec3 operator-(const vec3 &v1, const vec3 &v2){
        return vec3(v1.f[0] -v2.f[0], v1.f[1]-v2.f[1], v1.f[2]-v2.f[2]);
}
//vector multiplication
inline vec3 operator*(const vec3 &v1, const vec3 &v2){
        return vec3(v1.f[0] *v2.f[0], v1.f[1]*v2.f[1], v1.f[2]*v2.f[2]);
}
//vector division
inline vec3 operator/(const vec3 &v1, const vec3 &v2){
        return vec3(v1.f[0] /v2.f[0], v1.f[1]/v2.f[1], v1.f[2]/v2.f[2]);
}
//scalar multiplication
inline vec3 operator*(float t, const vec3 &v){
        return vec3(t*v.f[0], t*v.f[1], t*v.f[2]);
}
//scalar division
inline vec3 operator/(vec3 v, float t){
        return vec3(v.f[0]/t, v.f[1]/t, v.f[2]/t);
}
inline vec3 operator*(const vec3 &v, float t){
        return vec3(t*v.f[0], t*v.f[1], t*v.f[2]);
}

//dot product
inline float dot(const vec3 &v1, const vec3 &v2){
        return v1.f[0] *v2.f[0]+v1.f[1] *v2.f[1]+v1.f[2] *v2.f[2];
}
//cross product
inline vec3 cross(const vec3 &v1, const vec3 &v2){
        return vec3( (v1.f[1]*v2.f[2]-v1.f[2]*v2.f[1]),
                     (-(v1.f[0]*v2.f[2]-v1.f[2]*v2.f[0])),
                     (v1.f[0]*v2.f[1]-v1.f[1]*v2.f[0]));
}
//vector +=
inline vec3& vec3::operator+=(const vec3 &v){
        f[0]+=v.f[0];
        f[1]+=v.f[1];
        f[2]+=v.f[2];
        return *this;
}
//vector *=
inline vec3& vec3::operator*=(const vec3 &v){
        f[0]  *= v.f[0];
        f[1]  *= v.f[1];
        f[2]  *= v.f[2];
        return *this;
}
//vector /=
inline vec3& vec3::operator/=(const vec3 &v){
        f[0]/=v.f[0];
        f[1]/=v.f[1];
        f[2]/=v.f[2];
        return *this;
}
//vector -=
inline vec3& vec3::operator-=(const vec3& v){
        f[0]-=v.f[0];
        f[1]-=v.f[1];
        f[2]-=v.f[2];
        return *this;
}
//vector *= scalar
inline vec3& vec3::operator*=(const float t){
        f[0]*=t;
        f[1]*=t;
        f[2]*=t;
        return *this;
}
//vector /=scalar
inline vec3& vec3::operator/=(const float t){
        float k=1.0/t;
        f[0]*=k;
        f[1]*=k;
        f[2]*=k;
        return *this;
}
//create unit vector
inline vec3 unit_vector(vec3 v){
        return v/v.length();
}

/*==============================**RAY**================================*/
/* Class to hold ray object. ray holds the origin, direction, and time
 * parameters to calculate the ray equation p=P+t*Q where P is the origin
 * and Q is the direction of the ray.
 */
class ray{
public:
    ray(){}
    ray(const vec3& p, const vec3& q, float time=0.0) {P=p; Q=q; _time=time;}
    vec3 orig() const {return P;}
    vec3 dir() const {return Q;}
    float time() const {return _time;}
    vec3 pp(float t) const {return P+t*Q;}

    vec3 P;
    vec3 Q;
    float _time;
};

//get random suitable value for p and returns vec3 object
vec3 random_in_unit_disk(){
    vec3 p;
    do{
        p=2.0*vec3(drand48(), drand48(), 0)-vec3(1,1,0);
    }while(dot(p,p)>=1.0);
    return p;
}

/*============================**CAMERA**===============================*/
/*Class to hold camera functions. camera holds the point where we want to
 * look from, the point where we want to look at, the view up, the vfov,
 * the aspect, the aperture, the focus, and two time parameters.
 * Camera calculates where we are looking in the scene and returns a ray
 */
class cam{
public:
    cam(){}
    cam(vec3 camfrom, vec3 camat, vec3 vup, float vfov, float aspect, float aperture, float focus_dist, float t0, float t1){
        //time value for how long lense stays open
        time0=t0;
        time1=t1;
        origin=camfrom;
        lens_r=aperture/2;
        float theta=vfov*M_PI/180;
        float half_height=tan(theta/2);
        float half_width=aspect*half_height;
        //calculate the w,u,v vectors to create the axes for the camera
        w=unit_vector(camfrom-camat);
        u=unit_vector(cross(vup, w));
        v=cross(w,u);
        llc= origin - half_width*focus_dist*u -half_height*focus_dist*v - focus_dist*w;
        hori=2*half_width*focus_dist*u;
        vert=2*half_height*focus_dist*v;

    }
    ray get_ray(float s, float t){
        vec3 rd=lens_r*random_in_unit_disk();
        vec3 offset=u*rd.x()+v*rd.y();
        float time=time0+drand48()*(time1-time0);
        return ray(origin+offset, llc+s*hori+t*vert-origin-offset, time);}

    vec3 origin;
    vec3 llc;
    vec3 hori;
    vec3 vert;
    vec3 u,v,w;
    float lens_r;
    float time0, time1;
};

/*==========================**Hit_sphere CHECK**========================*/
/*check discriminant to see if sphere is hit or not
 * discriminant=0, not hit
 * discriminant=1, hit once
 * discriminant=2, hit twice
 */
float hit_sphere(const vec3& center, float radius, const ray& r){
    vec3 oc=r.orig()-center;
    float a=dot(r.dir(), r.dir());
    float b=2.0*dot(oc, r.dir());
    float c=dot(oc, oc)-radius*radius;
    float discr=b*b-4*a*c;
    if(discr<0){
        return -1.0;
    }else{
        return (-b-sqrt(discr))/(2.0*a);
    }
}

/*==============================**BOUNDS**=================================*/
inline float ffmin(float a, float b){return a<b?a:b;}
inline float ffmax(float a, float b){return a>b?a:b;}

/*creates boundary box to check if ray hits objects or not. Fuller description
 * in report
 */
class bound{
public:
    bound(){}
    bound(const vec3& p, const vec3& q){ _min=p; _max=q;}

    vec3 min() const {return _min;}
    vec3 max() const {return _max;}

    bool hit(const ray& r, float tmin, float tmax) const{
        for(int p=0; p<3; p++){
            float t0=ffmin((_min[p]-r.orig()[p])/r.dir()[p], (_max[p]-r.orig()[p])/r.dir()[p]);
            float t1=ffmax((_min[p]-r.orig()[p])/r.dir()[p], (_max[p]-r.orig()[p])/r.dir()[p]);
            tmin=ffmax(t0, tmin);
            tmax=ffmin(t1, tmax);
            if(tmax<=tmin)
                return false;
        }
        return true;
    }

    vec3 _min;
    vec3 _max;
};
//create box using min and max for checking if ray hits objects or not
bound surrounding_box(bound b0, bound b1){
    vec3 small(fmin(b0.min().x(),b1.min().x()), fmin(b0.min().y(), b1.min().y()), fmin(b0.min().z(), b1.min().z()));
    vec3 big(fmax(b0.max().x(), b1.max().x()), fmax(b0.max().y(), b1.max().y()), fmax(b0.max().z(), b1.max().z()));
    return bound(small, big);
}

/*============================**HITABLE**===============================*/
/*Struct to hold values for hitable. Holds three floats, point vector,
 * normal vector, and material type.
 */
struct hit_record{
    float t;
    float u;
    float v;
    vec3 p;
    vec3 normal;
    material *mat_ptr;
};

/*Hitable holds two virtual functions, hit and bounding_box
 */
class hitable{
public:
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const=0;
    virtual bool bounding_box(float t0, float t1, bound& box) const=0;
};

/*==========================**HITABLE_LIST**========================*/
/*hitable_list holds hitable objects and stores them in a list for easy access.
 */
class hitable_list: public hitable{
public:
    hitable_list(){}
    hitable_list(hitable **l, int n) {list=l; list_size=n;}
    virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, bound &box) const;
    hitable **list;
    int list_size;
};
//checks if rays hit objects in list via aabb box method
bool hitable_list::bounding_box(float t0, float t1, bound &box) const{
    if(list_size<1) return false;
    bound temp_box;
    bool first_ture=list[0]->bounding_box(t0, t1, temp_box);
    if(!first_ture)
        return false;
    else
        box=temp_box;
    for(int i=1; i<list_size; i++){
        if(list[0]->bounding_box(t0, t1, temp_box)){
            box=surrounding_box(box, temp_box);
        }
        else
            return false;
    }
    return true;
}
//checks if rays hit objects in list
bool hitable_list::hit(const ray& r, float t_min, float t_max, hit_record& rec) const{
    hit_record temp_rec;
    bool hityes=false;
    double closest=t_max;
    for(int i=0; i<list_size; i++){
        if(list[i]->hit(r, t_min, closest, temp_rec)){
            hityes=true;
            closest=temp_rec.t;
            rec=temp_rec;
        }
    }
    return hityes;
}


/* flip_normals is used in cornell functions. It inverts the normals so
 * we can see the background properly
 */
class flip_norm: public hitable{
public:
    flip_norm(hitable *p):ptr(p){}
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const{
        if(ptr->hit(r, t_min, t_max, rec)){
            rec.normal=-rec.normal;
            return true;
        }
        else
            return false;
    }
    virtual bool bounding_box(float t0, float t1, bound& box) const{
        return ptr->bounding_box(t0, t1, box);
    }
    hitable *ptr;
};

/*translate moves objects by the amount described by offset.
 */
class translate:public hitable{
public:
    translate(hitable *p, const vec3& displace):ptr(p), offset(displace) {}
    virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, bound& box) const;
    hitable *ptr;
    vec3 offset;
};
//check if rays hit object
bool translate::hit(const ray& r, float t_min, float t_max, hit_record& rec) const{
    ray moved_r(r.orig()-offset, r.dir(), r.time());
    if(ptr->hit(moved_r, t_min, t_max, rec)){
        rec.p+=offset;
        return true;
    }
    else
        return false;
}
//check if rays hit object via aabb box method
bool translate::bounding_box(float t0, float t1, bound& box) const{
    if(ptr->bounding_box(t0, t1, box)){
        box=bound(box.min()+offset, box.max()+offset);
        return true;
    }
    else
        return false;
}


/*============================**RECTANGLE**===============================*/
/*create a rectangle in either the x plane, y plane, or z plane depending
 * on which function the user chooses. Uses the box class to create a box with
 * one axis value at 0 to essentially form a two dimensional box, or a rectangle
 */
class xy_rectangle: public hitable{
public:
    xy_rectangle(){}
    xy_rectangle(float _x0, float _x1, float _y0, float _y1, float _k, material *mat): x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat){};
    virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, bound& box) const{
        box=bound(vec3(x0, y0, k-0.0001), vec3(x1, y1, k+0.0001));
        return true;
    }

    float x0, x1, y0, y1, k;
    material *mp;
};

class xz_rectangle: public hitable{
public:
    xz_rectangle(){}
    xz_rectangle(float _x0, float _x1, float _z0, float _z1, float _k, material *mat): x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat){};
    virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, bound& box) const{
        box=bound(vec3(x0, k-0.0001,z0), vec3(x1, k+0.001, z1));
        return true;
    }

    float x0, x1, z0, z1, k;
    material *mp;
};

class yz_rectangle: public hitable{
public:
    yz_rectangle(){}
    yz_rectangle(float _y0, float _y1, float _z0, float _z1, float _k, material *mat): y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {};
    virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, bound& box) const{
        box=bound(vec3(k-0.0001, y0, z0), vec3(k+0.0001, y1, z1));
        return true;
    }

    float y0, y1, z0, z1, k;
    material *mp;
};
//check if rays hit rectangles created
bool xy_rectangle::hit(const ray& r, float t0, float t1, hit_record& rec) const{
    float t=(k-r.orig().z())/r.dir().z();
    if(t<t0||t>t1)
        return false;
    float x=r.orig().x()+t*r.dir().x();
    float y=r.orig().y()+t*r.dir().y();
    if(x<x0||x>x1||y<y0||y>y1)
        return false;
    rec.u=(x-x0)/(x1-x0);
    rec.v=(y-y0)/(y1-y0);
    rec.t=t;
    rec.mat_ptr=mp;
    rec.p=r.pp(t);
    rec.normal=vec3(0,0,1);
    return true;
}

bool xz_rectangle::hit(const ray& r, float t0, float t1, hit_record& rec) const{
    float t=(k-r.orig().y())/r.dir().y();
    if(t<t0||t>t1)
        return false;
    float x=r.orig().x()+t*r.dir().x();
    float z=r.orig().z()+t*r.dir().z();
    if(x<x0||x>x1||z<z0||z>z1)
        return false;
    rec.u=(x-x0)/(x1-x0);
    rec.v=(z-z0)/(z1-z0);
    rec.t=t;
    rec.mat_ptr=mp;
    rec.p=r.pp(t);
    rec.normal=vec3(0,1,0);
    return true;
}

bool yz_rectangle::hit(const ray& r, float t0, float t1, hit_record& rec) const{
    float t=(k-r.orig().x())/r.dir().x();
    if(t<t0||t>t1)
        return false;
    float y=r.orig().y()+t*r.dir().y();
    float z=r.orig().z()+t*r.dir().z();
    if(y<y0||y>y1||z<z0||z>z1)
        return false;
    rec.u=(y-y0)/(y1-y0);
    rec.v=(z-z0)/(z1-z0);
    rec.t=t;
    rec.mat_ptr=mp;
    rec.p=r.pp(t);
    rec.normal=vec3(1,0,0);
    return true;
}


/*==========================**SPHERE**========================*/
/*sphere class holds the center and the radius of the sphere and calculates
 * the sphere equation. Checks the discriminant to see if the sphere is
 * hit by a ray or not.
 */
class sphere: public hitable{
public:
    sphere(){}
    sphere(vec3 cen, float r, material *m): cen(cen), rad(r), mat_ptr(m) {};
    virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, bound &box) const;
    vec3 cen;
    float rad;
    material *mat_ptr;
};
//check if ray hits sphere via aabb box method
bool sphere::bounding_box(float t0, float t1, bound &box) const{
    box=bound(cen-vec3(rad, rad, rad), cen+vec3(rad, rad, rad));
    return true;
}
//check if ray hits sphere
bool sphere:: hit(const ray& r, float t_min, float t_max, hit_record& rec) const{
    vec3 oc=r.orig()-cen;
    float a=dot(r.dir(), r.dir());
    float b=dot(oc, r.dir());
    float c=dot(oc, oc)-rad*rad;
    float discr=b*b-a*c;
    if(discr>0){
        float temp=(-b-sqrt(discr))/a;
        if(temp<t_max && temp>t_min){
            rec.t=temp;
            rec.p=r.pp(rec.t);
            rec.normal=(rec.p-cen)/rad;
            rec.mat_ptr=mat_ptr;
            return true;
        }
        temp=(-b+sqrt(discr))/a;
        if(temp<t_max && temp>t_min){
            rec.t=temp;
            rec.p=r.pp(rec.t);
            rec.normal=(rec.p-cen)/rad;
            rec.mat_ptr=mat_ptr;
            return true;
        }
    }
    return false;
}


/*============================**TEXTURE**===============================*/
/*texture class
 */
class texture{
public:
    virtual vec3 value(float u, float v, const vec3& p) const=0;
};
//create constant colour texture-nothing special just substitute for basic colour vec3
class plain_texture:public texture{
public:
    plain_texture(){}
    plain_texture(vec3 c): colour(c){}
    virtual vec3 value(float u, float v, const vec3& p) const{
        return colour;
    }

    vec3 colour;
};
//create checker texture
class checker_texture: public texture{
public:
    checker_texture(){}
    checker_texture(texture *t0, texture *t1): one(t0), two(t1){}
    virtual vec3 value(float u, float v, const vec3& p) const{
        /*can change 0.5 value depending on size of sphere
          the example by Shirley used a sphere with radius 1000 and value 5
          for my example since my sphere is of radius 100, I use 0.5
        */
        float signs=sin(0.5*p.x())*sin(0.5*p.y())*sin(0.5*p.z());
        //alternate signs to get checker effect
        if(signs<0)
            return two->value(u,v,p);
        else
            return one->value(u,v,p);
    }

    texture *one;
    texture *two;
};
//create noise texture-used for perlin

/*============================**BOX**===============================*/
/*box class uses aabb class to create a box using min and max values.
 * box constructor creates a hitable list holding rectangles and flip_normals
 * objects holding the coordinates of the box to be made.
 */
class box: public hitable{
public:
    box(){}
    box(const vec3& p0, const vec3& p1, material *ptr);
    virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, bound& box) const{
        box=bound(pmin, pmax);
        return true;
    }
    vec3 pmin, pmax;
    hitable *list_ptr;
};
//box constructor takes in two vectors and material type
//creates box based on min and max vector values and applies given material
box::box(const vec3 &p0, const vec3 &p1, material *ptr){
    pmin=p0;
    pmax=p1;
    hitable **list=new hitable*[6];
    list[0]=new xy_rectangle(p0.x(), p1.x(), p0.y(), p1.y(),p1.z(), ptr);
    list[1]=new flip_norm(new xy_rectangle(p0.x(),p1.x(),p0.y(),p1.y(),p0.z(), ptr));
    list[2]=new xy_rectangle(p0.x(),p1.x(),p0.z(),p1.z(),p1.y(), ptr);
    list[3]=new flip_norm(new xz_rectangle(p0.x(),p1.x(),p0.z(),p1.z(),p0.y(),ptr));
    list[4]=new yz_rectangle(p0.y(),p1.y(),p0.z(),p1.z(),p1.x(),ptr);
    list[5]=new flip_norm(new yz_rectangle(p0.y(),p1.y(),p0.z(),p1.z(),p0.x(),ptr));
    list_ptr=new hitable_list(list,6);
}
//check if ray hits box
bool box::hit(const ray& r, float t0, float t1, hit_record& rec) const{
    return list_ptr->hit(r, t0, t1, rec);
}

/*==========================**MATERIAL**========================*/
//calculates polynomial approximation of glass reflectivity based on angle
float schlickApprox(float cosine, float ref_idx){
    float r0=(1-ref_idx)/(1+ref_idx);
    r0=r0*r0;
    return r0+(1-r0)*pow((1-cosine),5);
}
//calculates amount to refract light
bool refract(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted){
    vec3 uv=unit_vector(v);
    float dt=dot(uv, n);

    float discr=1.0-ni_over_nt*ni_over_nt*(1-dt*dt);
    if(discr>0){
        refracted=ni_over_nt*(uv-n*dt)-n*sqrt(discr);
        return true;
    }else
        return false;
}
//calculates amount to reflect light
vec3 reflect(const vec3& v, const vec3& n){
    return v-2*dot(v,n)*n;
}
//finds suitable random p value
vec3 random_in_unit_sphere(){
    vec3 p;
    do{
        p=2.0*vec3(drand48(), drand48(), drand48())-vec3(1,1,1);
    }while(p.squared_length()>=1.0);
    return p;
}
//holds scattered and emitted functionality for other material classes
class material{
public:
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const=0;
    virtual vec3 emitted(float u, float v, const vec3& p) const{ return vec3(0,0,0);}
};
//creates light-like effect for scenes with light source
class diffuse_light: public material{
public:
    diffuse_light(texture *a):emit(a){}
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const{return false;}
    virtual vec3 emitted(float u, float v, const vec3& p) const{return emit->value(u,v,p);}
    texture *emit;
};
//creates diffuse quality for shapes
class lambertian: public material{
public:
    lambertian(texture *a): albedo(a){}
    //lambertian(const vec3& b): albedo2(b){}
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const{
        vec3 target=rec.p+rec.normal+random_in_unit_sphere();
        scattered=ray(rec.p, target-rec.p, r_in.time());

        attenuation=albedo->value(0,0,rec.p);
        //attenuation=albedo->value(rec.p, target-rec.p, r_in.time());
        //attenuation=albedo;

        return true;
    }

    texture *albedo;
    //vec3 albedo;
};
//creates metal quality for shapes
class metal: public material{
public:
    metal(const vec3& a, float f): albedo(a){if (f<1) fuzz=f; else fuzz=1;}
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const{
        vec3 reflected=reflect(unit_vector(r_in.dir()), rec.normal);
        scattered=ray(rec.p, reflected+fuzz*random_in_unit_sphere());
        attenuation=albedo;
        return (dot(scattered.dir(), rec.normal)>0);
    }

    vec3 albedo;
    float fuzz;
};
//creates glass quality for shapes
class dielectric: public material{
public:
    dielectric(float ri): ref_idx(ri) {}
    virtual bool scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered) const{
        vec3 outward_normal;
        vec3 reflected=reflect(r_in.dir(), rec.normal);
        float ni_over_nt;
        attenuation=vec3(1.0, 1.0, 1.0);
        vec3 refracted;
        float reflect_prob;
        float cosine;
        if(dot(r_in.dir(), rec.normal)>0){
            outward_normal=-rec.normal;
            ni_over_nt=ref_idx;
            cosine=dot(r_in.dir(), rec.normal)/r_in.dir().length();
            cosine=sqrt(1-ref_idx*ref_idx*(1-cosine*cosine));
        }else{
            outward_normal=rec.normal;
            ni_over_nt=1.0/ref_idx;
            cosine=-dot(r_in.dir(), rec.normal)/r_in.dir().length();
        }
        if(refract(r_in.dir(), outward_normal, ni_over_nt, refracted)){
            reflect_prob=schlickApprox(cosine, ref_idx);
        }else{
            reflect_prob=1.0;
        }
        if(drand48()<reflect_prob){
            scattered=ray(rec.p, reflected);
        }else{
            scattered=ray(rec.p, refracted);
        }
        return true;
    }
    float ref_idx;
};

/*============================**COLOUR & SCENES**===============================*/
/*colour function holds a ray, the world we are using, and the depth.
 * calculates a vector holding rgb values and sends back to inner loop in main
 *
 * Depending on what scene we are creating in main, some of the lines will need
 * to be uncommented and others commented. This program should currently be set up
 * for random_scene() and random_scene_move()
 */
vec3 colour(const ray& r, hitable *world, int depth){
    hit_record rec;
    if(world->hit(r, 0.001, MAXFLOAT, rec)){
        ray scattered;
        vec3 attenuation;
        //emitted allows for the "glow" of the light to be present in scenes
        vec3 emitted=rec.mat_ptr->emitted(rec.u,rec.v,rec.p);
        if(depth<50 && rec.mat_ptr->scatter(r, rec, attenuation,scattered)){
            return emitted+attenuation*colour(scattered,world,depth+1);
        }else{
            return emitted;
        }
    }else{
        vec3 unit_direction=unit_vector(r.dir());
        float t=0.5*(unit_direction.y()+1.0);
        //linear blend
        return (1.0-t)*vec3(1.0, 1.0, 1.0)+t*vec3(0.5, 0.7, 1.0);
    }
}

//two spheres on a plane with one light source
hitable *basic_req(){
    //colour and texture
    texture *plain=new plain_texture(vec3(0.3,0.2,0.4));
    material *white=new lambertian(new plain_texture(vec3(0.73,0.73, 0.73)));
    material *light=new diffuse_light(new plain_texture(vec3(4,4,4))); //above 4 for brightness

    //create list of items for scene
    int n=50;
    hitable **list=new hitable*[n+1];
    list[0]=new sphere(vec3(3,5,10),5,new lambertian(plain));
    list[1]=new sphere(vec3(10, 5, 10), 2, new lambertian(plain));
    list[2]=new xy_rectangle(10,20,5,10,-2, light);
    list[3]=new xz_rectangle(-555,555,-555,555,0,white);

    return new hitable_list(list, 4);
}

//multiple spheres with multiple textures in a cornell box
hitable *advanced_req(){
    hitable **list=new hitable*[10];
    int i=0;

    //colours and textures
    material *midnight=new lambertian(new plain_texture(vec3(0.0,0.0,0.8)));
    material *white=new lambertian(new plain_texture(vec3(0.73,0.73, 0.73)));
    material *teal=new lambertian(new plain_texture(vec3(0.0, 0.45, 0.45)));
    material *red=new lambertian(new plain_texture(vec3(0.3, 0, 0)));
    material *yellow=new lambertian(new plain_texture(vec3(0.3,0.3,0)));
    material *light=new diffuse_light(new plain_texture(vec3(4,4,4)));
    texture *checker=new checker_texture(new plain_texture(vec3(0.2,0.3,0.1)),new plain_texture(vec3(0.9,0.9,0.9)));

    //walls and light source cornell box
    list[i++]=new flip_norm(new yz_rectangle(0,555,0,555,555,teal));
    list[i++]=new yz_rectangle(0,555,0,555,0,midnight);
    list[i++]=new xz_rectangle(213,343,227,332,554,light);
    list[i++]=new flip_norm(new xz_rectangle(0,555,0,555,555,red));
    list[i++]=new xz_rectangle(0,555,0,555,0,white);
    list[i++]=new flip_norm(new xy_rectangle(0,555,0,555,555,yellow));

    //spheres
    hitable *boundary=new sphere(vec3(160, 100, 145), 100, new lambertian(checker));
    list[i++]=boundary;
    hitable *boundary2=new sphere(vec3(320, 100, 330), 100, new dielectric(1.5));
    list[i++]=boundary2;
    hitable *boundary3=new sphere(vec3(400, 450, 400), 100, new metal(vec3(0.3, 0.4, 0.6), 0.0));
    list[i++]=boundary3;
    hitable *boundary4=new sphere(vec3(200, 450, 400), 100, new lambertian(new plain_texture(vec3(0.45, 0.2, 0.3))));
    list[i++]=boundary4;

    return new hitable_list(list, i);
}

/*===================================**MAIN**=======================================*/
/*IMPORTANT: Make sure screen size and camera parameters are set properly for each scene!!
 * The code is currently set up to work for basic_req()
 * To use other scenes, comment and uncomment appropriate lines (as labelled below)
 * and make sure to change the camera lines in the camera method before running the code
 *
 * Outputs image to output.png in your current file directory
 */
int main()
{
    //set output image to file called "output.png" in current directory
    ofstream my_out_file;
    my_out_file.open("output.png");

    if(my_out_file.is_open()){
        cout <<"Producing image...May take a few mins. Please wait..." <<endl;

        //screen size
        int nx=800; //width
        int ny=800; //height
        int ns=100;


        my_out_file <<"P3\n" << nx << " " << ny << "\n255" << endl;

        /*uncomment one of the worlds below (comment the other) to see effects*/
        hitable *world=basic_req();           //relatively quick
        //hitable *world=advanced_req();      //takes longer to generate

        /*NOTE: it is super important to change camera parameters or you won't see anything...*/

        vec3 camfrom(-30,10,-40);          //basic scene -adjustable
        vec3 camat(0,0,0);

        //vec3 camfrom(278, 278, -800);    //advanced -adjustable
        //vec3 camat(278,278,0);

        float dist_to_focus=10.0;
        float aperture=0.0;
        float vfov=40.0;

        cam camera(camfrom, camat, vec3(0,1,0), vfov, float(nx)/float(ny), aperture, dist_to_focus, 0.0, 1.0);

        //for each pixel, get the ray, the point, the colour value, output rgb in array
        for(int j=ny-1; j>=0; j--){
            for(int i=0; i<nx; i++){
                vec3 col(0,0,0);
                for(int s=0; s<ns; s++){
                    float u=float(i+drand48())/float(nx);
                    float v=float(j+drand48())/float(ny);
                    ray r=camera.get_ray(u,v);
                    vec3 p=r.pp(2.0);
                    col+=colour(r, world, 0);
                }
                col/=float(ns);
                col=vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
                int ir=int(255.99*col[0]);
                int ig=int(255.99*col[1]);
                int ib=int(255.99*col[2]);
                //send values to the file
                my_out_file << ir << " " << ig << " " << ib << endl;
            }
        }
        cout << "Done generating image." << endl;
    }else{
        cout << "Sorry. No output was produced. Check file I/O." <<endl;
    }
}

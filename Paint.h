#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <utility>

const float PI = 3.1415926535897932384626433832795;
const float RAD_TO_DEG = 180/PI;
const float DEG_TO_RAD = PI/180;
const float defaultDepth = 2;

class Vector2;
class Vector3;
class Vector4;
class Quaternion;
class Matrix4x4;
class Transform;
class Mesh;
class Object;
class Scene;

class Vector2
{
public:
	float x, y;
	Vector2(float a, float b) { x = a; y = b; }
	Vector2() :Vector2(0, 0) {}
	float & operator[](int i) {
        return reinterpret_cast<float *>(this)[i];
	}
    const float & operator[](int i) const {
        return reinterpret_cast<const float *>(this)[i];
    }
	friend Vector2 operator+(const Vector2 & a, const Vector2 & b) { return {a.x + b.x, a.y + b.y}; }
    friend Vector2 operator-(const Vector2 & a, const Vector2 & b) { return {a.x - b.x, a.y - b.y}; }
	friend Vector2 operator*(const Vector2 & a, float b) { return {a.x*b, a.y*b}; }
	friend Vector2 operator*(float b, const Vector2 & a) { return {a.x*b, a.y*b}; }
    friend Vector2 operator/(const Vector2 & a, float b) { return {a.x/b, a.y/b}; }
	Vector2 operator+() const { return *this; }
	Vector2 operator-() const { return *this * -1; }

	static float dot(const Vector2 & a, const Vector2 & b) { return a.x * b.x + a.y * b.y; }
	float magnitude() const { return sqrt(x * x + y * y); }
	Vector2 normalized() const { float m = this->magnitude(); return *this / m;}
	void normalize() { float m = this->magnitude(); x/=m; y/=m; }
};

class Vector3
{
public:
	float x, y, z;
	Vector3(float a, float b, float c) { x = a; y = b; z = c; }
	Vector3() : Vector3(0, 0, 0) {}
    float & operator[](int i) {
        return reinterpret_cast<float *>(this)[i];
    }
    const float & operator[](int i) const {
        return reinterpret_cast<const float *>(this)[i];
    }
    const Vector3 & operator+= (const Vector3 & a) { x+= a.x; y += a.y; z += a.z; return *this; }
	friend Vector3 operator+(const Vector3 & a, const Vector3 & b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
    friend Vector3 operator-(const Vector3 & a, const Vector3 & b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
	friend Vector3 operator*(const Vector3 & a, float b) { return {a.x*b, a.y*b, a.z*b}; }
	friend Vector3 operator*(float b, const Vector3& a) { return {a.x*b, a.y*b, a.z*b}; }
	friend Vector3 operator/(const Vector3& a, float b) { return {a.x/b, a.y/b, a.z/b}; }
	Vector3 operator+() const { return *this; }
	Vector3 operator-() const { return *this * -1; }


	static Vector3 cross(const Vector3 & a, const Vector3 & b) { return {a.y * b.z - b.y * a.z, a.x * b.z - b.x * a.z, a.x * b.y - b.x * a.y}; }
	static float dot(const Vector3 & a, const Vector3 & b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
	float magnitude() const { return sqrt(x * x + y * y + z * z); }
	Vector3 normalized() const { float m = this->magnitude();  return *this * (1 / m); }
    void normalize() { float m = this->magnitude(); x/=m; y/=m; z/=m; }
	friend std::ostream& operator<< (std::ostream & out, const Vector3 & v)
	{
		out << v.x << " " << v.y << " " << v.z;
		return out;
	}
	friend std::istream& operator>> (std::istream & in, Vector3 & v)
	{
		in >> v.x >> v.y >> v.z;
		return in;
	}
};

class Vector4
{
public:
	float x, y, z, w;
	Vector4(float a, float b, float c, float d) { x = a; y = b; z = c; w = d; }
	Vector4() : Vector4(0, 0, 0, 0) {}
    float & operator[](int i) {
        return reinterpret_cast<float *>(this)[i];
    }
    const float & operator[](int i) const {
        return reinterpret_cast<const float *>(this)[i];
    }
	friend Vector4 operator+(const Vector4 & a, const Vector4& b) { return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w}; }
    friend Vector4 operator-(const Vector4 & a, const Vector4& b) { return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w}; }
	friend Vector4 operator*(const Vector4& a, float b) { return {a.x*b, a.y*b, a.z*b, a.w*b}; }
	friend Vector4 operator*(float b, const Vector4& a) { return {a.x*b, a.y*b, a.z*b, a.w*b}; }
    friend Vector4 operator/(const Vector4& a, float b) { return {a.x/b, a.y/b, a.z/b, a.w/b}; }
	Vector4 operator+() const { return *this; }
	Vector4 operator-() const { return *this * -1; }

	friend std::ostream& operator<< (std::ostream& out, const Vector4& v)
	{
		out << v.x << " " << v.y << " " << v.z << " " << v.w;
		return out;
	}

	float magnitude() const { return sqrt(x * x + y * y + z * z + w * w); }
	Vector4 normalized() { float m = this->magnitude();  return *this / m; }
	void normalize() { float m = this->magnitude(); x/=m; y /=m; z/=m; w/=m; }
};

class Quaternion {
private:
    static void makePositive(Vector3 &euler) {
        float num1 = -9.f / (500.f * PI);
        float num2 = 360.f + num1;
        if (euler.x < num1)
            euler.x += 360.f;
        else if (euler.x > num2)
            euler.x -= 360.f;
        if (euler.y < num1)
            euler.y += 360.f;
        else if (euler.y > num2)
            euler.y -= 360.f;
        if (euler.z < num1)
            euler.z += 360.f;
        else if (euler.z > num2)
            euler.z -= 360.f;
    }

    float magnitude() const {
        return w*w + x*x + y*y + z*z;
    }
public:
    float w, x, y, z;

    Quaternion() {
        x = 0;
        y = 0;
        z = 0;
        w = 1;
    }

    Quaternion(float d, float a, float b, float c) {
        x = a;
        y = b;
        z = c;
        w = d;
    }

    static Vector3 toEulerAngles(const Quaternion & q) {
        Vector3 angles;
        // x-axis rotation
        float sinr_cosp = 2 * (q.w * q.x + q.y * q.z);
        float cosr_cosp = 1 - 2 * (q.x * q.x + q.y * q.y);
        angles.x = std::atan2(sinr_cosp, cosr_cosp);
        // y-axis rotation
        float sinp = 2 * (q.w * q.y - q.z * q.x);
        if (std::abs(sinp) >= 1)
            angles.y = std::copysign(M_PI / 2, sinp); // use 90 degrees if out of range
        else
            angles.y = std::asin(sinp);
        // z-axis rotation
        float siny_cosp = 2 * (q.w * q.z + q.x * q.y);
        float cosy_cosp = 1 - 2 * (q.y * q.y + q.z * q.z);
        angles.z = std::atan2(siny_cosp, cosy_cosp);
        angles = angles * RAD_TO_DEG;
        makePositive(angles);
        return angles;
    }

    static Quaternion fromEulerAngles(const Vector3 & angles) {
        return Quaternion::fromEulerAngles(angles.x, angles.y, angles.z);
    }

    static Quaternion fromEulerAngles(float x, float y, float z) {
        x *= DEG_TO_RAD;
        y *= DEG_TO_RAD;
        z *= DEG_TO_RAD;
        // Abbreviations for the various angular functions
        float cy = cos(z * 0.5);
        float sy = sin(z * 0.5);
        float cp = cos(y * 0.5);
        float sp = sin(y * 0.5);
        float cr = cos(x * 0.5);
        float sr = sin(x * 0.5);
        Quaternion q;
        q.w = cy * cp * cr + sy * sp * sr;
        q.x = cy * cp * sr - sy * sp * cr;
        q.y = sy * cp * sr + cy * sp * cr;
        q.z = sy * cp * cr - cy * sp * sr;
        return q;
    }

    static Quaternion angleAxis(const Vector3 &axis, float angle) {
        Vector3 axisNormalized = axis.normalized();
        float c = cos(angle / 2), s = sin(angle / 2);
        return {c, axisNormalized.x * s, axisNormalized.y * s, axisNormalized.z * s};
    }

    static Quaternion inverse(const Quaternion & q) {
        float m = q.magnitude();
        return {q.w/m, -q.x/m, -q.y/m, -q.z/m};
    }

    static const Quaternion identity;
    const Quaternion & operator*=(const Quaternion & q) {
        this->w = w * q.w - x * q.x - y * q.y - z * q.z;
        this->x = w * q.x + x * q.w + y * q.z - z * q.y;
        this->y = w * q.y - x * q.z + y * q.w + z * q.x;
        this->z = w * q.z + x * q.y - y * q.x + z * q.w;
        return *this;
    }
    friend Quaternion operator*(const Quaternion & q1, const Quaternion & q2) {
        return {q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z,
                q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
                q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
                q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w,
        };
    }
};

class Matrix4x4
{
private:
	float m[16];
public:
	Matrix4x4()
	{
		m[0] = m[5] = m[10] = m[15] = 1.0;
		m[1] = m[2] = m[3] = m[4] = 0.0;
		m[6] = m[7] = m[8] = m[9] = 0.0;
		m[11] = m[12] = m[13] = m[14] = 0.0;
	}
	float & operator[](int i)
	{
		return m[i];
	}
    const float & operator[](int i) const
    {
        return m[i];
    }
	friend std::ostream& operator<< (std::ostream& out, const Matrix4x4 & m)
	{
		for (int i = 0; i < 4; i++)
		{
			out << m[0+i] << " " << m[4+i] << " " << m[8+i] << " " << m[12+i] << std::endl;
		}
		return out;
 	}
	friend Matrix4x4 operator* (const Matrix4x4 & m1, const Matrix4x4 & m2)
	{
		Matrix4x4 result;
		// Fisrt Column
		result[0] = m1[0] * m2[0] + m1[4] * m2[1] + m1[8] * m2[2] + m1[12] * m2[3];
		result[1] = m1[1] * m2[0] + m1[5] * m2[1] + m1[9] * m2[2] + m1[13] * m2[3];
		result[2] = m1[2] * m2[0] + m1[6] * m2[1] + m1[10] * m2[2] + m1[14] * m2[3];
		result[3] = m1[3] * m2[0] + m1[7] * m2[1] + m1[11] * m2[2] + m1[15] * m2[3];

		// Second Column
		result[4] = m1[0] * m2[4] + m1[4] * m2[5] + m1[8] * m2[6] + m1[12] * m2[7];
		result[5] = m1[1] * m2[4] + m1[5] * m2[5] + m1[9] * m2[6] + m1[13] * m2[7];
		result[6] = m1[2] * m2[4] + m1[6] * m2[5] + m1[10] * m2[6] + m1[14] * m2[7];
		result[7] = m1[3] * m2[4] + m1[7] * m2[5] + m1[11] * m2[6] + m1[15] * m2[7];

		// Third Column
		result[8] = m1[0] * m2[8] + m1[4] * m2[9] + m1[8] * m2[10] + m1[12] * m2[11];
		result[9] = m1[1] * m2[8] + m1[5] * m2[9] + m1[9] * m2[10] + m1[13] * m2[11];
		result[10] = m1[2] * m2[8] + m1[6] * m2[9] + m1[10] * m2[10] + m1[14] * m2[11];
		result[11] = m1[3] * m2[8] + m1[7] * m2[9] + m1[11] * m2[10] + m1[15] * m2[11];

		// Fourth Column
		result[12] = m1[0] * m2[12] + m1[4] * m2[13] + m1[8] * m2[14] + m1[12] * m2[15];
		result[13] = m1[1] * m2[12] + m1[5] * m2[13] + m1[9] * m2[14] + m1[13] * m2[15];
		result[14] = m1[2] * m2[12] + m1[6] * m2[13] + m1[10] * m2[14] + m1[14] * m2[15];
		result[15] = m1[3] * m2[12] + m1[7] * m2[13] + m1[11] * m2[14] + m1[15] * m2[15];

		return result;
	}
	friend Vector4 operator*(const Matrix4x4 & m1, const Vector4 & v)
	{
		Vector4 result;
		result[0] = m1[0] * v[0] + m1[4] * v[1] + m1[8] * v[2] + m1[12] * v[3];
		result[1] = m1[1] * v[0] + m1[5] * v[1] + m1[9] * v[2] + m1[13] * v[3];
		result[2] = m1[2] * v[0] + m1[6] * v[1] + m1[10] * v[2] + m1[14] * v[3];
		result[3] = m1[3] * v[0] + m1[7] * v[1] + m1[11] * v[2] + m1[15] * v[3];
		return result;
	}
	Matrix4x4 Transposed()
	{
		Matrix4x4 m2;

		//First Column
		m2[0] = m[0];
		m2[1] = m[4];
		m2[2] = m[8];;
		m2[3] = m[12];

		//Second Column
		m2[4] = m[1];
		m2[5] = m[4];
		m2[6] = m[9];;
		m2[7] = m[13];

		//Third Column
		m2[8] = m[2];
		m2[9] = m[6];
		m2[10] = m[10];
		m2[11] = m[14];

		//Fourth Column
		m2[12] = m[3];
		m2[13] = m[7];
		m2[14] = m[11];
		m2[15] = m[15];
		return m2;
	}

	static Matrix4x4 Translate(const Vector3 & pos)
	{
		Matrix4x4 matrix;
		matrix[12] = pos.x;
		matrix[13] = pos.y;
		matrix[14] = pos.z;
		return matrix;
	}
	static Matrix4x4 Rotate(const Quaternion & q)
	{
		Matrix4x4 matrix;
		// First Column
		matrix[0] = 1 - 2 * (q.y * q.y + q.z * q.z);
		matrix[1] = 2 * (q.x * q.y + q.z * q.w);
		matrix[2] = 2 * (q.x * q.z - q.y * q.w);
		matrix[3] = 0;

		// Second Column
		matrix[4] = 2 * (q.x * q.y - q.z * q.w);
		matrix[5] = 1 - 2 * (q.x * q.x + q.z * q.z);
		matrix[6] = 2 * (q.z * q.y + q.x * q.w);
		matrix[7] = 0;

		// Third Column
		matrix[8] = 2 * (q.x * q.z + q.y * q.w);
		matrix[9] = 2 * (q.y * q.z - q.x * q.w);
		matrix[10] = 1 - 2 * (q.x * q.x + q.y * q.y);
		matrix[11] = 0;

		// Fourth Column
		matrix[12] = 0;
		matrix[13] = 0;
		matrix[14] = 0;
		matrix[15] = 1;
		return matrix;
	}
	static Matrix4x4 Scale(const Vector3 & s)
	{
		Matrix4x4 matrix;
		matrix[0] = s.x;
		matrix[5] = s.y;
		matrix[10] = s.z;
		return matrix;
	}
	static Matrix4x4 TRS(const Vector3 & pos, const Quaternion & q, const Vector3 & s)
	{
		//fout << Matrix4x4::Scale(s);
		//fout << Matrix4x4::Rotate(q);
		//fout << Matrix4x4::Translate(pos);
		return Matrix4x4::Translate(pos) * Matrix4x4::Rotate(q) * Matrix4x4::Scale(s);
	}

	static Matrix4x4 Perspective(float fov = PI / 3, float aspect = 16. / 9, float zNear = 0.3, float zFar = 1000)
	{
		Matrix4x4 matrix;
		float tan = tanf(fov / 2.0);
		//float size = zNear * tanf(fov / 2.0);
		//float left = -size, right = size, bottom = -size / aspect, top = size / aspect;

		// First Column
		matrix[0] = 1/(aspect * tan);
		matrix[1] = 0.0;
		matrix[2] = 0.0;
		matrix[3] = 0.0;

		// Second Column
		matrix[4] = 0.0;
		matrix[5] = 1/tan;
		matrix[6] = 0.0;
		matrix[7] = 0.0;

		// Third Column
		matrix[8] = 0;
		matrix[9] = 0;
		matrix[10] = (zNear + zFar) / (zFar-zNear);
		matrix[11] = 1;

		// Fourth Column
		matrix[12] = 0.0;
		matrix[13] = 0.0;
		matrix[14] = (-2 * zFar * zNear) / (zFar-zNear);
		matrix[15] = 0.0;
		return matrix;
	}
};

class Color {
public:
    unsigned char r, g, b, a;
    Color(unsigned char r, unsigned char g, unsigned char b, unsigned char a) {
        this->r = r; this->g = g; this->b = b; this->a = a; }
    Color() : Color(0,0,0, 255) {}

    friend Color operator*(const Color & c, float k) {
        return {static_cast<unsigned char>(c.r*k),
                static_cast<unsigned char>(c.g*k),
                static_cast<unsigned char>(c.b*k),
                static_cast<unsigned char>(c.a*k)};
    }
    friend Color operator*(float k, const Color& c) {
        return {static_cast<unsigned char>(c.r*k),
                static_cast<unsigned char>(c.g*k),
                static_cast<unsigned char>(c.b*k),
                static_cast<unsigned char>(c.a*k)};
    }
};

class Transform
{
private:
public:
    Vector3 position = {0, 0, 0};
    Quaternion rotation;
    Vector3 scale = {1, 1, 1};

    Vector3 eulerAngles() const { return Quaternion::toEulerAngles(rotation); }
    void eulerAngles(const Vector3 & angles) { this->rotation = Quaternion::fromEulerAngles(angles); }
    void eulerAngles(float x, float y, float z) { this->rotation = Quaternion::fromEulerAngles(x, y, z); }

    void rotate(Vector3 eulerAngles, bool relativeToItself = true)
    {
        Quaternion quaternion = Quaternion::fromEulerAngles(eulerAngles.x, eulerAngles.y, eulerAngles.z);
        if (relativeToItself)
            this->rotation *= quaternion;
        else
            this->rotation = quaternion * this->rotation;
    }
};

class Mesh
{
public:
    std::vector<Vector3> vertices;
    std::vector<int> triangles;

    Mesh() = default;
    explicit Mesh(const char * str) : Mesh(std::string(str)) {}
    explicit Mesh(const std::string & str) {
        Vector3 vertex;
        std::ifstream fin;
        fin.open(str + "/vertices.txt");
        while (fin >> vertex) {
            this->vertices.push_back(vertex);
        }
        fin.close();
        fin.open(str + "/triangles.txt");
        int index;
        while ( fin >> index) {
            this->triangles.push_back(index);
        }
        fin.close();
    }
};

class Object
{
private:
    Color color = {200, 200, 200, 255};
    const Mesh *mesh = nullptr;
public:
    Transform transform;

    Object() = default;
    explicit Object(const Mesh & mesh) {
        this->mesh = &mesh;
    }
    explicit Object(const Mesh * mesh) {
        this->mesh = mesh;
    }
    void setColor(int r, int g, int b) {
        color.r = r;
        color.g = g;
        color.b = b;
    }
    friend class Camera;
};

class Scene {
private:
    std::vector<const Object*> objects;
public:
    void addToScene(const Object & obj) {
        this->objects.push_back(&obj);
    }
    void addToScene(const Object * obj) {
        this->objects.push_back(obj);
    }

    void removeFromScene(const Object & obj) {
        this->objects.erase(std::find(this->objects.begin(), this->objects.end(), &obj));
    }

    void removeFromScene(const Object * obj) {
        this->objects.erase(std::find(this->objects.begin(), this->objects.end(), obj));
    }

    friend class Camera;
};


class Camera
{
private:
    typedef std::function<void(int, int, unsigned char r, unsigned char g, unsigned char b)> SetPixel;
    SetPixel setPixel;
    int width = 1920;
    int height = 1080;
    Vector3 light_dir = {0, 1, 0};
    std::unique_ptr<float[]> zBuffer;

    void drawTriangle(Vector3 t0, Vector3 t1, Vector3 t2, Color color) {
        int ymax, ymin, ysc, e1, e, i;
        int x[3], y[3];
        float z[3];

        x[0] = t0.x;
        x[1] = t1.x;
        x[2] = t2.x;

        y[0] = t0.y;
        y[1] = t1.y;
        y[2] = t2.y;

        z[0] = t0.z;
        z[1] = t1.z;
        z[2] = t2.z;

        ymax = ymin = y[0];
        if (ymax < y[1]) ymax = y[1]; else if (ymin > y[1]) ymin = y[1];
        if (ymax < y[2]) ymax = y[2]; else if (ymin > y[2]) ymin = y[2];
        int ne;
        int x1, x2, xsc1, xsc2;
        float z1, z2, tc, z_;
        ymin = (ymin < 0) ? 0 : ymin;
        ymax = (ymax < height) ? ymax : height - 1;

        for (ysc = ymin; ysc < ymax; ysc++) {
            ne = 0;
            for (int e = 0; e < 3; e++) {
                e1 = e + 1;
                if (e1 == 3) {
                    e1 = 0;
                }
                if (y[e] < y[e1]) {
                    if (y[e1] <= ysc || ysc < y[e]) {
                        continue;
                    }
                } else {
                    if (y[e] > y[e1]) {
                        if (y[e1] > ysc || ysc >= y[e]) {
                            continue;
                        }
                    } else {
                        continue;
                    }
                }
                tc = float(y[e] - ysc) / (y[e] - y[e1]);
                if (ne) {
                    x2 = x[e] + int(tc * (x[e1] - x[e]));
                    z2 = z[e] + tc * (z[e1] - z[e]);
                } else {
                    x1 = x[e] + int(tc * (x[e1] - x[e]));
                    z1 = z[e] + tc * (z[e1] - z[e]);
                    ne = 1;
                }
            }
            if (x2 < x1) {
                e = x1;
                x1 = x2;
                x2 = e;
                tc = z1;
                z1 = z2;
                z2 = tc;
            }
            xsc1 = (x1 < 0) ? 0 : x1;
            xsc2 = (x2 < width) ? x2 : width - 1;
            for (int xsc = xsc1; xsc <= xsc2; xsc++) {
                tc = float(x1 - xsc) / (x1 - x2);
                z_ = z1 + tc * (z2 - z1);

                if (z_ < zBuffer[ysc * width + xsc]) {
                    setPixel(xsc, ysc, color.r, color.g, color.b);
                    zBuffer[ysc * width + xsc] = z_;
                }
            }
        }
    }

    void drawObject(const Object * obj) {

        Vector3 n, center;
        for (int i = 0; i < obj->mesh->vertices.size(); i++) {
            center = center + obj->mesh->vertices[i];
        }
        center = center / obj->mesh->vertices.size();

        Matrix4x4 m1 = Matrix4x4::Perspective() *
                       Matrix4x4::TRS(obj->transform.position, obj->transform.rotation, obj->transform.scale);
        Matrix4x4 m2 = Matrix4x4::TRS(obj->transform.position, obj->transform.rotation, obj->transform.scale);

        for (int i = 0; i < obj->mesh->triangles.size(); i += 3) {
            Vector3 localA, localB, localC;
            localA = obj->mesh->vertices[obj->mesh->triangles[i]];
            localB = obj->mesh->vertices[obj->mesh->triangles[i + 1]];
            localC = obj->mesh->vertices[obj->mesh->triangles[i + 2]];

            Vector4 CCVa, CCVb, CCVc;
            CCVa = Vector4(localA.x, localA.y, localA.z, 1);
            CCVb = Vector4(localB.x, localB.y, localB.z, 1);
            CCVc = Vector4(localC.x, localC.y, localC.z, 1);

            Vector4 an4, bn4, cn4;
            /*an4 = CCVa;
            bn4 = CCVb;
            cn4 = CCVc;*/

            an4 = m2 * CCVa;
            bn4 = m2 * CCVb;
            cn4 = m2 * CCVc;

            CCVa = m1 * CCVa;
            CCVb = m1 * CCVb;
            CCVc = m1 * CCVc;

            Vector3 CCV1 = {CCVa.x / CCVa.w, CCVa.y / CCVa.w, CCVa.z / CCVa.w};
            Vector3 CCV2 = {CCVb.x / CCVb.w, CCVb.y / CCVb.w, CCVa.z / CCVb.w};
            Vector3 CCV3 = {CCVc.x / CCVc.w, CCVc.y / CCVc.w, CCVa.z / CCVc.w};

            //std::cout << CCV1 << "|" << CCV2 << "|" << CCV3 << std::std::endl;

            if (CCVa.z / CCVa.w > 1 || CCVa.z / CCVa.w < -1 ||
                CCVb.z / CCVb.w > 1 || CCVb.z / CCVb.w < -1 ||
                CCVc.z / CCVc.w > 1 || CCVc.z / CCVc.w < -1) {
                continue;
            }
            n = Vector3::cross({bn4.x - an4.x, bn4.y - an4.y, bn4.z - an4.z},
                               {cn4.x - an4.x, cn4.y - an4.y, cn4.z - an4.z});
            n = n.normalized();

            float k = abs(Vector3::dot(n, light_dir));
            float depth = 1.f;

            Vector3 a1 = Vector3(width * 0.5f * (1 + CCVa.x / CCVa.w),height * 0.5f * (1 + CCVa.y / CCVa.w),
                                 CCVa.z * depth / CCVa.w);
            Vector3 b1 = Vector3(width * 0.5f * (1 + CCVb.x / CCVb.w), height * 0.5f * (1 + CCVb.y / CCVb.w),
                                 CCVb.z * depth / CCVb.w);
            Vector3 c1 = Vector3(width * 0.5f * (1 + CCVc.x / CCVc.w), height * 0.5f * (1 + CCVc.y / CCVc.w),
                                 CCVc.z * depth / CCVc.w);

            drawTriangle(a1, b1, c1, obj->color * k);
        }
    }
public:
	Camera(SetPixel setPixel, int width, int height) : zBuffer(new float[width * height]) {
	    this->setPixel = std::move(setPixel);
	    this->width = width;
	    this->height = height;
	}

	void setLightDirection(const Vector3 & v) {
        this->light_dir = v.normalized();
    }

    void Render(const Scene & scene) {
        std::fill(zBuffer.get(), zBuffer.get()+width*height, defaultDepth);
        for (auto obj : scene.objects) {
            drawObject(obj);
        }
    }

};
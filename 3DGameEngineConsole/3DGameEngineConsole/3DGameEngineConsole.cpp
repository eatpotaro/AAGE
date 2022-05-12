#include "Header/olcConsoleGameEngine.h"
#include <fstream>
#include <strstream>
#include <algorithm>
using namespace std;

struct Vec3D
{
    float x = 0;
    float y = 0;
    float z = 0;
    float w = 1;
};

struct Triangle
{
    Vec3D Points[3];

    wchar_t Sym;
    short col;
};

struct Mesh
{
    vector<Triangle> Tris;

    bool LoadFromObjectFile(string SFileName)
    {
        ifstream f(SFileName);
        if (!f.is_open())
            return false;

        vector<Vec3D> verts;

        while (!f.eof())
        {
            char line[128];
            f.getline(line, 128);

            strstream s;
            s << line;

            char junk;

            if (line[0] == 'v')
            {
                Vec3D v;
                s >> junk >> v.x >> v.y >> v.z;
                verts.push_back(v);
            }

            if (line[0] == 'f')
            {
                int f[3];
                s >> junk >> f[0] >> f[1] >> f[2];
                Tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
            }
        }

        return true;
    }
};

struct Mat4x4
{
    float M[4][4] = { 0 };
};

class AAGEEngine3D : public olcConsoleGameEngine
{
private:
    Mesh MeshCube;
    Mat4x4 ProjectionMatrix;
    float FTheta;

    Vec3D VCamera = { 0.0f, 0.0f, 0.0f };
    Vec3D VLookDir;

    float FYaw;


#pragma region Utilites
    void MultiplyMatrixVector(Vec3D& i, Vec3D& o, Mat4x4& m)
    {
        o.x = i.x * m.M[0][0] + i.y * m.M[1][0] + i.z * m.M[2][0] + m.M[3][0];
        o.y = i.x * m.M[0][1] + i.y * m.M[1][1] + i.z * m.M[2][1] + m.M[3][1];
        o.z = i.x * m.M[0][2] + i.y * m.M[1][2] + i.z * m.M[2][2] + m.M[3][2];
        float w = i.x * m.M[0][3] + i.y * m.M[1][3] + i.z * m.M[2][3] + m.M[3][3];

        if (w != 0.0f)
        {
            o.x /= w; o.y /= w; o.z /= w;
        }
    }

    Vec3D Matrix_MultiplyVector(Mat4x4& m, Vec3D& i)
    {
        Vec3D v;
        v.x = i.x * m.M[0][0] + i.y * m.M[1][0] + i.z * m.M[2][0] + i.w * m.M[3][0];
        v.y = i.x * m.M[0][1] + i.y * m.M[1][1] + i.z * m.M[2][1] + i.w * m.M[3][1];
        v.z = i.x * m.M[0][2] + i.y * m.M[1][2] + i.z * m.M[2][2] + i.w * m.M[3][2];
        v.w = i.x * m.M[0][3] + i.y * m.M[1][3] + i.z * m.M[2][3] + i.w * m.M[3][3];
        return v;
    }

    Mat4x4 Matrix_MakeIdentity()
    {
        Mat4x4 matrix;
        matrix.M[0][0] = 1.0f;
        matrix.M[1][1] = 1.0f;
        matrix.M[2][2] = 1.0f;
        matrix.M[3][3] = 1.0f;
        return matrix;
    }

    Mat4x4 Matrix_MakeRotationX(float fAngleRad)
    {
        Mat4x4 matrix;
        matrix.M[0][0] = 1.0f;
        matrix.M[1][1] = cosf(fAngleRad);
        matrix.M[1][2] = sinf(fAngleRad);
        matrix.M[2][1] = -sinf(fAngleRad);
        matrix.M[2][2] = cosf(fAngleRad);
        matrix.M[3][3] = 1.0f;
        return matrix;
    }

    Mat4x4 Matrix_MakeRotationY(float fAngleRad)
    {
        Mat4x4 matrix;
        matrix.M[0][0] = cosf(fAngleRad);
        matrix.M[0][2] = sinf(fAngleRad);
        matrix.M[2][0] = -sinf(fAngleRad);
        matrix.M[1][1] = 1.0f;
        matrix.M[2][2] = cosf(fAngleRad);
        matrix.M[3][3] = 1.0f;
        return matrix;
    }

    Mat4x4 Matrix_MakeRotationZ(float fAngleRad)
    {
        Mat4x4 matrix;
        matrix.M[0][0] = cosf(fAngleRad);
        matrix.M[0][1] = sinf(fAngleRad);
        matrix.M[1][0] = -sinf(fAngleRad);
        matrix.M[1][1] = cosf(fAngleRad);
        matrix.M[2][2] = 1.0f;
        matrix.M[3][3] = 1.0f;
        return matrix;
    }

    Mat4x4 Matrix_MakeTranslation(float x, float y, float z)
    {
        Mat4x4 matrix;
        matrix.M[0][0] = 1.0f;
        matrix.M[1][1] = 1.0f;
        matrix.M[2][2] = 1.0f;
        matrix.M[3][3] = 1.0f;
        matrix.M[3][0] = x;
        matrix.M[3][1] = y;
        matrix.M[3][2] = z;
        return matrix;
    }

    Mat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
    {
        float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
        Mat4x4 matrix;
        matrix.M[0][0] = fAspectRatio * fFovRad;
        matrix.M[1][1] = fFovRad;
        matrix.M[2][2] = fFar / (fFar - fNear);
        matrix.M[3][2] = (-fFar * fNear) / (fFar - fNear);
        matrix.M[2][3] = 1.0f;
        matrix.M[3][3] = 0.0f;
        return matrix;
    }

    Mat4x4 Matrix_MultiplyMatrix(Mat4x4& m1, Mat4x4& m2)
    {
        Mat4x4 matrix;
        for (int c = 0; c < 4; c++)
            for (int r = 0; r < 4; r++)
                matrix.M[r][c] = m1.M[r][0] * m2.M[0][c] + m1.M[r][1] * m2.M[1][c] + m1.M[r][2] * m2.M[2][c] + m1.M[r][3] * m2.M[3][c];
        return matrix;
    }

    Mat4x4 Matrix_PointAt(Vec3D& pos, Vec3D& target, Vec3D& up)
    {
        // Calculate new forward direction
        Vec3D newForward = Vector_Sub(target, pos);
        newForward = Vector_Normalise(newForward);

        // Calculate new Up direction
        Vec3D a = Vector_Mul(newForward, Vector_DotProduct(up, newForward));
        Vec3D newUp = Vector_Sub(up, a);
        newUp = Vector_Normalise(newUp);

        // New Right direction is easy, its just cross product
        Vec3D newRight = Vector_CrossProduct(newUp, newForward);

        // Construct Dimensioning and Translation Matrix	
        Mat4x4 matrix;
        matrix.M[0][0] = newRight.x;	matrix.M[0][1] = newRight.y;	matrix.M[0][2] = newRight.z;	matrix.M[0][3] = 0.0f;
        matrix.M[1][0] = newUp.x;		matrix.M[1][1] = newUp.y;		matrix.M[1][2] = newUp.z;		matrix.M[1][3] = 0.0f;
        matrix.M[2][0] = newForward.x;	matrix.M[2][1] = newForward.y;	matrix.M[2][2] = newForward.z;	matrix.M[2][3] = 0.0f;
        matrix.M[3][0] = pos.x;			matrix.M[3][1] = pos.y;			matrix.M[3][2] = pos.z;			matrix.M[3][3] = 1.0f;
        return matrix;

    }

    Mat4x4 Matrix_QuickInverse(Mat4x4& m) // Only for Rotation/Translation Matrices
    {
        Mat4x4 matrix;
        matrix.M[0][0] = m.M[0][0]; matrix.M[0][1] = m.M[1][0]; matrix.M[0][2] = m.M[2][0]; matrix.M[0][3] = 0.0f;
        matrix.M[1][0] = m.M[0][1]; matrix.M[1][1] = m.M[1][1]; matrix.M[1][2] = m.M[2][1]; matrix.M[1][3] = 0.0f;
        matrix.M[2][0] = m.M[0][2]; matrix.M[2][1] = m.M[1][2]; matrix.M[2][2] = m.M[2][2]; matrix.M[2][3] = 0.0f;
        matrix.M[3][0] = -(m.M[3][0] * matrix.M[0][0] + m.M[3][1] * matrix.M[1][0] + m.M[3][2] * matrix.M[2][0]);
        matrix.M[3][1] = -(m.M[3][0] * matrix.M[0][1] + m.M[3][1] * matrix.M[1][1] + m.M[3][2] * matrix.M[2][1]);
        matrix.M[3][2] = -(m.M[3][0] * matrix.M[0][2] + m.M[3][1] * matrix.M[1][2] + m.M[3][2] * matrix.M[2][2]);
        matrix.M[3][3] = 1.0f;
        return matrix;
    }

    Vec3D Vector_Add(Vec3D& v1, Vec3D& v2)
    {
        return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
    }

    Vec3D Vector_Sub(Vec3D& v1, Vec3D& v2)
    {
        return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
    }

    Vec3D Vector_Mul(Vec3D& v1, float k)
    {
        return { v1.x * k, v1.y * k, v1.z * k };
    }

    Vec3D Vector_Div(Vec3D& v1, float k)
    {
        return { v1.x / k, v1.y / k, v1.z / k };
    }

    float Vector_DotProduct(Vec3D& v1, Vec3D& v2)
    {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    float Vector_Length(Vec3D& v)
    {
        return sqrtf(Vector_DotProduct(v, v));
    }

    Vec3D Vector_Normalise(Vec3D& v)
    {
        float l = Vector_Length(v);
        return { v.x / l, v.y / l, v.z / l };
    }

    Vec3D Vector_CrossProduct(Vec3D& v1, Vec3D& v2)
    {
        Vec3D v;
        v.x = v1.y * v2.z - v1.z * v2.y;
        v.y = v1.z * v2.x - v1.x * v2.z;
        v.z = v1.x * v2.y - v1.y * v2.x;
        return v;
    }

    Vec3D Vector_IntersectPlane(Vec3D& plane_p, Vec3D& plane_n, Vec3D& lineStart, Vec3D& LineEnd)
    {
        plane_n = Vector_Normalise(plane_n);
        float plane_d = -Vector_DotProduct(plane_n, plane_p);
        float ad = Vector_DotProduct(lineStart, plane_n);
        float bd = Vector_DotProduct(LineEnd, plane_n);
        float t = (-plane_d - ad) / (bd - ad);
        Vec3D lineStartToEnd = Vector_Sub(LineEnd, lineStart);
        Vec3D lineToIntersect = Vector_Mul(lineStartToEnd, t);
        return Vector_Add(lineStart, lineToIntersect);
    }

    int Vector_ClipAgainstPlane(Vec3D plane_p, Vec3D plane_n, Triangle& in_tri, Triangle& out_tri1, Triangle& out_tri2)
    {
        //makle sure plane normal is normal;
        plane_n = Vector_Normalise(plane_n);

        //return signed shortest distance from point to plkant, plane normal must be normal;
        auto dist = [&](Vec3D &p)
        {
            Vec3D n = Vector_Normalise(p);
            return (plane_n.x * p.x + plane_n.y * p.y + plane_n.x * p.z - Vector_DotProduct(plane_n, plane_p));
        };

        //Create two temporaty strage arrays to classify poihtns either side of plane;
        //if distnace sign is positive, point lies on "inside" of trhe plane;
        Vec3D* inside_points[3]; int nInsidePointCount = 0;
        Vec3D* outside_points[3]; int nOutsidePointCount = 0;

        //Get sidned disnatce of each opint in trangle to plane;
        float d0 = dist(in_tri.Points[0]);
        float d1 = dist(in_tri.Points[1]);
        float d2 = dist(in_tri.Points[2]);

        if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.Points[0]; }
        else { outside_points[nOutsidePointCount++] = &in_tri.Points[0]; }
        if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.Points[1]; }
        else { outside_points[nOutsidePointCount++] = &in_tri.Points[1]; }
        if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.Points[2]; }
        else { outside_points[nOutsidePointCount++] = &in_tri.Points[2]; }

        //classify trangle points, and break the iunput triangle 
        //into smaller output triangles if required;
        //there are 4 outcomes;

        if (nInsidePointCount == 0)
        {
            //all points are outside;
            //stop living triangel1!!1;
            return 0; //no trangles;
        }

        if (nInsidePointCount == 3)
        {
            //all points in plane;
            out_tri1 = in_tri;

            return 1; //just the one original trangle is good to go;
        }

        if (nInsidePointCount == 1 && nOutsidePointCount == 2)
        {
            //trangle should be clipped. 2 points outside;
            //trangl e becomes a smaller triangle;

            //copy looks to trangle that is new;
            out_tri1.col = in_tri.col;
            out_tri1.Sym = in_tri.Sym;

            //inside points still god so we can use it;
            out_tri1.Points[0] = *inside_points[0];

            //but the two new points are at the location where the;
            //original side of rthe trangle intersect with the plane;
            out_tri1.Points[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
            out_tri1.Points[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

            return 1; //gives back the one good triangel;
        }

        if (nInsidePointCount == 2 && nOutsidePointCount == 1)
        {
            //Triangle shoudl be cliuped. two points in plane,
            //the clipped tiangle becomes S Q U A R E;
            //make S Q U A R E with 2 traignels;
            out_tri1.col = in_tri.col;
            out_tri1.Sym = in_tri.Sym;

            out_tri2.col = in_tri.col;
            out_tri2.Sym = in_tri.Sym;

            //first triangle is the 2 inside points and new point
            //chosen by the loaction where one side of teh tirangle
            //intersects with the polane;
            out_tri1.Points[0] = *inside_points[0];
            out_tri1.Points[1] = *inside_points[1];
            out_tri1.Points[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

            //secdon traignel is composed of on eof the inside points, 
            //a new point chosen byt he intersectrion of the other side of the 
            //traignel and the plane, and the newly created point above;
            out_tri2.Points[0] = *inside_points[1];
            out_tri2.Points[1] = out_tri1.Points[1];
            out_tri2.Points[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

            return 2; //return the 2 new good triagnle;
        }
    }

    // Taken From Command Line Webcam Video
    CHAR_INFO GetColour(float lum)
    {
        short bg_col, fg_col;
        wchar_t sym;
        int pixel_bw = (int)(13.0f * lum);
        switch (pixel_bw)
        {
        case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;

        case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
        case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
        case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
        case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

        case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
        case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
        case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
        case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

        case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
        case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
        case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
        case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;
        default:
            bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
        }

        CHAR_INFO c;
        c.Attributes = bg_col | fg_col;
        c.Char.UnicodeChar = sym;
        return c;
    }
#pragma endregion


public:
    AAGEEngine3D()
    {
        m_sAppName = L"3D Demo";
    }

public:
    bool OnUserCreate() override
    {
        MeshCube.LoadFromObjectFile("Resources\\Crate.obj");

        //Projetion Matrix
        float FNear = 0.1f;
        float FFar = 1000.0f;
        float FFov = 90.0f;
        float FAspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
        float FFovRad = 1.0f / tanf(FFov * 0.5f / 180.0 * 3.13159f);

        ProjectionMatrix.M[0][0] = FAspectRatio * FFovRad;
        ProjectionMatrix.M[1][1] = FFovRad;
        ProjectionMatrix.M[2][2] = FFar / (FFar - FNear);
        ProjectionMatrix.M[3][2] = (-FFar * FNear) / (FFar - FNear);
        ProjectionMatrix.M[2][3] = 1.0f;
        ProjectionMatrix.M[3][3] = 0.0f;

        return true;
    }

    bool OnUserUpdate(float fElapsedTime) override
    {
        if (GetKey(VK_UP).bHeld)
            VCamera.y += 8.0f * fElapsedTime;
        if (GetKey(VK_DOWN).bHeld)
            VCamera.y -= 8.0f * fElapsedTime;
        if (GetKey(VK_LEFT).bHeld)
            VCamera.x += 8.0f * fElapsedTime;
        if (GetKey(VK_RIGHT).bHeld)
            VCamera.x -= 8.0f * fElapsedTime;

        Vec3D VForward = Vector_Mul(VLookDir, 8.0f * fElapsedTime);

        if (GetKey(L'W').bHeld)
            VCamera = Vector_Add(VCamera, VForward);
        if (GetKey(L'S').bHeld)
            VCamera = Vector_Sub(VCamera, VForward);

        if (GetKey(L'A').bHeld)
            FYaw -= 2.0f * fElapsedTime;
        if (GetKey(L'D').bHeld)
            FYaw += 2.0f * fElapsedTime;

        //clear screen;
        Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

        Mat4x4 MatRotx, MatRotz;
        //FTheta += 1.0f * fElapsedTime;

        //RotationMatriXZ;
        MatRotz = Matrix_MakeRotationZ(FTheta);
        MatRotx = Matrix_MakeRotationX(FTheta * 0.5);

        Mat4x4 MatTrans;
        MatTrans = Matrix_MakeTranslation(0.0f, 0.0f, 16.0f);

        Mat4x4 MatWorld;
        MatWorld = Matrix_MakeIdentity();
        MatWorld = Matrix_MultiplyMatrix(MatRotz, MatRotx);
        MatWorld = Matrix_MultiplyMatrix(MatWorld, MatTrans);

        Vec3D VUp = { 0.0f, 1.0f, 0.0f };
        Vec3D VTarget = { 0, 0, 1 };
        Mat4x4 MatCameraRot = Matrix_MakeRotationY(FYaw);
        VLookDir = Matrix_MultiplyVector(MatCameraRot, VTarget);
        VTarget = Vector_Add(VCamera, VLookDir);

        Mat4x4 MatCamera = Matrix_PointAt(VCamera, VTarget, VUp);
        Mat4x4 MatView = Matrix_QuickInverse(MatCamera);

        vector<Triangle> VecTrianglesToRaster;

        //draw triangle;
        for (auto Tri : MeshCube.Tris)
        {
            Triangle TriProjected, TriTransformed, TriViewed;

            //Transorm triangle;
            TriTransformed.Points[0] = Matrix_MultiplyVector(MatWorld, Tri.Points[0]);
            TriTransformed.Points[1] = Matrix_MultiplyVector(MatWorld, Tri.Points[1]);
            TriTransformed.Points[2] = Matrix_MultiplyVector(MatWorld, Tri.Points[2]);

            //calculate normal of rtriangle;
            Vec3D Normal, Line1, Line2;
            Line1 = Vector_Sub(TriTransformed.Points[1], TriTransformed.Points[0]);
            Line2 = Vector_Sub(TriTransformed.Points[2], TriTransformed.Points[0]);

            //Take Cross product of the linbes to get the transgle's suface normal;
            Normal = Vector_CrossProduct(Line1, Line2);

            //Normalise;
            Normal = Vector_Normalise(Normal);

            /*
            float l = sqrtf((Normal.x * Normal.x) + (Normal.y * Normal.y) + (Normal.z * Normal.z));
            Normal.x /= l; Normal.y /= l; Normal.z /= l;
            */

            Vec3D VCameraRay = Vector_Sub(TriTransformed.Points[0], VCamera);

            if (Vector_DotProduct(Normal, VCameraRay) < 0.0f) //buggy prevension of culling (probably shoiudlnt do); 
            {
                //Illumintaionj;
                Vec3D LightDirection = { 0.0f, 1.0f, -1.0f };
                LightDirection = Vector_Normalise(LightDirection);

                float DP = max(0.1f, Vector_DotProduct(LightDirection, Normal));

                CHAR_INFO c = GetColour(DP);
                TriTransformed.col = c.Attributes;
                TriTransformed.Sym = c.Char.UnicodeChar;

                //world space to view space;
                TriViewed.Points[0] = Matrix_MultiplyVector(MatView, TriTransformed.Points[0]);
                TriViewed.Points[1] = Matrix_MultiplyVector(MatView, TriTransformed.Points[1]);
                TriViewed.Points[2] = Matrix_MultiplyVector(MatView, TriTransformed.Points[2]);

                //clip viewed tanrgle against near plane, this could make two mroe traingles;
                /*int nClippedTriangles = 0;
                Triangle clipped[3];
                nClippedTriangles = Vector_ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, TriViewed, clipped[0], clipped[1]);

                for (int n = 0; n < nClippedTriangles; n++)
                {*/
                //3D to 2D;
                TriProjected.Points[0] = Matrix_MultiplyVector(ProjectionMatrix, TriViewed.Points[0]);
                TriProjected.Points[1] = Matrix_MultiplyVector(ProjectionMatrix, TriViewed.Points[1]);
                TriProjected.Points[2] = Matrix_MultiplyVector(ProjectionMatrix, TriViewed.Points[2]);
                TriProjected.col = TriTransformed.col;
                TriProjected.Sym = c.Char.UnicodeChar;

                TriProjected.Points[0] = Vector_Div(TriProjected.Points[0], TriProjected.Points[0].w);
                TriProjected.Points[1] = Vector_Div(TriProjected.Points[1], TriProjected.Points[1].w);
                TriProjected.Points[2] = Vector_Div(TriProjected.Points[2], TriProjected.Points[2].w);

                TriProjected.Points[0].x *= -1;
                TriProjected.Points[0].y *= -1;
                TriProjected.Points[1].x *= -1;
                TriProjected.Points[1].y *= -1;
                TriProjected.Points[2].x *= -1;
                TriProjected.Points[2].y *= -1;

                //scale into view;
                Vec3D VOffsetView = { 1, 1, 0 };
                TriProjected.Points[0] = Vector_Add(TriProjected.Points[0], VOffsetView);
                TriProjected.Points[1] = Vector_Add(TriProjected.Points[1], VOffsetView);
                TriProjected.Points[2] = Vector_Add(TriProjected.Points[2], VOffsetView);

                TriProjected.Points[0].x *= 0.5 * (float)ScreenWidth();
                TriProjected.Points[0].y *= 0.5 * (float)ScreenHeight();
                TriProjected.Points[1].x *= 0.5 * (float)ScreenWidth();
                TriProjected.Points[1].y *= 0.5 * (float)ScreenHeight();
                TriProjected.Points[2].x *= 0.5 * (float)ScreenWidth();
                TriProjected.Points[2].y *= 0.5 * (float)ScreenHeight();

                VecTrianglesToRaster.push_back(TriProjected);
                //}
            }

        }

        //sort back to front;

        sort(VecTrianglesToRaster.begin(), VecTrianglesToRaster.end(), [](Triangle& t1, Triangle& t2)
        {
            float z1 = (t1.Points[0].z + t1.Points[1].z + t1.Points[2].z) / 3.0f;
            float z2 = (t2.Points[0].z + t2.Points[1].z + t2.Points[2].z) / 3.0f;

            return z1 > z2;
        });

        for (auto& TriProjected : VecTrianglesToRaster)
        {
            //rasterize triangle
            FillTriangle(TriProjected.Points[0].x, TriProjected.Points[0].y,
                TriProjected.Points[1].x, TriProjected.Points[1].y,
                TriProjected.Points[2].x, TriProjected.Points[2].y,
                TriProjected.Sym, TriProjected.col);


            DrawTriangle(TriProjected.Points[0].x, TriProjected.Points[0].y,
                TriProjected.Points[1].x, TriProjected.Points[1].y,
                TriProjected.Points[2].x, TriProjected.Points[2].y,
                PIXEL_SOLID, FG_GREEN);
        }

        return true;
    }

};


int main()
{
    AAGEEngine3D demo;
    if (demo.ConstructConsole(450, 225, 4, 4))
        demo.Start();
    return 0;
}
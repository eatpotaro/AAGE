#include "Header/olcConsoleGameEngine.h"
#include <fstream>
#include <strstream>
#include <algorithm>
using namespace std;

struct Vec3D
{
    float x, y, z;
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

    Vec3D vCamera = { 0.0f, 0.0f, 0.0f };

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


public: 
    AAGEEngine3D()
    {
        m_sAppName = L"3D Demo"; 
    }

public:
    bool OnUserCreate() override
    {
        /*MeshCube.Tris = {
            //south;
            {0.0f, 0.0f, 0.0f,  0.0f, 1.0f, 0.0f,   1.0f, 1.0f, 0.0f},
            {0.0f, 0.0f, 0.0f,  1.0f, 1.0f, 0.0f,   1.0f, 0.0f, 0.0f},

            //east;
            {1.0f, 0.0f, 0.0f,  1.0f, 1.0f, 0.0f,   1.0f, 1.0f, 1.0f},
            {1.0f, 0.0f, 0.0f,  1.0f, 1.0f, 1.0f,   1.0f, 0.0f, 1.0f},

            //north;
            {1.0f, 0.0f, 1.0f,  1.0f, 1.0f, 1.0f,   0.0f, 1.0f, 1.0f},
            {1.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 0.0f, 1.0f},

            //west;
            {0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 1.0f,   0.0f, 1.0f, 0.0f},
            {0.0f, 0.0f, 1.0f,  0.0f, 1.0f, 0.0f,   0.0f, 0.0f, 0.0f},

            //top;
            {0.0f, 1.0f, 0.0f,  0.0f, 1.0f, 1.0f,   1.0f, 1.0f, 1.0f},
            {0.0f, 1.0f, 0.0f,  1.0f, 1.0f, 1.0f,   1.0f, 1.0f, 0.0f},

            //bottom;
            {1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 1.0f,   0.0f, 0.0f, 0.0f},
            {1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 0.0f,   1.0f, 0.0f, 0.0f},
        };*/
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
        //clear screen;
        Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);
        
        Mat4x4 MatRotx, MatRotz;
        FTheta += 1.0f * fElapsedTime;

        //RotationMatrixZ;
        MatRotz.M[0][0] = cosf(FTheta);
        MatRotz.M[0][1] = sinf(FTheta);
        MatRotz.M[1][0] = -sinf(FTheta);
        MatRotz.M[1][1] = cosf(FTheta);
        MatRotz.M[2][2] = 1;
        MatRotz.M[3][3] = 1;

        //RotationMatrixX;
        MatRotx.M[0][0] = 1;
        MatRotx.M[1][1] = cosf(FTheta * 0.5f);
        MatRotx.M[1][2] = sinf(FTheta * 0.5f);
        MatRotx.M[2][1] = -sinf(FTheta * 0.5f);
        MatRotx.M[2][2] = cosf(FTheta * 0.5f);
        MatRotx.M[3][3] = 1;

        vector<Triangle> VecTrianglesToRaster;

        //draw triangle;
        for (auto Tri : MeshCube.Tris)
        {
            Triangle TriProjected, TriTranslated, TriRotatedZ, TriRotatedZX;

            //rotate z axis;
            MultiplyMatrixVector(Tri.Points[0], TriRotatedZ.Points[0], MatRotz);
            MultiplyMatrixVector(Tri.Points[1], TriRotatedZ.Points[1], MatRotz);
            MultiplyMatrixVector(Tri.Points[2], TriRotatedZ.Points[2], MatRotz);

            //rotate x axis;
            MultiplyMatrixVector(TriRotatedZ.Points[0], TriRotatedZX.Points[0], MatRotx);
            MultiplyMatrixVector(TriRotatedZ.Points[1], TriRotatedZX.Points[1], MatRotx);
            MultiplyMatrixVector(TriRotatedZ.Points[2], TriRotatedZX.Points[2], MatRotx);

            //Translate into screen;
            TriTranslated = TriRotatedZX;
            TriTranslated.Points[0].z = TriRotatedZX.Points[0].z + 3.0f;
            TriTranslated.Points[1].z = TriRotatedZX.Points[1].z + 3.0f;
            TriTranslated.Points[2].z = TriRotatedZX.Points[2].z + 3.0f;

            Vec3D normal, line1, line2;
            line1.x = TriTranslated.Points[1].x - TriTranslated.Points[0].x;
            line1.y = TriTranslated.Points[1].y - TriTranslated.Points[0].y;
            line1.z = TriTranslated.Points[1].z - TriTranslated.Points[0].z;

            line2.x = TriTranslated.Points[2].x - TriTranslated.Points[0].x;
            line2.y = TriTranslated.Points[2].y - TriTranslated.Points[0].y;
            line2.z = TriTranslated.Points[2].z - TriTranslated.Points[0].z;

            normal.x = line1.y * line2.z - line1.z * line2.y;
            normal.y = line1.z * line2.x - line1.x * line2.z;
            normal.z = line1.x * line2.y - line1.y * line2.x;

            float l = sqrtf((normal.x * normal.x) + (normal.y * normal.y) + (normal.z * normal.z));
            normal.x /= l; normal.y /= l; normal.z /= l;

            if ((normal.x * (TriTranslated.Points[0].x - vCamera.x) + 
                normal.y * (TriTranslated.Points[0].y - vCamera.y) + 
                normal.z * (TriTranslated.Points[0].z - vCamera.z)) < 0.0f)
            {
                //Illumination;
                Vec3D LightDirection = { 0.0f, 0.0f, -1.0f };
                float Lld = sqrtf(LightDirection.x * LightDirection.x +
                    LightDirection.y * LightDirection.y +
                    LightDirection.z * LightDirection.z);
                LightDirection.x /= Lld; LightDirection.y /= Lld; LightDirection.z /= Lld;

                float DP = normal.x * LightDirection.x + normal.y * LightDirection.y + normal.z * LightDirection.z;

                CHAR_INFO c = GetColour(DP);
                TriTranslated.col = c.Attributes;
                TriTranslated.Sym = c.Char.UnicodeChar;


                //3D to 2D;
                MultiplyMatrixVector(TriTranslated.Points[0], TriProjected.Points[0], ProjectionMatrix);
                MultiplyMatrixVector(TriTranslated.Points[1], TriProjected.Points[1], ProjectionMatrix);
                MultiplyMatrixVector(TriTranslated.Points[2], TriProjected.Points[2], ProjectionMatrix);
                TriProjected.col = TriTranslated.col;
                TriProjected.Sym = TriTranslated.Sym;

                //scale into view;
                TriProjected.Points[0].x += 1.0f;   TriProjected.Points[0].y += 1.0f;
                TriProjected.Points[1].x += 1.0f;   TriProjected.Points[1].y += 1.0f;
                TriProjected.Points[2].x += 1.0f;   TriProjected.Points[2].y += 1.0f;

                TriProjected.Points[0].x *= 0.5 * (float)ScreenWidth();
                TriProjected.Points[0].y *= 0.5 * (float)ScreenHeight();
                TriProjected.Points[1].x *= 0.5 * (float)ScreenWidth();
                TriProjected.Points[1].y *= 0.5 * (float)ScreenHeight();
                TriProjected.Points[2].x *= 0.5 * (float)ScreenWidth();
                TriProjected.Points[2].y *= 0.5 * (float)ScreenHeight();
                 
                VecTrianglesToRaster.push_back(TriProjected);                
            }

        }

        //sort back to front;

        sort(VecTrianglesToRaster.begin(), VecTrianglesToRaster.end(), [](Triangle &t1, Triangle &t2)
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


            /*DrawTriangle(TriProjected.Points[0].x, TriProjected.Points[0].y,
                TriProjected.Points[1].x, TriProjected.Points[1].y,
                TriProjected.Points[2].x, TriProjected.Points[2].y,
                PIXEL_SOLID, FG_GREEN);*/
        }

        return true;
    }

};


int main()
{
    AAGEEngine3D demo;
    if (demo.ConstructConsole(250, 250, 4, 4))
        demo.Start();
    return 0;
}

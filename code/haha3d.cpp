#include "haha3d.h"
#include "haha3d_math.cpp"
#include "haha3d_render_command_buffer.cpp"
#include <vector>

mat4 
camera::GetRotationMatrix(void)
{ 
    mat4 Result = LookAt(P, P + Dir); 

    return(Result);
} 

internal b32
Contains(vec3 *Vertices, u32 Count, vec3 Check)
{
    b32 Result = false;

    for(u32 Index = 0;
        Index < Count;
        Index++)
    {
        vec3 Vertex = Vertices[Index];
        b32 XIsEqual = Absolute(Vertex.x - Check.x) <= Epsilon;
        b32 YIsEqual = Absolute(Vertex.y - Check.y) <= Epsilon;
        b32 ZIsEqual = Absolute(Vertex.z - Check.z) <= Epsilon;
        if(XIsEqual && YIsEqual && ZIsEqual)
        {
            Result = true;
            break;
        }
    }

    return(Result);
}

internal b32
VectorsAreEqual(vec3 A, vec3 B)
{
    b32 XIsEqual = Absolute(A.x - B.x) <= Epsilon;
    b32 YIsEqual = Absolute(A.y - B.y) <= Epsilon;
    b32 ZIsEqual = Absolute(A.z - B.z) <= Epsilon;

    b32 Result = XIsEqual && YIsEqual && ZIsEqual;
    return(Result);
}

internal vec3
Support(game_object *Obj, vec3 D)
{
    vec3 Result;
    // TODO(georgy): It is not the best idea to do all these computations
    vec3 AxisX = vec3(Obj->RigidBody.Orientation.a11, Obj->RigidBody.Orientation.a21, Obj->RigidBody.Orientation.a31);
    vec3 AxisY = vec3(Obj->RigidBody.Orientation.a12, Obj->RigidBody.Orientation.a22, Obj->RigidBody.Orientation.a32);
    vec3 AxisZ = vec3(Obj->RigidBody.Orientation.a13, Obj->RigidBody.Orientation.a23, Obj->RigidBody.Orientation.a33);

    vec3 Points[8] = 
    {
        Obj->RigidBody.P + 0.5f*Obj->Width*AxisX + 0.5f*Obj->Height*AxisY + 0.5f*Obj->Depth*AxisZ,
        Obj->RigidBody.P + 0.5f*Obj->Width*AxisX + 0.5f*Obj->Height*AxisY - 0.5f*Obj->Depth*AxisZ,
        Obj->RigidBody.P - 0.5f*Obj->Width*AxisX + 0.5f*Obj->Height*AxisY - 0.5f*Obj->Depth*AxisZ,
        Obj->RigidBody.P - 0.5f*Obj->Width*AxisX + 0.5f*Obj->Height*AxisY + 0.5f*Obj->Depth*AxisZ,
        Obj->RigidBody.P - 0.5f*Obj->Width*AxisX - 0.5f*Obj->Height*AxisY + 0.5f*Obj->Depth*AxisZ,
        Obj->RigidBody.P + 0.5f*Obj->Width*AxisX - 0.5f*Obj->Height*AxisY + 0.5f*Obj->Depth*AxisZ,
        Obj->RigidBody.P + 0.5f*Obj->Width*AxisX - 0.5f*Obj->Height*AxisY - 0.5f*Obj->Depth*AxisZ,
        Obj->RigidBody.P - 0.5f*Obj->Width*AxisX - 0.5f*Obj->Height*AxisY - 0.5f*Obj->Depth*AxisZ,
    };
    u32 Count = ArrayCount(Points);

    Result = Points[0];
    r32 Max = Dot(D, Points[0]);
    for(u32 Point = 1;  
        Point < Count;
        Point++)
    {
        r32 Value = Dot(D, Points[Point]);
        if(Value > Max)
        {
            Max = Value;
            Result = Points[Point];
        }
    }

    return(Result);
}

internal void
DoSimplexTetrahedronSubset(closest_voronoi_region_triangle_point TetrahedronVR, 
                           vec3 A, vec3 B, vec3 C, vec3 *TriangleFromSimplexA, vec3 *TriangleFromSimplexB,
                           u32 *SimplexCount, vec3 *Simplex, vec3 *Dir, vec3 *SimplexA, vec3 *SimplexB)
{
    switch(TetrahedronVR.VR)
    {
        case VoronoiRegion_Vertex:
        {
            *SimplexCount = 1;
            Simplex[0] = TetrahedronVR.P; 
            SimplexA[0] = TriangleFromSimplexA[TetrahedronVR.VertIndex0];
            SimplexB[0] = TriangleFromSimplexB[TetrahedronVR.VertIndex0];
            *Dir = -TetrahedronVR.P;
        } break;

        case VoronoiRegion_Edge:
        {
            *SimplexCount = 2;
            Simplex[0] = TetrahedronVR.EdgeP0;
            Simplex[1] = TetrahedronVR.EdgeP1;
            SimplexA[0] = TriangleFromSimplexA[TetrahedronVR.VertIndex0];
            SimplexA[1] = TriangleFromSimplexA[TetrahedronVR.VertIndex1];
            SimplexB[0] = TriangleFromSimplexB[TetrahedronVR.VertIndex0];
            SimplexB[1] = TriangleFromSimplexB[TetrahedronVR.VertIndex1];
            *Dir = -TetrahedronVR.P;
        } break;

        case VoronoiRegion_Triangle:
        {
            *SimplexCount = 3;
            Simplex[0] = C;
            Simplex[1] = B;
            Simplex[2] = A;
            SimplexA[0] = TriangleFromSimplexA[2];
            SimplexA[1] = TriangleFromSimplexA[1];
            SimplexA[2] = TriangleFromSimplexA[0];
            SimplexB[0] = TriangleFromSimplexB[2];
            SimplexB[1] = TriangleFromSimplexB[1];
            SimplexB[2] = TriangleFromSimplexB[0];
            *Dir = -TetrahedronVR.P;
        } break;
    }
}

internal b32
DoSimplex(u32 *SimplexCount, vec3 *Simplex, vec3 *Dir, vec3 *SimplexA, vec3 *SimplexB)
{
    b32 Result = false;

    if(*SimplexCount == 2)
    {
        // NOTE(georgy): Line
        vec3 B = Simplex[0];
        vec3 A = Simplex[1];

        vec3 AB = B - A;
        vec3 AO = -A;

        if(Dot(AB, AO) > 0.0f)
        {
            *Dir = Cross(Cross(AB, AO), AB);
        }
        else
        {
            *SimplexCount = 1;
            Simplex[0] = A;
            SimplexA[0] = SimplexA[1];
            SimplexB[0] = SimplexB[1];
            *Dir = AO;
        }
    }
    else if(*SimplexCount == 3)
    {
        // NOTE(georgy): Triangle
        vec3 C = Simplex[0];
        vec3 B = Simplex[1]; 
        vec3 A = Simplex[2];
        vec3 O = vec3(0.0f, 0.0f, 0.0f);

        closest_voronoi_region_triangle_point TriangleVR = ClosestPointInTriangleVR(O, A, B, C);
        switch(TriangleVR.VR)
        {
            case VoronoiRegion_Vertex:
            {
                *SimplexCount = 1;
                Simplex[0] = TriangleVR.P;
                SimplexA[0] = SimplexA[2 - TriangleVR.VertIndex0];
                SimplexB[0] = SimplexB[2 - TriangleVR.VertIndex0];
                *Dir = -TriangleVR.P;
            } break;

            case VoronoiRegion_Edge:
            {
                *SimplexCount = 2;
                Simplex[0] = TriangleVR.EdgeP0;
                Simplex[1] = TriangleVR.EdgeP1;
                SimplexA[0] = SimplexA[2 - TriangleVR.VertIndex0];
                SimplexA[1] = SimplexA[2 - TriangleVR.VertIndex1];
                SimplexB[0] = SimplexB[2 - TriangleVR.VertIndex0];
                SimplexB[1] = SimplexB[2 - TriangleVR.VertIndex1];
                *Dir = -TriangleVR.P;
            } break;

            case VoronoiRegion_Triangle:
            {
                *Dir = -TriangleVR.P;
            } break;
        }
    }
    else if(*SimplexCount == 4)
    {
        // NOTE(georgy): Tetrahedron
        vec3 D = Simplex[0];
        vec3 C = Simplex[1];
        vec3 B = Simplex[2];
        vec3 A = Simplex[3];

        vec3 AB = B - A;
        vec3 AC = C - A;
        vec3 AD = D - A;
        vec3 AO = -A;
        vec3 O = vec3(0.0f, 0.0f, 0.0f);

        vec3 ABCNormal = Cross(AB, AC);
        vec3 ACDNormal = Cross(AC, AD);
        vec3 ADBNormal = Cross(AD, AB);
        vec3 CBDNormal = Cross(B - C, D - C);

        b32 InsideABCPlane = (Dot(ABCNormal, AO) * Dot(ABCNormal, AD) >= 0.0f);
        b32 InsideACDPlane = (Dot(ACDNormal, AO) * Dot(ACDNormal, AB) >= 0.0f);
        b32 InsideADBPlane = (Dot(ADBNormal, AO) * Dot(ADBNormal, AC) >= 0.0f);
        b32 InsideCBDPlane = (Dot(CBDNormal, O - C) * Dot(CBDNormal, A - C) >= 0.0f);
        if(InsideABCPlane && InsideACDPlane && InsideADBPlane && InsideCBDPlane) return(true);

        // NOTE(georgy): 0 - ABC, 1 - ACD, 2 - ADB, 3 - CBD
        u32 TriangleIndex = 0;
        r32 ClosestDistSq = FLT_MAX;

        closest_voronoi_region_triangle_point ABC_VR, ACD_VR, ADB_VR, CBD_VR;
        if(!InsideABCPlane)
        {
            ABC_VR = ClosestPointInTriangleVR(O, A, B, C);
            r32 DistSq = LengthSq(ABC_VR.P);
            if(DistSq < ClosestDistSq)
            {
                ClosestDistSq = DistSq;
                TriangleIndex = 0;
            }
        }
        
        if(!InsideACDPlane)
        {
            ACD_VR = ClosestPointInTriangleVR(O, A, C, D);
            r32 DistSq = LengthSq(ACD_VR.P);
            if(DistSq < ClosestDistSq)
            {
                ClosestDistSq = DistSq;
                TriangleIndex = 1;
            }
        }

        if(!InsideADBPlane)
        {
            ADB_VR = ClosestPointInTriangleVR(O, A, D, B);
            r32 DistSq = LengthSq(ADB_VR.P);
            if(DistSq < ClosestDistSq)
            {
                ClosestDistSq = DistSq;
                TriangleIndex = 2;
            }
        }

        if(!InsideCBDPlane)
        {
            CBD_VR = ClosestPointInTriangleVR(O, C, B, D);
            r32 DistSq = LengthSq(CBD_VR.P);
            if(DistSq < ClosestDistSq)
            {
                ClosestDistSq = DistSq;
                TriangleIndex = 3;
            }
        }

        switch(TriangleIndex)
        {
            case 0:
            {
                vec3 TriangleFromSimplexA[3] = { SimplexA[3], SimplexA[2], SimplexA[1] }; 
                vec3 TriangleFromSimplexB[3] = { SimplexB[3], SimplexB[2], SimplexB[1] }; 
                DoSimplexTetrahedronSubset(ABC_VR, A, B, C, TriangleFromSimplexA, TriangleFromSimplexB, 
                                           SimplexCount, Simplex, Dir, SimplexA, SimplexB);
            } break;

            case 1:
            {
                vec3 TriangleFromSimplexA[3] = { SimplexA[3], SimplexA[1], SimplexA[0] }; 
                vec3 TriangleFromSimplexB[3] = { SimplexB[3], SimplexB[1], SimplexB[0] }; 
                DoSimplexTetrahedronSubset(ACD_VR, A, C, D, TriangleFromSimplexA, TriangleFromSimplexB, 
                                          SimplexCount, Simplex, Dir, SimplexA, SimplexB);
            } break;

            case 2:
            {
                vec3 TriangleFromSimplexA[3] = { SimplexA[3], SimplexA[0], SimplexA[2] }; 
                vec3 TriangleFromSimplexB[3] = { SimplexB[3], SimplexB[0], SimplexB[2] }; 
                DoSimplexTetrahedronSubset(ADB_VR, A, D, B, TriangleFromSimplexA, TriangleFromSimplexB, 
                                           SimplexCount, Simplex, Dir, SimplexA, SimplexB);
            } break;
            
            case 3:
            {
                vec3 TriangleFromSimplexA[3] = { SimplexA[1], SimplexA[2], SimplexA[0] }; 
                vec3 TriangleFromSimplexB[3] = { SimplexB[1], SimplexB[2], SimplexB[0] }; 
                DoSimplexTetrahedronSubset(CBD_VR, C, B, D, TriangleFromSimplexA, TriangleFromSimplexB, 
                                           SimplexCount, Simplex, Dir, SimplexA, SimplexB);
            } break;
        }
    }

    Assert(*SimplexCount != 4);

    if(LengthSq(*Dir) <= Square(0.01f))
    {
        return(true);
    }

    return(Result);
}

global_variable r32 PersistentThreshold = 0.1f;
global_variable r32 PersistentThresholdSq = Square(PersistentThreshold);
internal b32
Intersect(collision_manifold *Manifold, game_object *ObjA, game_object *ObjB)
{
    vec3 Sa = Support(ObjB, vec3(1.0f, 0.0f, 0.0f));
    vec3 Sb = Support(ObjA, -vec3(1.0f, 0.0f, 0.0f));
    vec3 S = Sa - Sb;
    u32 SimplexCount = 1;
    vec3 Simplex[4] = {S, vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0)};
    vec3 SimplexA[4] = {Sa, vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0)};
    vec3 SimplexB[4] = {Sb, vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0)};
    vec3 Dir = -S;

    // NOTE(georgy): GJK pass
    b32 Result;
    while(true)
    {
        vec3 Aa = Support(ObjB, Dir);
        vec3 Ab = Support(ObjA, -Dir);
        vec3 A = Aa - Ab;
        r32 DirDotA = Dot(Dir, A);
        // if(DirDotA < 0.0f)
        r32 MaxSqMagnitude = -FLT_MAX;
        for(u32 I = 0; I < SimplexCount; I++)
        {
            if(LengthSq(Simplex[I]) > MaxSqMagnitude)
            {
                MaxSqMagnitude = LengthSq(Simplex[I]);
            }
        }
        if(Dot(A + Dir, Dir) <= 2.0f*Epsilon * MaxSqMagnitude)
        {
            Result = false;
            break;
        } 
        if(Contains(Simplex, SimplexCount, A))
        {
            Result = false;
            break;
        }
        SimplexA[SimplexCount] = Aa;
        SimplexB[SimplexCount] = Ab;
        Simplex[SimplexCount++] = A;
        Assert(VectorsAreEqual(Simplex[SimplexCount-1], SimplexA[SimplexCount-1] - SimplexB[SimplexCount-1]));
        if(DoSimplex(&SimplexCount, Simplex, &Dir, SimplexA, SimplexB))
        {
            Result = true;
            break;
        } 
    }

    // NOTE(georgy): EPA pass
    if(Result)
    {
        vec3 O = vec3(0.0f, 0.0f, 0.0f);
        r32 Tolerance = 0.00001f;
        r32 ToleranceSq = Square(Tolerance);

        switch(SimplexCount)
        {
            // TODO(georgy): If our first point is origin, GJK returns true immediately
            case 1:
            {
                vec3 SearchDirections[] = 
                {
                    vec3(1.0f, 0.0f, 0.0f),
                    vec3(-1.0f, 0.0f, 0.0f),
                    vec3(0.0f, 1.0f, 0.0f),
                    vec3(0.0f, -1.0f, 0.0f),
                    vec3(0.0f, 0.0f, 1.0f),
                    vec3(0.0f, 0.0f, -1.0f),
                };

                for(vec3 SearchDir : SearchDirections)
                {
                    SimplexA[1] = Support(ObjB, SearchDir);
                    SimplexB[1] = Support(ObjA, -SearchDir);
                    Simplex[1] = SimplexA[1] - SimplexB[1];

                    if(LengthSq(Simplex[1] - Simplex[0]) >= ToleranceSq)
                    {
                        SimplexCount++;
                        Assert(VectorsAreEqual(Simplex[SimplexCount-1], SimplexA[SimplexCount-1] - SimplexB[SimplexCount-1]));
                        break;
                    }
                }

                Assert(SimplexCount == 2);
            }

            case 2:
            {
                vec3 Axes[] = 
                {
                    vec3(1.0f, 0.0f, 0.0f),
                    vec3(0.0f, 1.0f, 0.0f),
                    vec3(0.0f, 0.0f, 1.0f),
                };

                vec3 LineVec = Simplex[1] - Simplex[0];

                u32 LeastSignificantAxis = 0;
                if(Absolute(LineVec.x) <= Absolute(LineVec.y))
                {
                    if(Absolute(LineVec.x) <= Absolute(LineVec.z))
                    {
                        LeastSignificantAxis = 0;
                    }
                    else
                    {
                        LeastSignificantAxis = 2;
                    }
                }
                else if(Absolute(LineVec.y) <= Absolute(LineVec.z))
                {
                    LeastSignificantAxis = 1;
                }
                else
                {
                    LeastSignificantAxis = 2;
                }

                vec3 SearchDir = Cross(LineVec, Axes[LeastSignificantAxis]);

                mat3 Rot = Rotation3x3(60.0f, LineVec);

                for(u32 Search = 0;
                    Search < 6;
                    Search++)
                {
                    SimplexA[2] = Support(ObjB, SearchDir);
                    SimplexB[2] = Support(ObjA, -SearchDir);
                    Simplex[2] = SimplexA[2] - SimplexB[2];

                    if(Dot(Simplex[2], SearchDir) >= ToleranceSq)
                    {
                        SimplexCount++;
                        Assert(VectorsAreEqual(Simplex[SimplexCount-1], SimplexA[SimplexCount-1] - SimplexB[SimplexCount-1]));
                        break;
                    }

                    SearchDir = Rot * SearchDir;
                }

                Assert(SimplexCount == 3);
            } 

            case 3:
            {
                vec3 SearchDir = Cross(Simplex[1] - Simplex[0], Simplex[2] - Simplex[0]);
                SimplexA[3] = Support(ObjB, SearchDir);
                SimplexB[3] = Support(ObjA, -SearchDir);
                Simplex[3] = SimplexA[3] - SimplexB[3];
                
                if(Dot(Simplex[3], SearchDir) < ToleranceSq)
                {
                    SearchDir = -SearchDir;
                    SimplexA[3] = Support(ObjB, SearchDir);
                    SimplexB[3] = Support(ObjA, -SearchDir);
                    Simplex[3] = SimplexA[3] - SimplexB[3];
                }

                SimplexCount++;
                Assert(VectorsAreEqual(Simplex[SimplexCount-1], SimplexA[SimplexCount-1] - SimplexB[SimplexCount-1]));
                Assert(SimplexCount == 4);
            }
        }

        Assert(SimplexCount == 4);

        vec3 D = Simplex[0];
        vec3 C = Simplex[1];
        vec3 B = Simplex[2];
        vec3 A = Simplex[3];

        // NOTE(georgy): ABC, ACD, ADB, CBD have the same winding order. I want them to be CCW.
        vec3 ABCNormal = Cross(B - A, C - A);
        if(Dot(ABCNormal, D - A) > 0.0f)
        {
            vec3 Temp = Simplex[1];
            Simplex[1] = Simplex[2];
            Simplex[2] = Temp;

            Temp = SimplexA[1];
            SimplexA[1] = SimplexA[2];
            SimplexA[2] = Temp;

            Temp = SimplexB[1];
            SimplexB[1] = SimplexB[2];
            SimplexB[2] = Temp;

            C = Simplex[1];
            B = Simplex[2];
                
            Assert(VectorsAreEqual(Simplex[0], SimplexA[0] - SimplexB[0]));
            Assert(VectorsAreEqual(Simplex[1], SimplexA[1] - SimplexB[1]));
            Assert(VectorsAreEqual(Simplex[2], SimplexA[2] - SimplexB[2]));
            Assert(VectorsAreEqual(Simplex[3], SimplexA[3] - SimplexB[3]));
        }

        Assert(VectorsAreEqual(Simplex[0], SimplexA[0] - SimplexB[0]));
        Assert(VectorsAreEqual(Simplex[1], SimplexA[1] - SimplexB[1]));
        Assert(VectorsAreEqual(Simplex[2], SimplexA[2] - SimplexB[2]));
        Assert(VectorsAreEqual(Simplex[3], SimplexA[3] - SimplexB[3]));

        struct triangle
        {
            union
            {
                struct
                {
                    vec3 A, B, C;
                    vec3 aA, aB, aC;
                    vec3 bA, bB, bC;
                };
                vec3 E[9];
            };
        };

        // NOTE(georgy): Push first 4 triangles
        std::vector<triangle> Triangles;
        Triangles.push_back({A, B, C, SimplexA[3], SimplexA[2], SimplexA[1], SimplexB[3], SimplexB[2], SimplexB[1]});
        Triangles.push_back({A, C, D, SimplexA[3], SimplexA[1], SimplexA[0], SimplexB[3], SimplexB[1], SimplexB[0]}); 
        Triangles.push_back({A, D, B, SimplexA[3], SimplexA[0], SimplexA[2], SimplexB[3], SimplexB[0], SimplexB[2]});
        Triangles.push_back({C, B, D, SimplexA[1], SimplexA[2], SimplexA[0], SimplexB[1], SimplexB[2], SimplexB[0]});

        vec3 PenetrationVector = vec3(0.0f, 0.0f, 0.0f);
        while(true)
        {
            r32 ClosestDist = FLT_MAX;
            vec3 ClosestNormal = vec3(0.0f, 0.0f, 0.0f);
            u32 ClosestTriangleIndex = 0;
            for(u32 FaceIndex = 0;
                FaceIndex < Triangles.size();
                FaceIndex++)
            {
                triangle Triangle = Triangles[FaceIndex];
                vec3 FaceNormal = Normalize(Cross(Triangle.B - Triangle.A, Triangle.C - Triangle.A));
                r32 Dist = -Dot(-Triangle.A, FaceNormal);
                if(Dist < ClosestDist)
                {
                    ClosestDist = Dist;
                    ClosestNormal = FaceNormal;
                    ClosestTriangleIndex = FaceIndex;
                }
            }

            vec3 aP = Support(ObjB, ClosestNormal);
            vec3 bP = Support(ObjA, -ClosestNormal);
            vec3 P = aP - bP;
            r32 D = Dot(P, ClosestNormal);
            if(Absolute(D - ClosestDist) < Tolerance)
            {
                triangle Triangle = Triangles[ClosestTriangleIndex];
                vec3 TriangleNormal = Normalize(Cross(Triangle.B - Triangle.A, Triangle.C - Triangle.A));
                r32 Dist = Dot(TriangleNormal, O - Triangle.A);
                vec3 ProjO = O - TriangleNormal*(Dist / LengthSq(TriangleNormal)); 

                r32 ProjOBC = Dot(TriangleNormal, Cross(Triangle.B - ProjO, Triangle.C - ProjO));
                r32 ProjOCA = Dot(TriangleNormal, Cross(Triangle.C - ProjO, Triangle.A - ProjO));
                r32 ProjOAB = Dot(TriangleNormal, Cross(Triangle.A - ProjO, Triangle.B - ProjO));

                r32 u = ProjOBC / (ProjOAB + ProjOCA + ProjOBC);
                r32 v = ProjOCA / (ProjOAB + ProjOCA + ProjOBC);
                r32 w = 1.0f - u - v;

                Assert(VectorsAreEqual(Triangle.A, Triangle.aA - Triangle.bA));
                Assert(VectorsAreEqual(Triangle.B, Triangle.aB - Triangle.bB));
                Assert(VectorsAreEqual(Triangle.C, Triangle.aC - Triangle.bC));

                // NOTE(georgy): 'A' object here is one that is A here "Support(A) - Support(B)"
                contact_point NewContact;
                NewContact.Persistent = false;
                NewContact.GlobalA = u*Triangle.bA + v*Triangle.bB + w*Triangle.bC;
                NewContact.LocalA = ObjA->RigidBody.InverseOrientation*(NewContact.GlobalA - ObjA->RigidBody.P);
                NewContact.GlobalB = u*Triangle.aA + v*Triangle.aB + w*Triangle.aC;
                NewContact.LocalB = ObjB->RigidBody.InverseOrientation*(NewContact.GlobalB - ObjB->RigidBody.P);
                NewContact.Normal = ClosestNormal;
                NewContact.Penetration = D;
                NewContact.AccumulatedImpulse = 0.0f;
                NewContact.AccumulatedImpulseTangent1 = 0.0f;
                NewContact.AccumulatedImpulseTangent2 = 0.0f;

                b32 FarEnough = true;
                for(u32 ContactIndex = 0;
                    ContactIndex < Manifold->ContactCount;
                    ContactIndex++)
                {
                    contact_point *Contact = Manifold->Contacts + ContactIndex;
                    vec3 rA = NewContact.GlobalA - Contact->GlobalA;
                    vec3 rB = NewContact.GlobalB - Contact->GlobalB;

                    b32 rATooClose = (LengthSq(rA) <= PersistentThresholdSq);
                    b32 rBTooClose = (LengthSq(rB) <= PersistentThresholdSq);

                    if(rATooClose || rBTooClose)
                    {
                        Contact->GlobalA = NewContact.GlobalA;
                        Contact->LocalA = NewContact.LocalA;
                        Contact->GlobalB = NewContact.GlobalB;
                        Contact->LocalB = NewContact.LocalB;
                        Contact->Normal = NewContact.Normal;
                        Contact->Penetration = NewContact.Penetration;

                        FarEnough = false;
                        break;
                    }
                }
                if(FarEnough)
                {
                    Manifold->Contacts[Manifold->ContactCount++] = NewContact;
                    
                    if(Manifold->ContactCount == ArrayCount(Manifold->Contacts))
                    {
                        contact_point Deepest;
                        r32 Penetration = -FLT_MAX;
                        for(u32 ContactIndex = 0;
                            ContactIndex < Manifold->ContactCount;
                            ContactIndex++)
                        {
                            contact_point *Contact = Manifold->Contacts + ContactIndex;

                            if(Contact->Penetration > Penetration)
                            {
                                Penetration = Contact->Penetration;
                                Deepest = *Contact;
                            }
                        }

                        contact_point Furthest1;
                        r32 DistanceSq1 = -FLT_MAX;
                        for(u32 ContactIndex = 0;
                            ContactIndex < Manifold->ContactCount;
                            ContactIndex++)
                        {
                            contact_point *Contact = Manifold->Contacts + ContactIndex;
                            r32 Dist = LengthSq(Contact->GlobalA - Deepest.GlobalA);

                            if(Dist > DistanceSq1)
                            {
                                DistanceSq1 = Dist;
                                Furthest1 = *Contact;
                            }
                        }

                        contact_point Furthest2;
                        r32 Distance2 = -FLT_MAX;
                        for(u32 ContactIndex = 0;
                            ContactIndex < Manifold->ContactCount;
                            ContactIndex++)
                        {
                            contact_point *Contact = Manifold->Contacts + ContactIndex;
                            r32 Dist = DistanceFromLineSegment(Contact->GlobalA, Deepest.GlobalA, Furthest1.GlobalA);

                            if(Dist > Distance2)
                            {
                                Distance2 = Dist;
                                Furthest2 = *Contact;
                            }
                        }

                        contact_point Furthest3;
                        r32 Distance3 = -FLT_MAX;
                        for(u32 ContactIndex = 0;
                            ContactIndex < Manifold->ContactCount;
                            ContactIndex++)
                        {
                            contact_point *Contact = Manifold->Contacts + ContactIndex;
                            r32 Dist = DistanceFromTriangle(Contact->GlobalA, Deepest.GlobalA, Furthest1.GlobalA, Furthest2.GlobalA);

                            if(Dist > Distance3)
                            {
                                Distance3 = Dist;
                                Furthest3 = *Contact;
                            }
                        }

                        Manifold->ContactCount = 0;
                        Manifold->Contacts[Manifold->ContactCount++] = Deepest;
                        Manifold->Contacts[Manifold->ContactCount++] = Furthest1;
                        Manifold->Contacts[Manifold->ContactCount++] = Furthest2;
                        b32 FourthPointInTriangle = PointInTriangle(Furthest3.GlobalA, Deepest.GlobalA, Furthest1.GlobalA, Furthest2.GlobalA);
                        if(!FourthPointInTriangle)
                        {
                            Manifold->Contacts[Manifold->ContactCount++] = Furthest3;
                        }
                    }
                }

                break;
            }
            else
            {
                struct edge
                {
                    vec3 Start;
                    vec3 End;

                    vec3 aStart;
                    vec3 aEnd;

                    vec3 bStart;
                    vec3 bEnd;
                };

                std::vector<edge> Edges;
                for(u32 FaceIndex = 0;
                    FaceIndex < Triangles.size();
                    )
                {
                    triangle Triangle = Triangles[FaceIndex];
                    vec3 FaceNormal = Normalize(Cross(Triangle.B - Triangle.A, Triangle.C - Triangle.A));   
                    if(Dot(FaceNormal, P - Triangle.A) >= 0.0f)
                    {
                        // TODO(georgy): Erase is slow! Change it!
                        Triangles.erase(Triangles.begin() + FaceIndex);
                        for(u32 I0 = 0, I1 = 2;
                            I0 < 3;
                            I1 = I0, I0++)
                        {
                            edge Edge = { Triangle.E[I1], Triangle.E[I0], Triangle.E[I1 + 3], Triangle.E[I0 + 3], Triangle.E[I1 + 6], Triangle.E[I0 + 6]};

                            b32 ContainsReverseEdge = false;
                            u32 ReverseEdgeIndex = 0;
                            for(u32 EdgeIndex = 0;
                                EdgeIndex < Edges.size();
                                EdgeIndex++)
                            {
                                edge ReverseEdge = { Edge.End, Edge.Start };
                                if(VectorsAreEqual(ReverseEdge.Start, Edges[EdgeIndex].Start) && 
                                   VectorsAreEqual(ReverseEdge.End, Edges[EdgeIndex].End))
                                {
                                    ContainsReverseEdge = true;
                                    ReverseEdgeIndex = EdgeIndex;
                                    break;
                                }
                            }

                            if(ContainsReverseEdge)
                            {
                                Edges.erase(Edges.begin() + ReverseEdgeIndex);
                            }
                            else
                            {
                                Edges.push_back(Edge);
                            }
                        }
                    }
                    else
                    {
                        FaceIndex++;
                    }
                }

                for(u32 EdgeIndex = 0;
                    EdgeIndex < Edges.size();
                    EdgeIndex++)
                {
                    edge Edge = Edges[EdgeIndex];
                    Triangles.push_back({Edge.Start, Edge.End, P, Edge.aStart, Edge.aEnd, aP, Edge.bStart, Edge.bEnd, bP});
                }
            }
        }
    }

    return(Result);
}

internal model
InitSphereMesh(void)
{
    u32 ParallelCount = 20;
    u32 MeridianCount = 20;
    r32 Radius = 0.5f;

    u32 VerticesCount = 0;
    vec3 Vertices[4096];

    u32 IndicesCount = 0;
    u32 Indices[4096];

    Vertices[VerticesCount++] = Radius*vec3(0.0f, 1.0f, 0.0f);
    for(u32 Parallel = 0;
        Parallel < ParallelCount;
        Parallel++)
    {
        r32 Phi = -((r32)(Parallel + 1) / (ParallelCount + 1))*PI + 0.5f*PI;
        for(u32 Meridian = 0;
            Meridian < MeridianCount;
            Meridian++)
        {
            r32 Theta = ((r32)Meridian / MeridianCount) * 2.0f*PI;
            r32 X = -Sin(Theta)*Cos(Phi);
            r32 Y = Sin(Phi);
            r32 Z = -Cos(Phi)*Cos(Theta);
            vec3 P = Radius*Normalize(vec3(X, Y, Z));
            Vertices[VerticesCount++] = P;
        }
    }
    Vertices[VerticesCount++] = Radius*vec3(0.0f, -1.0f, 0.0f);

    for(u32 Meridian = 0;
        Meridian < MeridianCount;
        Meridian++)
    {
        Indices[IndicesCount++] = 0;
        Indices[IndicesCount++] = Meridian + 1;
        Indices[IndicesCount++] = ((Meridian + 1) % MeridianCount) + 1;
    }

    for(u32 Parallel = 0;
        Parallel < ParallelCount - 1;
        Parallel++)
    {
        for(u32 Meridian = 0;
            Meridian < MeridianCount;
            Meridian++)
        {
            u32 A = (Parallel*MeridianCount + 1) + Meridian;
            u32 B = ((Parallel + 1)*MeridianCount + 1) + Meridian;
            u32 C = ((Parallel + 1)*MeridianCount + 1) + ((Meridian + 1) % MeridianCount);
            u32 D = (Parallel*MeridianCount + 1) + ((Meridian + 1) % MeridianCount);

            Indices[IndicesCount++] = A;
            Indices[IndicesCount++] = B;
            Indices[IndicesCount++] = C;
            Indices[IndicesCount++] = A;
            Indices[IndicesCount++] = C;
            Indices[IndicesCount++] = D;
        }
    }

    for(u32 Meridian = 0;
        Meridian < MeridianCount;
        Meridian++)
    {
        Indices[IndicesCount++] = ((ParallelCount - 1)*MeridianCount + 1) + ((Meridian + 1) % MeridianCount);
        Indices[IndicesCount++] = ((ParallelCount - 1)*MeridianCount + 1) + Meridian;
        Indices[IndicesCount++] = VerticesCount - 1;
    }

    model Result;
    Result.IndexCount = IndicesCount;
    Platform.InitBuffersWithEBO(&Result.Handle, sizeof(vec3)*VerticesCount, (r32 *)Vertices, 3*sizeof(r32), sizeof(u32)*IndicesCount, Indices);

    return(Result);
}

extern "C" 
GAME_UPDATE_AND_RENDER(GameUpdateAndRender)
{
    Platform = Memory->PlatformAPI;

    Assert(sizeof(game_state) <= Memory->PermanentStorageSize);
    game_state *GameState = (game_state *)Memory->PermanentStorage;
    if(!GameState->IsInitialized)
    {
        Platform.CompileShader(&GameState->Shader.Handle, "shaders/shader.glsl");
        Platform.CompileShader(&GameState->DebugShader.Handle, "shaders/DebugShader.glsl");
        
        r32 CubeVertices[] = {
            // back face
            -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f, 
            0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f, 
            0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f, 
            0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f, 
            -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f, 
            -0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f, 
            // front face
            -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f, 
            0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f, 
            0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f, 
            0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f, 
            -0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f, 
            -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f, 
            // left face
            -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f, 
            -0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f, 
            -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f, 
            -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f, 
            -0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f, 
            -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f, 
            // right face
            0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f, 
            0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f, 
            0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f, 
            0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f, 
            0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f, 
            0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f, 
            // bottom face
            -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f, 
            0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f, 
            0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f, 
            0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f, 
            -0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f, 
            -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f, 
            // top face
            -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f, 
            0.5f,  0.5f , 0.5f,  0.0f,  1.0f,  0.0f, 
            0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f, 
            0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f, 
            -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f, 
            -0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f, 
	    };

        InitModel(&GameState->Cube, sizeof(CubeVertices), CubeVertices, 36, 6*sizeof(r32));

        r32 QuadVertices[] = 
        {
            -0.5f, 0.5f, 0.0f,
            -0.5f, -0.5f, 0.0f,
            0.5f, 0.5f, 0.0f,

            0.5f, 0.5f, 0.0f,
            -0.5f, -0.5f, 0.0f,
            0.5f, -0.5f, 0.0f,
        };

        InitModel(&GameState->Quad, sizeof(QuadVertices), QuadVertices, 6, 3*sizeof(r32));

        GameState->Sphere = InitSphereMesh();

        GameState->GameObjectCount = 1;
        GameState->Hero = GameState->GameObjects + 0;
        GameState->Hero->Type = GameObject_Cube;
        GameState->Hero->Model = &GameState->Cube;
        GameState->Hero->Width = 1.0f;
        GameState->Hero->Height = 1.0f;
        GameState->Hero->Depth = 1.0f;
        GameState->Hero->RigidBody.P = vec3(4.0f, 2.0f, 4.0f);
        GameState->Hero->RigidBody.Mass = 6.0f;
        r32 OneOverTwelve = 1.0f / 12.0f;
        r32 Mass = GameState->Hero->RigidBody.Mass;
        r32 Width = GameState->Hero->Width;
        r32 Height = GameState->Hero->Height;
        r32 Depth = GameState->Hero->Depth;
        vec3 InertiaTensorMainDiagonal = vec3(OneOverTwelve*Mass*(Square(Height) + Square(Depth)), 
                                              OneOverTwelve*Mass*(Square(Width) + Square(Depth)),
                                              OneOverTwelve*Mass*(Square(Width) + Square(Height)));
        GameState->Hero->RigidBody.InertiaTensor = Scaling3x3(InertiaTensorMainDiagonal);
        GameState->Hero->RigidBody.InverseInertiaTensor = Inverse3x3(GameState->Hero->RigidBody.InertiaTensor);
        GameState->Hero->RigidBody.Orientation = Identity3x3();
        // GameState->Hero->RigidBody.Orientation = Rotation3x3(25.0f, vec3(1.0f, 0.0f, 0.0f));
        GameState->Hero->RigidBody.InverseOrientation = Transpose3x3(GameState->Hero->RigidBody.Orientation);
        GameState->Hero->RigidBody.GlobalInverseInertiaTensor = GameState->Hero->RigidBody.Orientation *  
                                                                GameState->Hero->RigidBody.InverseInertiaTensor * 
                                                                GameState->Hero->RigidBody.InverseOrientation;
        GameState->Hero->RigidBody.CoeffOfRestitution = 0.5f;
        GameState->Hero->RigidBody.CoeffOfFriction = 0.0f;

#if 1
        for(i32 ZOffset = -1;
            ZOffset <= 1;
            ZOffset++)
        {
            for(i32 YOffset = 0;
                YOffset <= 10;
                YOffset++)
            {
                for(i32 XOffset = -1;
                    XOffset <= 1;
                    XOffset++)
                {
                    game_object *GameObject = &GameState->GameObjects[GameState->GameObjectCount++];

                    GameObject->Type = GameObject_Cube;
                    GameObject->Model = &GameState->Cube;
                    GameObject->Width = 1.0f;
                    GameObject->Height = 1.0f;
                    GameObject->Depth = 1.0f;
                    GameObject->RigidBody.P = vec3(1.0f + 1.5f * XOffset, 4.0f*(YOffset) + 1.0f, 1.2f*ZOffset);
                    GameObject->RigidBody.Mass = 6.0f;
                    Mass = GameObject->RigidBody.Mass;
                    Width = GameObject->Width;
                    Height = GameObject->Height;
                    Depth = GameObject->Depth;
                    InertiaTensorMainDiagonal = vec3(OneOverTwelve*Mass*(Square(Height) + Square(Depth)), 
                                                    OneOverTwelve*Mass*(Square(Width) + Square(Depth)),
                                                    OneOverTwelve*Mass*(Square(Width) + Square(Height)));
                    GameObject->RigidBody.InertiaTensor = Scaling3x3(InertiaTensorMainDiagonal);
                    GameObject->RigidBody.InverseInertiaTensor = Inverse3x3(GameObject->RigidBody.InertiaTensor);
                    if(YOffset == 0)
                    {
                        GameObject->RigidBody.Orientation = Rotation3x3(10.0f, vec3(1.0f, 0.3f, 0.6f));
                    }
                    else
                    {
                        GameObject->RigidBody.Orientation = Rotation3x3(32.0f, vec3(1.0f, 0.3f, 0.0f));
                    }

                    // GameObject->RigidBody.Orientation = Identity3x3();
                    GameObject->RigidBody.InverseOrientation = Transpose3x3(GameObject->RigidBody.Orientation);
                    GameObject->RigidBody.GlobalInverseInertiaTensor = GameObject->RigidBody.Orientation *  
                                                                    GameObject->RigidBody.InverseInertiaTensor * 
                                                                    GameObject->RigidBody.InverseOrientation;
                    GameObject->RigidBody.CoeffOfRestitution = 0.5f;
                    GameObject->RigidBody.CoeffOfFriction = 0.125f;
                }
            }
        }
#else
        for(i32 YOffset = 0;
            YOffset <= 5;
            YOffset++)
        {
            game_object *GameObject = &GameState->GameObjects[GameState->GameObjectCount++];

            GameObject->Type = GameObject_Cube;
            GameObject->Model = &GameState->Cube;
            GameObject->Width = 1.0f;
            GameObject->Height = 1.0f;
            GameObject->Depth = 1.0f;
            GameObject->RigidBody.P = vec3(0.0f, 0.0001f + 1.0f*(YOffset), 0.0f);
            GameObject->RigidBody.Mass = 6.0f;
            Mass = GameObject->RigidBody.Mass;
            Width = GameObject->Width;
            Height = GameObject->Height;
            Depth = GameObject->Depth;
            InertiaTensorMainDiagonal = vec3(OneOverTwelve*Mass*(Square(Height) + Square(Depth)), 
                                            OneOverTwelve*Mass*(Square(Width) + Square(Depth)),
                                            OneOverTwelve*Mass*(Square(Width) + Square(Height)));
            GameObject->RigidBody.InertiaTensor = Scaling3x3(InertiaTensorMainDiagonal);
            GameObject->RigidBody.InverseInertiaTensor = Inverse3x3(GameObject->RigidBody.InertiaTensor);

            GameObject->RigidBody.Orientation = Identity3x3();
            GameObject->RigidBody.InverseOrientation = Transpose3x3(GameObject->RigidBody.Orientation);
            GameObject->RigidBody.GlobalInverseInertiaTensor = GameObject->RigidBody.Orientation *  
                                                            GameObject->RigidBody.InverseInertiaTensor * 
                                                            GameObject->RigidBody.InverseOrientation;
            GameObject->RigidBody.CoeffOfRestitution = 0.5f;
            GameObject->RigidBody.CoeffOfFriction = 0.125f;
        }
#endif

        game_object *GameObject = &GameState->GameObjects[GameState->GameObjectCount++];
        GameObject->Type = GameObject_Wall;
        GameObject->Model = &GameState->Cube;
        GameObject->Width = 1000.0f;
        GameObject->Height = 10.0f;
        GameObject->Depth = 1000.0f;
        GameObject->RigidBody.P = vec3(0.0f, -5.5f, 0.0f);
        GameObject->RigidBody.Mass = 0.0f;
        GameObject->RigidBody.Orientation = Identity3x3();
        GameObject->RigidBody.InverseOrientation = Transpose3x3(GameObject->RigidBody.Orientation);
        GameObject->RigidBody.CoeffOfRestitution = 0.5f;
        GameObject->RigidBody.CoeffOfFriction = 0.125f;   

        GameState->IsInitialized = true;
    }

    r32 dt = Input->dtForFrame;

    camera *Camera = &GameState->Camera;
    r32 CameraRotationSensetivity = 0.1f;
    Camera->CameraPitch -= Input->MouseYDisplacement*CameraRotationSensetivity;
    Camera->CameraHead -= Input->MouseXDisplacement*CameraRotationSensetivity;

    Camera->CameraPitch = Camera->CameraPitch > 89.0f ? 89.0f : Camera->CameraPitch;
    Camera->CameraPitch = Camera->CameraPitch < -89.0f ? -89.0f : Camera->CameraPitch;

    r32 PitchRadians = Radians(Camera->CameraPitch);
    r32 HeadRadians = Radians(Camera->CameraHead);

    r32 CameraDistanceFromHero = 10.0f;
    r32 FloorDistanceFromHero = CameraDistanceFromHero * Cos(-PitchRadians);
    
    r32 XOffsetFromHero = FloorDistanceFromHero * Sin(HeadRadians);
    r32 YOffsetFromHero = CameraDistanceFromHero * Sin(-PitchRadians);
    r32 ZOffsetFromHero = FloorDistanceFromHero * Cos(HeadRadians);
    vec3 CameraOffsetFromHero = vec3(XOffsetFromHero, YOffsetFromHero, ZOffsetFromHero);

    vec3 CameraForward = Normalize(-CameraOffsetFromHero);
    vec3 Forward = Normalize(vec3(CameraForward.x, 0.0f, CameraForward.z));
    vec3 CameraRight = Normalize(Cross(Forward, vec3(0.0f, 1.0f, 1.0f)));
    r32 Theta = Degrees(ATan2(Forward.z, Forward.x)) - 90.0f;
    GameState->Hero->RigidBody.ForceAccumulated = vec3(0.0f, 0.0f, 0.0f);
    if(Input->MoveForward.EndedDown)
    {
        GameState->Hero->RigidBody.ForceAccumulated += 10000.0f*dt*Forward;
    }
    if(Input->MoveBack.EndedDown)
    {
        
        GameState->Hero->RigidBody.ForceAccumulated += 10000.0f*dt*-Forward;
    }
    if(Input->MoveRight.EndedDown)
    {
        
        GameState->Hero->RigidBody.ForceAccumulated += 10000.0f*dt*CameraRight;
    }
    if(Input->MoveLeft.EndedDown)
    {
        
        GameState->Hero->RigidBody.ForceAccumulated += 10000.0f*dt*-CameraRight;
    }

    // Camera->P = GameState->Hero->RigidBody.P + CameraOffsetFromHero;
    // Camera->Dir = CameraForward;
    Camera->P = vec3(7.0f, 14.0f, 14.0f);
    Camera->Dir = -Normalize(Camera->P);

    mat4 Projection = Perspective(45.0f, (r32)WindowWidth/(r32)WindowHeight, 0.1f, 50.0f);
    mat4 View = Camera->GetRotationMatrix();

    Clear(RenderCommandBuffer, vec3(1.0f, 0.0f, 0.0f));

    PushShader(RenderCommandBuffer, GameState->DebugShader);
    PushMat4(RenderCommandBuffer, "Projection", &Projection);
    PushMat4(RenderCommandBuffer, "View", &View);

    PushShader(RenderCommandBuffer, GameState->Shader);
    PushMat4(RenderCommandBuffer, "Projection", &Projection);
    PushMat4(RenderCommandBuffer, "View", &View);
    PushVec3(RenderCommandBuffer, "CamP", Camera->P);
    PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 0.0f, 1.0f));

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < GameState->GameObjectCount;
        GameObjectIndex++)
    {
        if(GameObjectIndex != (GameState->GameObjectCount-1))
        {
            game_object *GameObject = GameState->GameObjects + GameObjectIndex;
            rigid_body *RigidBody = &GameObject->RigidBody;
        
            // RigidBody->ForceAccumulated = vec3(0.0f, 0.0f, 0.0f);
            vec3 ddP = RigidBody->ForceAccumulated * (1.0f / RigidBody->Mass);
            ddP += vec3(0.0f, -9.8f, 0.0f);
            RigidBody->dP += dt*ddP;

            RigidBody->AngularMomentum += dt*RigidBody->TorqueAccumulated;
            RigidBody->AngularVelocity = RigidBody->GlobalInverseInertiaTensor*RigidBody->AngularMomentum;
        }
    }

    for(u32 CollisionIndex = 0;
        CollisionIndex < GameState->CollisionCount;
        CollisionIndex++)
    {
        collision_data *Collision = GameState->Collisions + CollisionIndex;

        for(u32 ContactIndex = 0;
            ContactIndex < Collision->Manifold.ContactCount;
            )
        {
            contact_point *Contact = Collision->Manifold.Contacts + ContactIndex;

            vec3 NewGlobalA = Collision->A->RigidBody.P + Collision->A->RigidBody.Orientation*Contact->LocalA;
            vec3 NewGlobalB = Collision->B->RigidBody.P + Collision->B->RigidBody.Orientation*Contact->LocalB;
            Contact->GlobalA = NewGlobalA;
            Contact->GlobalB = NewGlobalB;
            
            vec3 rAB = NewGlobalB - NewGlobalA;
            if(Dot(rAB, Contact->Normal) <= -PersistentThreshold)
            {
                Collision->Manifold.Contacts[ContactIndex] = Collision->Manifold.Contacts[--Collision->Manifold.ContactCount];
            }
            else
            {
                vec3 ProjectedP = NewGlobalA + Contact->Normal*Dot(rAB, Contact->Normal);
                float DistSq = LengthSq(NewGlobalB - ProjectedP);
                if(DistSq > PersistentThresholdSq)
                {
                    Collision->Manifold.Contacts[ContactIndex] = Collision->Manifold.Contacts[--Collision->Manifold.ContactCount];
                }
                else
                {
                    ContactIndex++;
                }
            }
        }
    }

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < GameState->GameObjectCount - 1;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;

        for(u32 TestObjectIndex = GameObjectIndex + 1;
            TestObjectIndex < GameState->GameObjectCount;
            TestObjectIndex++)
        {
            game_object *TestGameObject = GameState->GameObjects + TestObjectIndex;

            collision_data *FoundCollision = 0;
            for(u32 CollisionIndex = 0;
                CollisionIndex < GameState->CollisionCount;
                CollisionIndex++)
            {
                collision_data *Collision = GameState->Collisions + CollisionIndex;

                if((Collision->A == GameObject) && (Collision->B == TestGameObject))
                {   
                    FoundCollision = Collision;
                    break;
                }
            }
            if(!FoundCollision)
            {
                FoundCollision = GameState->Collisions + GameState->CollisionCount;
                FoundCollision->A = GameObject;
                FoundCollision->B = TestGameObject;
                GameState->CollisionCount++;
            }

            if(!Intersect(&FoundCollision->Manifold, GameObject, TestGameObject))
            {
                FoundCollision->Manifold = {};
            }
        }
    }

    for(u32 CollisionIndex = 0;
        CollisionIndex < GameState->CollisionCount;
        CollisionIndex++)
    {
        collision_data *Collision = GameState->Collisions + CollisionIndex;

        rigid_body *A = &Collision->A->RigidBody;
        rigid_body *B = &Collision->B->RigidBody;

        for(u32 ContactIndex = 0;
            ContactIndex < Collision->Manifold.ContactCount;
            ContactIndex++)
        {
            contact_point *Contact = Collision->Manifold.Contacts + ContactIndex;
            if(Contact->Persistent)
            {
                r32 ImpulseMagnitude = Contact->AccumulatedImpulse;
                r32 ImpulseTangent1 = Contact->AccumulatedImpulseTangent1;
                r32 ImpulseTangent2 = Contact->AccumulatedImpulseTangent2;

                vec3 N = Contact->Normal;
                r32 AbsX = Absolute(N.x);
                r32 AbsY = Absolute(N.y);
                r32 AbsZ = Absolute(N.z);
                
                vec3 Axis = vec3(0.0f, 0.0f, 0.0f);
                if(AbsX < AbsY)
                {
                    if(AbsX < AbsZ)
                    {
                        Axis.x = 1.0f;
                    }
                    else
                    {
                        Axis.z = 1.0f;
                    }
                }
                else if(AbsY < AbsZ)
                {
                    Axis.y = 1.0f;
                }
                else
                {
                    Axis.z = 1.0f;
                }
                vec3 Tangent1 = Normalize(Cross(N, Axis));
                vec3 Tangent2 = Normalize(Cross(N, Tangent1));

                vec3 aP = A->Orientation*Contact->LocalA;
                r32 OneOverMassA = 1.0f / A->Mass;
                A->dP += (ImpulseMagnitude * OneOverMassA)*N +
                         (ImpulseTangent1 * OneOverMassA)*Tangent1 +
                         (ImpulseTangent2 * OneOverMassA)*Tangent2;
                A->AngularMomentum += Cross(aP, ImpulseMagnitude*N) +
                                      Cross(aP, ImpulseTangent1*Tangent1) + 
                                      Cross(aP, ImpulseTangent2*Tangent2);
                A->AngularVelocity = A->GlobalInverseInertiaTensor*A->AngularMomentum;

                if(B->Mass != 0.0f)
                {
                    vec3 bP = B->Orientation*Contact->LocalB;
                    r32 OneOverMassB = 1.0f / B->Mass;
                    B->dP += (-ImpulseMagnitude * OneOverMassB)*N + 
                             (-ImpulseTangent1 * OneOverMassB)*Tangent1 +
                             (-ImpulseTangent2 * OneOverMassB)*Tangent2;;
                    B->AngularMomentum += Cross(bP, -ImpulseMagnitude*N) +
                                          Cross(bP, -ImpulseTangent1*Tangent1) + 
                                          Cross(bP, -ImpulseTangent2*Tangent2);
                    B->AngularVelocity = B->GlobalInverseInertiaTensor*B->AngularMomentum;
                }
            }
        }
    }

    for(u32 Iter = 0;
        Iter < 20;
        Iter++)
    {
        for(u32 CollisionIndex = 0;
            CollisionIndex < GameState->CollisionCount;
            CollisionIndex++)
        {
            collision_data *Collision = GameState->Collisions + CollisionIndex;

            collision_manifold *Manifold = &Collision->Manifold;

            rigid_body *A = &Collision->A->RigidBody;
            rigid_body *B = &Collision->B->RigidBody;
            
            r32 CoeffOfRestitution = Min(A->CoeffOfRestitution, B->CoeffOfRestitution);
            r32 CoeffOfFriction = SquareRoot(A->CoeffOfFriction*B->CoeffOfFriction);

            for(u32 ContactIndex = 0;
                ContactIndex < Collision->Manifold.ContactCount;
                ContactIndex++)
            {
                contact_point *Contact = Collision->Manifold.Contacts + ContactIndex;
                vec3 N = Contact->Normal;

                r32 Slop = 0.01f;
                r32 BiasFactor = 0.2f;
                r32 BiasVelocity = (BiasFactor / dt)*Max(Contact->Penetration - Slop, 0.0f);

                r32 AbsX = Absolute(N.x);
                r32 AbsY = Absolute(N.y);
                r32 AbsZ = Absolute(N.z);
                
                vec3 Axis = vec3(0.0f, 0.0f, 0.0f);
                if(AbsX < AbsY)
                {
                    if(AbsX < AbsZ)
                    {
                        Axis.x = 1.0f;
                    }
                    else
                    {
                        Axis.z = 1.0f;
                    }
                }
                else if(AbsY < AbsZ)
                {
                    Axis.y = 1.0f;
                }
                else
                {
                    Axis.z = 1.0f;
                }

                vec3 Tangent1 = Normalize(Cross(N, Axis));
                vec3 Tangent2 = Normalize(Cross(N, Tangent1));

                vec3 aP = Collision->A->RigidBody.Orientation*Contact->LocalA;
                vec3 VelocityA = A->dP + Cross(A->AngularVelocity, aP);
                r32 OneOverMassA = 1.0f / A->Mass;

                r32 ImpulseMagnitude = 0.0f;
                if(B->Mass == 0.0f)
                {
                    r32 ImpulseNom = -(1.0f + CoeffOfRestitution)*Dot(VelocityA, N) + BiasVelocity;
                    r32 ImpulseDenom = OneOverMassA + Dot(N, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, N), aP)));
                    ImpulseMagnitude = ImpulseNom / ImpulseDenom;
                    
                    r32 Temp = Contact->AccumulatedImpulse;
                    Contact->AccumulatedImpulse = Max(Contact->AccumulatedImpulse + ImpulseMagnitude, 0.0f);
                    ImpulseMagnitude = Contact->AccumulatedImpulse - Temp;
                }
                else
                {
                    vec3 bP = Collision->B->RigidBody.Orientation*Contact->LocalB;
                    vec3 VelocityB = B->dP + Cross(B->AngularVelocity, bP);
                    r32 OneOverMassB = 1.0f / B->Mass;

                    vec3 RelativeVelocity = VelocityA - VelocityB;
                    r32 ImpulseNom = -(1.0f + CoeffOfRestitution)*Dot(RelativeVelocity, N) + BiasVelocity;
                    r32 ImpulseDenom = (OneOverMassA + OneOverMassB)+
                                    Dot(N, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, N), aP))) + 
                                    Dot(N, (Cross(B->GlobalInverseInertiaTensor*Cross(bP, N), bP)));
                    ImpulseMagnitude = ImpulseNom / ImpulseDenom;

                    r32 Temp = Contact->AccumulatedImpulse;
                    Contact->AccumulatedImpulse = Max(Contact->AccumulatedImpulse + ImpulseMagnitude, 0.0f);
                    ImpulseMagnitude = Contact->AccumulatedImpulse - Temp;

                    B->dP += (-ImpulseMagnitude * OneOverMassB)*N;
                    B->AngularMomentum += Cross(bP, -ImpulseMagnitude*N);
                    B->AngularVelocity = B->GlobalInverseInertiaTensor*B->AngularMomentum;
                }
                
                A->dP += (ImpulseMagnitude * OneOverMassA)*N;
                A->AngularMomentum += Cross(aP, ImpulseMagnitude*N);
                A->AngularVelocity = A->GlobalInverseInertiaTensor*A->AngularMomentum;

                r32 ImpulseTangent1 = 0.0f;
                r32 ImpulseTangent2 = 0.0f;
                VelocityA = A->dP + Cross(A->AngularVelocity, aP);
                if(B->Mass == 0.0f)
                {
                    r32 ImpulseNomTangent1 = -(1.0f + CoeffOfRestitution)*Dot(VelocityA, Tangent1);
                    r32 ImpulseDenomTangent1 = OneOverMassA + Dot(Tangent1, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, Tangent1), aP)));   
                    ImpulseTangent1 = ImpulseNomTangent1 / ImpulseDenomTangent1;

                    r32 Temp = Contact->AccumulatedImpulseTangent1;
                    Contact->AccumulatedImpulseTangent1 = Clamp(Contact->AccumulatedImpulseTangent1 + ImpulseTangent1,
                                                                -CoeffOfFriction*Contact->AccumulatedImpulse,
                                                                CoeffOfFriction*Contact->AccumulatedImpulse);
                    ImpulseTangent1 = Contact->AccumulatedImpulseTangent1 - Temp;

                    r32 ImpulseNomTangent2 = -(1.0f + CoeffOfRestitution)*Dot(VelocityA, Tangent2);
                    r32 ImpulseDenomTangent2 = OneOverMassA + Dot(Tangent2, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, Tangent2), aP)));   
                    ImpulseTangent2 = ImpulseNomTangent2 / ImpulseDenomTangent2;

                    Temp = Contact->AccumulatedImpulseTangent2;
                    Contact->AccumulatedImpulseTangent2 = Clamp(Contact->AccumulatedImpulseTangent2 + ImpulseTangent2,
                                                                -CoeffOfFriction*Contact->AccumulatedImpulse,
                                                                CoeffOfFriction*Contact->AccumulatedImpulse);
                    ImpulseTangent2 = Contact->AccumulatedImpulseTangent2 - Temp;
                }
                else
                {
                    vec3 bP = Collision->B->RigidBody.Orientation*Contact->LocalB;
                    vec3 VelocityB = B->dP + Cross(B->AngularVelocity, bP);
                    r32 OneOverMassB = 1.0f / B->Mass;

                    vec3 RelativeVelocity = VelocityA - VelocityB;

                    r32 ImpulseNomTangent1 = -(1.0f + CoeffOfRestitution)*Dot(VelocityA - VelocityB, Tangent1);
                    r32 ImpulseDenomTangent1 = (OneOverMassA + OneOverMassB)+
                                                Dot(Tangent1, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, Tangent1), aP))) + 
                                                Dot(Tangent1, (Cross(B->GlobalInverseInertiaTensor*Cross(bP, Tangent1), bP)));
                    ImpulseTangent1 = ImpulseNomTangent1 / ImpulseDenomTangent1;

                    r32 Temp = Contact->AccumulatedImpulseTangent1;
                    Contact->AccumulatedImpulseTangent1 = Clamp(Contact->AccumulatedImpulseTangent1 + ImpulseTangent1,
                                                                -CoeffOfFriction*Contact->AccumulatedImpulse,
                                                                CoeffOfFriction*Contact->AccumulatedImpulse);
                    ImpulseTangent1 = Contact->AccumulatedImpulseTangent1 - Temp;

                    r32 ImpulseNomTangent2 = -(1.0f + CoeffOfRestitution)*Dot(VelocityA - VelocityB, Tangent2);
                    r32 ImpulseDenomTangent2 = (OneOverMassA + OneOverMassB)+
                                                Dot(Tangent2, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, Tangent2), aP))) + 
                                                Dot(Tangent2, (Cross(B->GlobalInverseInertiaTensor*Cross(bP, Tangent2), bP)));   
                    ImpulseTangent2 = ImpulseNomTangent2 / ImpulseDenomTangent2;

                    Temp = Contact->AccumulatedImpulseTangent2;
                    Contact->AccumulatedImpulseTangent2 = Clamp(Contact->AccumulatedImpulseTangent2 + ImpulseTangent2,
                                                                -CoeffOfFriction*Contact->AccumulatedImpulse,
                                                                CoeffOfFriction*Contact->AccumulatedImpulse);
                    ImpulseTangent2 = Contact->AccumulatedImpulseTangent2 - Temp;

                    B->dP += (-ImpulseTangent1 * OneOverMassB)*Tangent1 +
                             (-ImpulseTangent2 * OneOverMassB)*Tangent2;
                    B->AngularMomentum += Cross(bP, -ImpulseTangent1*Tangent1) +
                                          Cross(bP, -ImpulseTangent2*Tangent2);
                    B->AngularVelocity = B->GlobalInverseInertiaTensor*B->AngularMomentum;
                }
                
                A->dP += (ImpulseTangent1 * OneOverMassA)*Tangent1 + 
                         (ImpulseTangent2 * OneOverMassA)*Tangent2;
                         
                A->AngularMomentum += Cross(aP, ImpulseTangent1*Tangent1) +
                                      Cross(aP, ImpulseTangent2*Tangent2);
                A->AngularVelocity = A->GlobalInverseInertiaTensor*A->AngularMomentum;
            }
        }
    }
    
    for(u32 GameObjectIndex = 0;
        GameObjectIndex < GameState->GameObjectCount;
        GameObjectIndex++)
    {
        if(GameObjectIndex != (GameState->GameObjectCount-1))
        {
            game_object *GameObject = GameState->GameObjects + GameObjectIndex;
            rigid_body *RigidBody = &GameObject->RigidBody;
            
            RigidBody->P += dt*RigidBody->dP;
            RigidBody->Orientation = Rotation3x3(dt*Degrees(Length(RigidBody->AngularVelocity)), 
                                                 RigidBody->AngularVelocity) *
                                     RigidBody->Orientation;

            // NOTE(georgy): Reorthogonalize orientation matrix 
            vec3 XAxis = vec3(RigidBody->Orientation.a11, RigidBody->Orientation.a21, RigidBody->Orientation.a31);
            vec3 YAxis = vec3(RigidBody->Orientation.a12, RigidBody->Orientation.a22, RigidBody->Orientation.a32);
            vec3 ZAxis = vec3(RigidBody->Orientation.a13, RigidBody->Orientation.a23, RigidBody->Orientation.a33);

            XAxis = Normalize(XAxis);
            YAxis = Normalize(YAxis - XAxis*Dot(YAxis, XAxis));
            ZAxis = Normalize(ZAxis - XAxis*Dot(ZAxis, XAxis) - YAxis*Dot(ZAxis, YAxis));
            
            RigidBody->Orientation.a11 = XAxis.x;
            RigidBody->Orientation.a21 = XAxis.y;
            RigidBody->Orientation.a31 = XAxis.z;
            RigidBody->Orientation.a12 = YAxis.x;
            RigidBody->Orientation.a22 = YAxis.y;
            RigidBody->Orientation.a33 = YAxis.z;
            RigidBody->Orientation.a13 = ZAxis.x;
            RigidBody->Orientation.a23 = ZAxis.y;
            RigidBody->Orientation.a33 = ZAxis.z;
            
            RigidBody->InverseOrientation = Transpose3x3(RigidBody->Orientation);
            RigidBody->GlobalInverseInertiaTensor = RigidBody->Orientation * 
                                                    RigidBody->InverseInertiaTensor * 
                                                    RigidBody->InverseOrientation;
        }
    }

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < GameState->GameObjectCount;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;

        if(GameObjectIndex == (GameState->GameObjectCount - 1))
        {
            PushVec3(RenderCommandBuffer, "Color", vec3(1.0f, 1.0f, 0.0f));
        }
        else
        {
            PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 0.0f, 1.0f));
        }

        mat4 Model = Translation(GameObject->RigidBody.P) * 
                        Mat4(GameObject->RigidBody.Orientation) * 
                        Scaling(vec3(GameObject->Width, GameObject->Height, GameObject->Depth));
        PushMat4(RenderCommandBuffer, "Model", &Model);
        DrawModel(RenderCommandBuffer, GameObject->Model);
    }

    PushShader(RenderCommandBuffer, GameState->DebugShader);
    DisableDepthTest(RenderCommandBuffer);
    for(u32 CollisionIndex = 0;
        CollisionIndex < GameState->CollisionCount;
        CollisionIndex++)
    {
        collision_data *Collision = GameState->Collisions + CollisionIndex;

        for(u32 ContactIndex = 0;
            ContactIndex < Collision->Manifold.ContactCount;
            ContactIndex++)
        {
            contact_point *Contact = Collision->Manifold.Contacts + ContactIndex;

            mat4 Model = Translation(Contact->GlobalA) * Scaling(0.1f);
            PushMat4(RenderCommandBuffer, "Model", &Model);
            PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 1.0f, 0.0f));
            DrawModelEBO(RenderCommandBuffer, &GameState->Sphere);

            // Model = Translation(Contact->GlobalB) * Scaling(0.1f);
            // PushMat4(RenderCommandBuffer, "Model", &Model);
            // PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 1.0f, 1.0f));
            // DrawModelEBO(RenderCommandBuffer, &GameState->Sphere);
        }
    }
    EnableDepthTest(RenderCommandBuffer);
}
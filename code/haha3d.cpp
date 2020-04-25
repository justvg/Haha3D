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

internal u32
CalculateMinkowskiDiff(vec3 *Diff, u32 ACount, vec3 *AVertices, u32 BCount, vec3 *BVertices)
{
    u32 VertexCount = 0;

    for(u32 AVertex = 0;
        AVertex < ACount;
        AVertex++)
    {
        for(u32 BVertex = 0;
            BVertex < BCount;
            BVertex++)
        {
            vec3 NewVertex = AVertices[AVertex] - BVertices[BVertex];
            if(!Contains(Diff, VertexCount, NewVertex))
            {
                Diff[VertexCount++] = NewVertex;
            }
        }
    }

    return(VertexCount);
}

internal u32
ConstructConvexHull(vec2 *CH, u32 Count, vec2 *Vertices)
{
    u32 UpperChainCount = 2;
    vec2 UpperChain[32];
    UpperChain[0] = Vertices[0];
    UpperChain[1] = Vertices[1];
    vec2 LastUpperChainEdge = UpperChain[1] - UpperChain[0];
    for(u32 Index = 2;
        Index < Count;
        Index++)
    {
        vec2 Vertex = Vertices[Index];

        while((UpperChainCount > 1) && 
              (Cross2D(Vertex - UpperChain[UpperChainCount - 2], LastUpperChainEdge) < 0.0f))
        {
            // NOTE(georgy): Vertex lies to the left of the LastUpperChainEdge
            UpperChainCount--;
            if(UpperChainCount > 1)
            {
                LastUpperChainEdge = UpperChain[UpperChainCount - 1] - UpperChain[UpperChainCount - 2];
            }
        }

        UpperChain[UpperChainCount++] = Vertex;
        LastUpperChainEdge = UpperChain[UpperChainCount - 1] - UpperChain[UpperChainCount - 2];
    }

    u32 LowerChainCount = 2;
    vec2 LowerChain[32];
    LowerChain[0] = Vertices[0];
    LowerChain[1] = Vertices[1];
    vec2 LastLowerChainEdge = LowerChain[1] - LowerChain[0];
    for(u32 Index = 2;
        Index < Count;
        Index++)
    {
        vec2 Vertex = Vertices[Index];

        while((LowerChainCount > 1) && 
              (Cross2D(Vertex - LowerChain[LowerChainCount - 2], LastLowerChainEdge) > 0.0f))
        {
            // NOTE(georgy): Vertex lies to the left of the LastLowerChainEdge
            LowerChainCount--;
            if(LowerChainCount > 1)
            {
                LastLowerChainEdge = LowerChain[LowerChainCount - 1] - LowerChain[LowerChainCount - 2];
            }
        }

        LowerChain[LowerChainCount++] = Vertex;
        LastLowerChainEdge = LowerChain[LowerChainCount - 1] - LowerChain[LowerChainCount - 2];
    }

    u32 CHCount = 0;
    for(u32 Index = 0;
        Index < LowerChainCount;
        Index++)
    {
        CH[CHCount++] = LowerChain[Index];
    }

    Assert(UpperChainCount >= 2);
    for(u32 Index = UpperChainCount - 2;
        Index > 0;
        Index--)
    {
        CH[CHCount++] = UpperChain[Index];
    }

    return(CHCount);
}

struct contact_configuration
{
    r32 Min, Max;
    u32 Index[2];
    char Type[2];
};

internal b32 
NoIntersection(contact_configuration ConfigA, contact_configuration ConfigB, 
               r32 Speed, r32 &tFirst, r32 &tLast, i32 &Side,
               contact_configuration *FirstContactConfigA, contact_configuration *FirstContactConfigB)
{
    if(ConfigA.Max < ConfigB.Min)
    {
        // NOTE(georgy): Interval of A is initially on 'left' of interval B
        
        if(Speed <= 0.0f) { return(true); }

        r32 t = (ConfigB.Min - ConfigA.Max) / Speed;
        if(t > tFirst)
        {
            tFirst = t;
            Side = -1;
            *FirstContactConfigA = ConfigA;
            *FirstContactConfigB = ConfigB;
        }

        t = (ConfigB.Max - ConfigA.Min) / Speed;
        if(t < tLast)
        {
            tLast = t;
        }
    }
    else if(ConfigB.Max < ConfigA.Min)
    {
        // NOTE(georgy): Interval of A is initially on 'right' of interval B

        if(Speed >= 0) { return(true); }

        r32 t = (ConfigB.Max - ConfigA.Min) / Speed;
        if(t > tFirst)
        {
            tFirst = t;
            Side = 1;
            *FirstContactConfigA = ConfigA;
            *FirstContactConfigB = ConfigB;
        }

        t = (ConfigB.Min - ConfigA.Max) / Speed;
        if(t < tLast)
        {
            tLast = t;
        }
    }
    else
    {
        // NOTE(georgy): Interval of A and interval of B initially overlap

        if(Speed > 0)
        {
            r32 t = (ConfigB.Max - ConfigA.Min) / Speed;
            if(t < tLast)
            {
                tLast = t;
            }
        }
        else if(Speed < 0)
        {
            r32 t = (ConfigB.Min - ConfigA.Max) / Speed;
            if(t < tLast)
            {
                tLast = t;
            }
        }
    }

    b32 NoIntersect = (tFirst > tLast);
    return(NoIntersect);
}

// NOTE(georgy): This function exploits the fact that it projects the polygon on one of its edge normals
internal void
ComputeProjectionIntervalForNormal(vec2 *Vertices, u32 Count, u32 EdgeIndex, vec2 Normal, contact_configuration *Config)
{
    Config->Max = Dot(Normal, Vertices[EdgeIndex]);
    Config->Index[1] = EdgeIndex;
    Config->Type[1] = 'E'; 
    Config->Type[0] = 'V'; 
 
    Config->Min = Config->Max;
    for(u32 I = 0, VertexIndex = (EdgeIndex + 2) % Count;
        I < (Count - 2);
        I++, VertexIndex = (VertexIndex + 1) % Count)
    {
        r32 Value = Dot(Normal, Vertices[VertexIndex]);
        if(Absolute(Value - Config->Min) <= Epsilon)
        {
            // NOTE(georgy): Found an edge parallel to initial projected edge
            Config->Type[0] = 'E';
            break;
        }
        else if(Value < Config->Min)
        {
            Config->Min = Value;
            Config->Index[0] = VertexIndex;
        }
        else
        {
            // NOTE(georgy): When dot product becomes larger than the min., we are walking back towards initial edge
            break;
        }
    }
}

internal void
ComputeProjectionIntervalGeneral(vec2 *Vertices, vec2 *Edges, u32 Count, vec2 V, contact_configuration *Config)
{
    Config->Min = Config->Max = Dot(V, Vertices[0]);
    Config->Index[0] = Config->Index[1] = 0;
    Config->Type[0] = Config->Type[1] = 'V';

    for(u32 Vertex = 1;
        Vertex < Count;
        Vertex++)
    {
        r32 Value = Dot(V, Vertices[Vertex]);
        if(Value < Config->Min)
        {
            Config->Min = Value;
            Config->Index[0] = Vertex;
        }
        else if(Value > Config->Max)
        {
            Config->Max = Value;
            Config->Index[1] = Vertex;
        }
    }

    for(u32 I = 0; I < 2; I++)
    {
        u32 FirstTestEdgeIndex = (Config->Index[I] == 0) ? (Count - 1) : (Config->Index[I] - 1);
        u32 SecondTestEdgeIndex = Config->Index[I];
        if(Dot(V, Edges[FirstTestEdgeIndex]) == 0.0f)
        {
            Config->Index[I] = FirstTestEdgeIndex;
            Config->Type[I] = 'E';
        }
        else if(Dot(V, Edges[SecondTestEdgeIndex]) == 0.0f)
        {
            Config->Type[I] = 'E';
        }
    }
}

#if 0
internal b32
TestIntersection(game_object *A, vec2 DeltaP, game_object *B, r32 *t, vec2 &PointOfContact, vec2 &CollisionNormal,
                 u32 MaxDepth = 100)
{
    b32 Result = false;

    if(MaxDepth > 0)
    {
        // NOTE(georgy): Process as if B is stationary, A is moving
        r32 tFirst = 0.0f;
        r32 tLast = FLT_MAX;

        vec2 AAxisX = vec2(Cos(Radians(A->RigidBody.Orientation)), Sin(Radians(A->RigidBody.Orientation)));
        vec2 AAxisY = Perp(AAxisX);
        // TODO(georgy): This function works for any convex polygon!
        // NOTE(georgy): CCW order
        vec2 APolygonVertices[4] = 
        {
            vec2(A->RigidBody.P + 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY),
            vec2(A->RigidBody.P - 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY),
            vec2(A->RigidBody.P - 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY),
            vec2(A->RigidBody.P + 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY),
        };

        vec2 APolygonEdges[4]
        {
            APolygonVertices[1] - APolygonVertices[0],
            APolygonVertices[2] - APolygonVertices[1],
            APolygonVertices[3] - APolygonVertices[2],
            APolygonVertices[0] - APolygonVertices[3],
        };

        vec2 BAxisX = vec2(Cos(Radians(B->RigidBody.Orientation)), Sin(Radians(B->RigidBody.Orientation)));
        vec2 BAxisY = Perp(BAxisX);
        // TODO(georgy): This function works for any convex polygon!
        // NOTE(georgy): CCW order
        vec2 BPolygonVertices[4] = 
        {
            vec2(B->RigidBody.P + 0.5f*B->Width*BAxisX + 0.5f*B->Height*BAxisY),
            vec2(B->RigidBody.P - 0.5f*B->Width*BAxisX + 0.5f*B->Height*BAxisY),
            vec2(B->RigidBody.P - 0.5f*B->Width*BAxisX - 0.5f*B->Height*BAxisY),
            vec2(B->RigidBody.P + 0.5f*B->Width*BAxisX - 0.5f*B->Height*BAxisY),
        };

        vec2 BPolygonEdges[4]
        {
            BPolygonVertices[1] - BPolygonVertices[0],
            BPolygonVertices[2] - BPolygonVertices[1],
            BPolygonVertices[3] - BPolygonVertices[2],
            BPolygonVertices[0] - BPolygonVertices[3],
        };

        contact_configuration FirstContactConfigA = {};
        contact_configuration FirstContactConfigB = {};
        i32 Side = 0;

        // NOTE(georgy): Test normals of A's edges for separation
        for(u32 EdgeIndex = 0;
            EdgeIndex < ArrayCount(APolygonEdges);
            EdgeIndex++)
        {
            vec2 Edge = APolygonEdges[EdgeIndex];
            vec2 Normal = Normalize(vec2(Edge.y, -Edge.x));
            r32 Speed = Dot(DeltaP, Normal);

            contact_configuration ConfigA, ConfigB;
            ComputeProjectionIntervalForNormal(APolygonVertices, ArrayCount(APolygonVertices), EdgeIndex, Normal, &ConfigA);
            ComputeProjectionIntervalGeneral(BPolygonVertices, BPolygonEdges, ArrayCount(BPolygonVertices), Normal, &ConfigB);

            if(NoIntersection(ConfigA, ConfigB, Speed, tFirst, tLast, Side, &FirstContactConfigA, &FirstContactConfigB))
            {
                return(false);
            }
        }

        // NOTE(georgy): Test normals of B's edges for separation
        for(u32 EdgeIndex = 0;
            EdgeIndex < ArrayCount(BPolygonEdges);
            EdgeIndex++)
        {
            vec2 Edge = BPolygonEdges[EdgeIndex];
            vec2 Normal = Normalize(vec2(Edge.y, -Edge.x));
            r32 Speed = Dot(DeltaP, Normal);

            contact_configuration ConfigA, ConfigB;
            ComputeProjectionIntervalGeneral(APolygonVertices, APolygonEdges, ArrayCount(APolygonVertices), Normal, &ConfigA);
            ComputeProjectionIntervalForNormal(BPolygonVertices, ArrayCount(BPolygonVertices), EdgeIndex, Normal, &ConfigB);

            if(NoIntersection(ConfigA, ConfigB, Speed, tFirst, tLast, Side, &FirstContactConfigA, &FirstContactConfigB))
            {
                return(false);
            }
        }

        Result = true;
        *t = tFirst;

        if(Side == 1)
        {
            // NOTE(georgy): A-min meets B-max

            if(FirstContactConfigA.Type[0] == 'V')
            {
                // NOTE(georgy): vertex-vertex or vertex-edge intersection
                PointOfContact = APolygonVertices[FirstContactConfigA.Index[0]] + tFirst*DeltaP;

                if(FirstContactConfigB.Type[1] == 'E')
                {
                    vec2 BContactEdge = BPolygonEdges[FirstContactConfigB.Index[1]];
                    CollisionNormal = Normalize(vec2(BContactEdge.y, -BContactEdge.x));
                }
                else
                {
                    u32 BContactEdgeIndex = (FirstContactConfigB.Index[1] == 0) ? ArrayCount(BPolygonEdges) - 1 : FirstContactConfigB.Index[1] - 1;
                    vec2 BContactEdge = BPolygonEdges[BContactEdgeIndex];
                    CollisionNormal = Normalize(vec2(BContactEdge.y, -BContactEdge.x));
                }
            }
            else if(FirstContactConfigB.Type[1] == 'V')
            {
                // NOTE(georgy): edge-vertex intersection
                PointOfContact = BPolygonVertices[FirstContactConfigB.Index[1]];

                vec2 AContactEdge = APolygonEdges[FirstContactConfigA.Index[0]];
                CollisionNormal = Normalize(Perp(AContactEdge));
            }
            else
            {
                // NOTE(georgy): edge-edge intersection
                vec2 P = BPolygonVertices[FirstContactConfigB.Index[1]];
                vec2 E = BPolygonEdges[FirstContactConfigB.Index[1]];
                vec2 U0 = APolygonVertices[FirstContactConfigA.Index[0]];
                vec2 U1 = APolygonVertices[(FirstContactConfigA.Index[0] + 1) % ArrayCount(APolygonVertices)];

                r32 S0 = Dot(E, U1 - P) / LengthSq(E);
                r32 S1 = Dot(E, U0 - P) / LengthSq(E);

                // NOTE(georgy): Find interval intersection
                r32 IntervalMinT = Max(0, S0);
                r32 IntervalMaxT = Min(1, S1);
                r32 Middle = 0.5f*(IntervalMinT + IntervalMaxT);

                PointOfContact = P + Middle*E;

                CollisionNormal = Normalize(vec2(E.y, -E.x));
            }
        }
        else if(Side == -1)
        {
            // NOTE(georgy): A-max meets B-min

            if(FirstContactConfigA.Type[1] == 'V')
            {
                // NOTE(georgy): vertex-vertex or vertex-edge intersection
                PointOfContact = APolygonVertices[FirstContactConfigA.Index[1]] + tFirst*DeltaP; 

                if(FirstContactConfigB.Type[0] == 'E')
                {
                    vec2 BContactEdge = BPolygonEdges[FirstContactConfigB.Index[0]];
                    CollisionNormal = Normalize(vec2(BContactEdge.y, -BContactEdge.x));
                }
                else
                {
                    u32 BContactEdgeIndex = (FirstContactConfigB.Index[0] == 0) ? ArrayCount(BPolygonEdges) - 1 : FirstContactConfigB.Index[0] - 1;
                    vec2 BContactEdge = BPolygonEdges[BContactEdgeIndex];
                    CollisionNormal = Normalize(vec2(BContactEdge.y, -BContactEdge.x));
                }
            }
            else if(FirstContactConfigB.Type[0] == 'V')
            {
                // NOTE(georgy): edge-vertex intersection
                PointOfContact = BPolygonVertices[FirstContactConfigB.Index[0]];

                vec2 AContactEdge = APolygonEdges[FirstContactConfigA.Index[1]];
                CollisionNormal = Normalize(Perp(AContactEdge));
            }
            else
            {
                // NOTE(georgy): edge-edge intersection
                vec2 P = APolygonVertices[FirstContactConfigA.Index[1]];
                vec2 E = APolygonEdges[FirstContactConfigA.Index[1]];
                vec2 U0 = BPolygonVertices[FirstContactConfigB.Index[0]];
                vec2 U1 = BPolygonVertices[(FirstContactConfigB.Index[0] + 1) % ArrayCount(BPolygonVertices)];

                r32 S0 = Dot(E, U1 - P) / LengthSq(E);
                r32 S1 = Dot(E, U0 - P) / LengthSq(E);

                // NOTE(georgy): Find interval intersection
                r32 IntervalMinT = Max(0, S0);
                r32 IntervalMaxT = Min(1, S1);
                r32 Middle = 0.5f*(IntervalMinT + IntervalMaxT);

                PointOfContact = P + Middle*E;

                CollisionNormal = Normalize(vec2(E.y, -E.x));
            }
        }
        else
        {
            vec2 MinkowskiDiff[32];
            u32 MinkowskiDiffVertexCount = CalculateMinkowskiDiff(MinkowskiDiff,
                                                                ArrayCount(BPolygonVertices), BPolygonVertices, 
                                                                ArrayCount(APolygonVertices), APolygonVertices);

            // NOTE(georgy): Sort by X. If Xs are equal, sort by Y
            for(u32 Outer = 0;
                Outer < MinkowskiDiffVertexCount; 
                Outer++)
            {
                for(u32 Inner = 0; 
                    Inner < (MinkowskiDiffVertexCount - 1); 
                    Inner++)
                {
                    vec2 *A = &MinkowskiDiff[Inner];
                    vec2 *B = &MinkowskiDiff[Inner + 1];
                    b32 XIsEqual = (Absolute(A->x - B->x) <= Epsilon);
                    if(XIsEqual)
                    {
                        if(A->y > B->y)
                        {
                            vec2 Temp = *A;
                            *A = *B;
                            *B = Temp;
                        }
                    }
                    else if(A->x > B->x)
                    {
                        vec2 Temp = *A;
                        *A = *B;
                        *B = Temp;
                    }
                }
            }

            vec2 CH[32];
            u32 CHCount = ConstructConvexHull(CH, MinkowskiDiffVertexCount, MinkowskiDiff);
            r32 MinLength = FLT_MAX;
            vec2 MinimalTranslationalDistance = {};
            for(u32 I0 = 0, I1 = CHCount - 1;
                I0 < CHCount;
                I1 = I0, I0++)
            {
                vec2 Edge = CH[I0] - CH[I1];
                vec2 EdgeInwardNormal = Perp(Edge);
                vec2 PointOnEdge = CH[I1];
                vec2 ClosestPointOnEdgeFromOrigin = -(Dot(EdgeInwardNormal, -PointOnEdge) / LengthSq(EdgeInwardNormal))*EdgeInwardNormal;
                
                r32 Len = Length(ClosestPointOnEdgeFromOrigin);
                if(Len < MinLength)
                {
                    MinLength = Len;
                    MinimalTranslationalDistance = ClosestPointOnEdgeFromOrigin;
                }
            }

            A->RigidBody.P = A->RigidBody.P + 1.2f*MinimalTranslationalDistance;
            Result = TestIntersection(A, DeltaP, B, t, PointOfContact, CollisionNormal, MaxDepth - 1);
        }
    }

    return(Result);
}
#endif

internal vec3
Support(vec3 D, vec3 SphereP, r32 SphereRadius)
{
    vec3 Result = SphereP + SphereRadius*Normalize(D);

    return(Result);
}

internal vec3 
Support(vec3 D, u32 Count, vec3 *Points)
{
    Assert(Count > 0);

    vec3 Result = Points[0];
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
                           vec3 A, vec3 B, vec3 C,
                           u32 *SimplexCount, vec3 *Simplex, vec3 *Dir)
{
    switch(TetrahedronVR.VR)
    {
        case VoronoiRegion_Vertex:
        {
            *SimplexCount = 1;
            Simplex[0] = TetrahedronVR.P; 
            *Dir = -TetrahedronVR.P;
        } break;

        case VoronoiRegion_Edge:
        {
            *SimplexCount = 2;
            Simplex[0] = TetrahedronVR.EdgeP0;
            Simplex[1] = TetrahedronVR.EdgeP1;
            *Dir = -TetrahedronVR.P;
        } break;

        case VoronoiRegion_Triangle:
        {
            *SimplexCount = 3;
            Simplex[0] = C;
            Simplex[1] = B;
            Simplex[2] = A;
            *Dir = -TetrahedronVR.P;
        } break;
    }
}

internal b32
DoSimplex(u32 *SimplexCount, vec3 *Simplex, vec3 *Dir)
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
                *Dir = -TriangleVR.P;
            } break;

            case VoronoiRegion_Edge:
            {
                *SimplexCount = 2;
                Simplex[0] = TriangleVR.EdgeP0;
                Simplex[1] = TriangleVR.EdgeP1;
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

        b32 InsideABCPlane = (Dot(ABCNormal, AO) * Dot(ABCNormal, AD) >= 0.0f);
        b32 InsideACDPlane = (Dot(ACDNormal, AO) * Dot(ACDNormal, AB) >= 0.0f);
        b32 InsideADBPlane = (Dot(ADBNormal, AO) * Dot(ADBNormal, AC) >= 0.0f);
        if(InsideABCPlane && InsideACDPlane && InsideADBPlane) return(true);

        // NOTE(georgy): 0 - ABC, 1 - ACD, 2 - ADB
        u32 TriangleIndex = 0;
        r32 ClosestDistSq = FLT_MAX;

        closest_voronoi_region_triangle_point ABC_VR, ACD_VR, ADB_VR;
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

        switch(TriangleIndex)
        {
            case 0:
            {
                DoSimplexTetrahedronSubset(ABC_VR, A, B, C, SimplexCount, Simplex, Dir);
            } break;

            case 1:
            {
                DoSimplexTetrahedronSubset(ACD_VR, A, C, D, SimplexCount, Simplex, Dir);
            } break;

            case 2:
            {
                DoSimplexTetrahedronSubset(ADB_VR, A, D, B, SimplexCount, Simplex, Dir);
            } break;
        }
    }

    Assert(*SimplexCount != 4);

    if(LengthSq(*Dir) <= Square(0.001f))
    {
        return(true);
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

internal b32
Intersect(game_object *ObjA, game_object *ObjB)
{
    vec3 AAxisX = vec3(Cos(Radians(ObjA->RigidBody.Orientation)), 0.0f, -Sin(Radians(ObjA->RigidBody.Orientation)));
    vec3 AAxisZ = vec3(Sin(Radians(ObjA->RigidBody.Orientation)), 0.0f, Cos(Radians(ObjA->RigidBody.Orientation)));
    vec3 AAxisY = vec3(0.0f, 1.0f, 0.0f);

    vec3 APolygonVertices[8] = 
    {
        ObjA->RigidBody.P + 0.5f*ObjA->Width*AAxisX + 0.5f*ObjA->Height*AAxisY + 0.5f*ObjA->Depth*AAxisZ,
        ObjA->RigidBody.P + 0.5f*ObjA->Width*AAxisX + 0.5f*ObjA->Height*AAxisY - 0.5f*ObjA->Depth*AAxisZ,
        ObjA->RigidBody.P - 0.5f*ObjA->Width*AAxisX + 0.5f*ObjA->Height*AAxisY - 0.5f*ObjA->Depth*AAxisZ,
        ObjA->RigidBody.P - 0.5f*ObjA->Width*AAxisX + 0.5f*ObjA->Height*AAxisY + 0.5f*ObjA->Depth*AAxisZ,
        ObjA->RigidBody.P - 0.5f*ObjA->Width*AAxisX - 0.5f*ObjA->Height*AAxisY + 0.5f*ObjA->Depth*AAxisZ,
        ObjA->RigidBody.P + 0.5f*ObjA->Width*AAxisX - 0.5f*ObjA->Height*AAxisY + 0.5f*ObjA->Depth*AAxisZ,
        ObjA->RigidBody.P + 0.5f*ObjA->Width*AAxisX - 0.5f*ObjA->Height*AAxisY - 0.5f*ObjA->Depth*AAxisZ,
        ObjA->RigidBody.P - 0.5f*ObjA->Width*AAxisX - 0.5f*ObjA->Height*AAxisY - 0.5f*ObjA->Depth*AAxisZ,
    };

#if 0
    vec3 BAxisX = vec3(Cos(Radians(ObjB->RigidBody.Orientation)), 0.0f, -Sin(Radians(ObjB->RigidBody.Orientation)));
    vec3 BAxisZ = vec3(Sin(Radians(ObjB->RigidBody.Orientation)), 0.0f, Cos(Radians(ObjB->RigidBody.Orientation)));
    vec3 BAxisY = vec3(0.0f, 1.0f, 0.0f);

    vec3 BPolygonVertices[8] = 
    {
        ObjB->RigidBody.P + 0.5f*ObjB->Width*BAxisX + 0.5f*ObjB->Height*BAxisY + 0.5f*ObjB->Depth*BAxisZ,
        ObjB->RigidBody.P + 0.5f*ObjB->Width*BAxisX + 0.5f*ObjB->Height*BAxisY - 0.5f*ObjB->Depth*BAxisZ,
        ObjB->RigidBody.P - 0.5f*ObjB->Width*BAxisX + 0.5f*ObjB->Height*BAxisY - 0.5f*ObjB->Depth*BAxisZ,
        ObjB->RigidBody.P - 0.5f*ObjB->Width*BAxisX + 0.5f*ObjB->Height*BAxisY + 0.5f*ObjB->Depth*BAxisZ,
        ObjB->RigidBody.P - 0.5f*ObjB->Width*BAxisX - 0.5f*ObjB->Height*BAxisY + 0.5f*ObjB->Depth*BAxisZ,
        ObjB->RigidBody.P + 0.5f*ObjB->Width*BAxisX - 0.5f*ObjB->Height*BAxisY + 0.5f*ObjB->Depth*BAxisZ,
        ObjB->RigidBody.P + 0.5f*ObjB->Width*BAxisX - 0.5f*ObjB->Height*BAxisY - 0.5f*ObjB->Depth*BAxisZ,
        ObjB->RigidBody.P - 0.5f*ObjB->Width*BAxisX - 0.5f*ObjB->Height*BAxisY - 0.5f*ObjB->Depth*BAxisZ,
    };

    vec3 MD[64];
    u32 MDCount = CalculateMinkowskiDiff(MD, ArrayCount(BPolygonVertices), BPolygonVertices, ArrayCount(APolygonVertices), APolygonVertices);
#else
    vec3 SphereP = ObjB->RigidBody.P;
    r32 SphereRadius = 1.0f;
#endif

    // vec3 S = Support(vec3(1.0f, 0.0f, 0.0f), MDCount, MD);
    vec3 S = Support(vec3(1.0f, 0.0f, 0.0f), SphereP, SphereRadius) - Support(-vec3(1.0f, 0.0f, 0.0f), ArrayCount(APolygonVertices), APolygonVertices);
    u32 SimplexCount = 1;
    vec3 Simplex[4] = {S, vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0)};
    vec3 Dir = -S;

    // NOTE(georgy): GJK pass
    b32 Result;
    while(true)
    {
        vec3 A = Support(Dir, SphereP, SphereRadius) - Support(-Dir, ArrayCount(APolygonVertices), APolygonVertices);
        // vec3 A = Support(Dir, MDCount, MD);
        if(Dot(Dir, A) < 0.0f)
        {
            Result = false;
            break;
        } 
        Simplex[SimplexCount++] = A;
        if(DoSimplex(&SimplexCount, Simplex, &Dir))
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
                    // Simplex[1] = Support(SearchDir, MDCount, MD);
                    Simplex[1] = Support(SearchDir, SphereP, SphereRadius) - Support(-SearchDir, ArrayCount(APolygonVertices), APolygonVertices);

                    if(LengthSq(Simplex[1] - Simplex[0]) >= ToleranceSq)
                    {
                        SimplexCount++;
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
                    Simplex[2] = Support(SearchDir, SphereP, SphereRadius) - Support(-SearchDir, ArrayCount(APolygonVertices), APolygonVertices);

                    if(Dot(Simplex[2], SearchDir) >= ToleranceSq)
                    {
                        SimplexCount++;
                        break;
                    }

                    SearchDir = Rot * SearchDir;
                }

                Assert(SimplexCount == 3);
            } 

            case 3:
            {
                vec3 SearchDir = Cross(Simplex[1] - Simplex[0], Simplex[2] - Simplex[0]);
                Simplex[3] = Support(SearchDir, SphereP, SphereRadius) - Support(-SearchDir, ArrayCount(APolygonVertices), APolygonVertices);
                
                if(Dot(Simplex[3], SearchDir) < ToleranceSq)
                {
                    SearchDir = -SearchDir;
                    Simplex[3] = Support(SearchDir, SphereP, SphereRadius) - Support(-SearchDir, ArrayCount(APolygonVertices), APolygonVertices);
                }

                SimplexCount++;
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

            C = Simplex[1];
            B = Simplex[2];
        }

        struct triangle
        {
            union
            {
                struct
                {
                    vec3 A, B, C;
                };
                vec3 E[3];
            };
        };

        // NOTE(georgy): Push first 4 triangles
        std::vector<triangle> Triangles;
        Triangles.push_back({A, B, C});
        Triangles.push_back({A, C, D}); 
        Triangles.push_back({A, D, B});
        Triangles.push_back({C, B, D});

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

            vec3 P = Support(ClosestNormal, SphereP, SphereRadius) - Support(-ClosestNormal, ArrayCount(APolygonVertices), APolygonVertices);
            r32 D = Dot(P, ClosestNormal);
            if(Absolute(D - ClosestDist) < Tolerance)
            {
                PenetrationVector = D*ClosestNormal;
                break;
            }
            else
            {
                struct edge
                {
                    vec3 Start;
                    vec3 End;
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
                            edge Edge = { Triangle.E[I1], Triangle.E[I0] };

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
                    Triangles.push_back({Edge.Start, Edge.End, P});
                }
            }
        }

        ObjA->RigidBody.P += PenetrationVector;
    }

    return(Result);
}

internal model
InitSphereMesh()
{
    u32 ParallelCount = 10;
    u32 MeridianCount = 10;
    r32 Radius = 1.0f;

    u32 VerticesCount = 0;
    vec3 Vertices[2048];

    u32 IndicesCount = 0;
    u32 Indices[2048];

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

        r32 CubeVertices[] = {
            // Back face
            -0.5f, -0.5f, -0.5f,
            0.5f,  0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f,  0.5f, -0.5f,
            -0.5f, -0.5f, -0.5f,
            -0.5f,  0.5f, -0.5f,
            // Front face
            -0.5f, -0.5f,  0.5f,
            0.5f, -0.5f,  0.5f,
            0.5f,  0.5f,  0.5f,
            0.5f,  0.5f,  0.5f,
            -0.5f,  0.5f,  0.5f,
            -0.5f, -0.5f,  0.5f,
            // Left face
            -0.5f,  0.5f,  0.5f,
            -0.5f,  0.5f, -0.5f,
            -0.5f, -0.5f, -0.5f,
            -0.5f, -0.5f, -0.5f,
            -0.5f, -0.5f,  0.5f,
            -0.5f,  0.5f,  0.5f,
            // Right face
            0.5f,  0.5f,  0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f,  0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f,  0.5f,  0.5f,
            0.5f, -0.5f,  0.5f,
            // Bottom face
            -0.5f, -0.5f, -0.5f,
            0.5f, -0.5f, -0.5f,
            0.5f, -0.5f,  0.5f,
            0.5f, -0.5f,  0.5f,
            -0.5f, -0.5f,  0.5f,
            -0.5f, -0.5f, -0.5f,
            // Top face
            -0.5f,  0.5f, -0.5f,
            0.5f,  0.5f , 0.5f,
            0.5f,  0.5f, -0.5f,
            0.5f,  0.5f,  0.5f,
            -0.5f,  0.5f, -0.5f,
            -0.5f,  0.5f,  0.5f,
        };

        InitModel(&GameState->Cube, sizeof(CubeVertices), CubeVertices, 36, 3*sizeof(r32));

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

        GameState->GameObjectCount = 2;
        GameState->Hero = GameState->GameObjects;
        GameState->Hero->Type = GameObject_Cube;
        GameState->Hero->Model = &GameState->Cube;
        GameState->Hero->Width = 1.0f;
        GameState->Hero->Height = 1.0f;
        GameState->Hero->Depth = 1.0f;
        GameState->Hero->RigidBody.P = vec3(0.0f, 0.0f, 0.0f);
        GameState->Hero->RigidBody.Mass = 6.0f;
        GameState->Hero->RigidBody.MomentOfInertia = (1.0f / 12.0f) * GameState->Hero->RigidBody.Mass * 
                                                      (Square(GameState->Hero->Width) + Square(GameState->Hero->Height));

        game_object *GameObject = &GameState->GameObjects[1];
        GameObject->Type = GameObject_Cube;
        GameObject->Model = &GameState->Sphere;
        GameObject->Width = 1.0f;
        GameObject->Height = 1.0f;
        GameObject->Depth = 1.0f;
        GameObject->RigidBody.P = vec3(0.0f, 0.0f, -3.0f);
        GameObject->RigidBody.Mass = 6.0f;
        GameObject->RigidBody.MomentOfInertia = (1.0f / 12.0f) * GameState->Hero->RigidBody.Mass * 
                                                      (Square(GameState->Hero->Width) + Square(GameState->Hero->Height));

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

    r32 CameraDistanceFromHero = 5.0f;
    r32 FloorDistanceFromHero = CameraDistanceFromHero * Cos(-PitchRadians);
    
    r32 XOffsetFromHero = FloorDistanceFromHero * Sin(HeadRadians);
    r32 YOffsetFromHero = CameraDistanceFromHero * Sin(-PitchRadians);
    r32 ZOffsetFromHero = FloorDistanceFromHero * Cos(HeadRadians);
    vec3 CameraOffsetFromHero = vec3(XOffsetFromHero, YOffsetFromHero, ZOffsetFromHero);

    vec3 CameraForward = Normalize(-CameraOffsetFromHero);
    vec3 Forward = Normalize(vec3(CameraForward.x, 0.0f, CameraForward.z));
    vec3 CameraRight = Normalize(Cross(Forward, vec3(0.0f, 1.0f, 1.0f)));
    r32 Theta = Degrees(ATan2(Forward.z, Forward.x)) - 90.0f;
    if(Input->MoveForward.EndedDown)
    {
        GameState->Hero->RigidBody.Orientation = Theta;
        GameState->Hero->RigidBody.P += 20.0f*dt*CameraForward;
    }
    if(Input->MoveBack.EndedDown)
    {
        GameState->Hero->RigidBody.Orientation = Theta + 180.0f;
        GameState->Hero->RigidBody.P += 20.0f*dt*-CameraForward;
    }
    if(Input->MoveRight.EndedDown)
    {
        GameState->Hero->RigidBody.Orientation = Theta - 90.0f;
        GameState->Hero->RigidBody.P += 20.0f*dt*CameraRight;
    }
    if(Input->MoveLeft.EndedDown)
    {
        GameState->Hero->RigidBody.Orientation = Theta + 90.0f;
        GameState->Hero->RigidBody.P += 20.0f*dt*-CameraRight;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        GameState->Hero->RigidBody.Orientation = Theta - 45.0f;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        GameState->Hero->RigidBody.Orientation = Theta + 45.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        GameState->Hero->RigidBody.Orientation = Theta - 135.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        GameState->Hero->RigidBody.Orientation = Theta + 135.0f;
    }

    Camera->P = GameState->Hero->RigidBody.P + CameraOffsetFromHero;
    Camera->Dir = CameraForward;

    mat4 Projection = Perspective(45.0f, (r32)WindowWidth/(r32)WindowHeight, 0.1f, 50.0f);
    mat4 View = Camera->GetRotationMatrix();

    Clear(RenderCommandBuffer, vec3(1.0f, 0.0f, 0.0f));

    PushShader(RenderCommandBuffer, GameState->Shader);
    PushMat4(RenderCommandBuffer, "Projection", &Projection);
    PushMat4(RenderCommandBuffer, "View", &View);
#if 1
    if(Intersect(&GameState->GameObjects[0], &GameState->GameObjects[1]))
    {
        PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 1.0f, 0.0f));
    }
    else
    {
        PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 0.0f, 1.0f));
    }

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < GameState->GameObjectCount;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;

        mat4 Model = Translation(GameObject->RigidBody.P) * Rotation(GameObject->RigidBody.Orientation, vec3(0.0f, 1.0f, 0.0f));
        PushMat4(RenderCommandBuffer, "Model", &Model);
        if(GameObjectIndex == 1)
        {
            DrawModelEBO(RenderCommandBuffer, GameObject->Model);
        }
        else
        {
            DrawModel(RenderCommandBuffer, GameObject->Model);
        }
    }

#else
    vec3 Triangle[3] = 
    {
        vec3(0.4f, -0.4f, 0.0f),
        vec3(-0.4f, -0.4f, 0.0f),
        vec3(0.0f, 0.4f, 0.0f),
    };
    PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 0.0f, 0.0f));

    mat4 Model = Identity();
    PushMat4(RenderCommandBuffer, "Model", &Model);
    for(u32 I0 = 0, I1 = ArrayCount(Triangle) - 1;
        I0 < ArrayCount(Triangle);
        I1 = I0, I0++)
    {
        DrawLine(RenderCommandBuffer, Triangle[I1], Triangle[I0]);
    }

    static vec3 Point = vec3(0.0f, 0.0f, 0.0f);
    if(Input->MoveForward.EndedDown)
    {
        // GameState->Hero->RigidBody.Orientation = Theta;
        Point += 10.0f*dt*vec3(0.0f, 1.0f, 0.0f);
    }
    if(Input->MoveBack.EndedDown)
    {
        Point += 10.0f*dt*vec3(0.0f, -1.0f, 0.0f);
    }
    if(Input->MoveRight.EndedDown)
    {
        Point += 10.0f*dt*vec3(1.0f, 0.0f, 0.0f);;
    }
    if(Input->MoveLeft.EndedDown)
    {
        Point += 10.0f*dt*vec3(-1.0f, 0.0f, 0.0f);;
    }

    Model = Translation(Point) * Scaling(0.05f);
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &GameState->Quad);

    closest_voronoi_region_triangle_point P_VR = ClosestPointInTriangleVR(Point, Triangle[0], Triangle[1], Triangle[2]);
    Model = Translation(P_VR.P) * Scaling(0.05f);
    PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 1.0f, 0.0f));
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &GameState->Quad);
#endif
}
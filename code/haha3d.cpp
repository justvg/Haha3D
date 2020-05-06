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

internal b32
VectorsAreEqual(vec3 A, vec3 B)
{
    b32 XIsEqual = Absolute(A.x - B.x) <= Epsilon;
    b32 YIsEqual = Absolute(A.y - B.y) <= Epsilon;
    b32 ZIsEqual = Absolute(A.z - B.z) <= Epsilon;

    b32 Result = XIsEqual && YIsEqual && ZIsEqual;
    return(Result);
}

struct collision_info
{
    vec3 ContactP;
    vec3 Normal;
    vec3 PenetrationVector;
};

internal b32
Intersect(collision_info *CollisionInfo, game_object *ObjA, game_object *ObjB)
{
    vec3 AAxisX = vec3(ObjA->RigidBody.Orientation.a11, ObjA->RigidBody.Orientation.a21, ObjA->RigidBody.Orientation.a31);
    vec3 AAxisY = vec3(ObjA->RigidBody.Orientation.a12, ObjA->RigidBody.Orientation.a22, ObjA->RigidBody.Orientation.a32);
    vec3 AAxisZ = vec3(ObjA->RigidBody.Orientation.a13, ObjA->RigidBody.Orientation.a23, ObjA->RigidBody.Orientation.a33);

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

#if 1
    vec3 BAxisX = vec3(ObjB->RigidBody.Orientation.a11, ObjB->RigidBody.Orientation.a21, ObjB->RigidBody.Orientation.a31);
    vec3 BAxisY = vec3(ObjB->RigidBody.Orientation.a12, ObjB->RigidBody.Orientation.a22, ObjB->RigidBody.Orientation.a32);
    vec3 BAxisZ = vec3(ObjB->RigidBody.Orientation.a13, ObjB->RigidBody.Orientation.a23, ObjB->RigidBody.Orientation.a33);

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
#else
    vec3 SphereP = ObjB->RigidBody.P;
    r32 SphereRadius = 1.0f;
#endif

    vec3 Sa = Support(vec3(1.0f, 0.0f, 0.0f), ArrayCount(BPolygonVertices), BPolygonVertices);
    vec3 Sb = Support(-vec3(1.0f, 0.0f, 0.0f), ArrayCount(APolygonVertices), APolygonVertices);
    vec3 S = Sa - Sb;
    u32 SimplexCount = 1;
    vec3 Simplex[4] = {S, vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0)};
    vec3 SimplexA[4] = {Sa, vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0)};
    vec3 SimplexB[4] = {Sb, vec3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0)};
    vec3 Dir = -S;

    Assert(VectorsAreEqual(Simplex[0], SimplexA[0] - SimplexB[0]));

    // NOTE(georgy): GJK pass
    b32 Result;
    while(true)
    {
        vec3 Aa = Support(Dir, ArrayCount(BPolygonVertices), BPolygonVertices);
        vec3 Ab = Support(-Dir, ArrayCount(APolygonVertices), APolygonVertices);
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
                    // Simplex[1] = Support(SearchDir, MDCount, MD);
                    SimplexA[1] = Support(SearchDir, ArrayCount(BPolygonVertices), BPolygonVertices);
                    SimplexB[1] = Support(-SearchDir, ArrayCount(APolygonVertices), APolygonVertices);
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
                    SimplexA[2] = Support(SearchDir, ArrayCount(BPolygonVertices), BPolygonVertices);
                    SimplexB[2] = Support(-SearchDir, ArrayCount(APolygonVertices), APolygonVertices);
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
                SimplexA[3] = Support(SearchDir, ArrayCount(BPolygonVertices), BPolygonVertices);
                SimplexB[3] = Support(-SearchDir, ArrayCount(APolygonVertices), APolygonVertices);
                Simplex[3] = SimplexA[3] - SimplexB[3];
                
                if(Dot(Simplex[3], SearchDir) < ToleranceSq)
                {
                    SearchDir = -SearchDir;
                    SimplexA[3] = Support(SearchDir, ArrayCount(BPolygonVertices), BPolygonVertices);
                    SimplexB[3] = Support(-SearchDir, ArrayCount(APolygonVertices), APolygonVertices);
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

            vec3 aP = Support(ClosestNormal, ArrayCount(BPolygonVertices), BPolygonVertices);
            vec3 bP = Support(-ClosestNormal, ArrayCount(APolygonVertices), APolygonVertices);
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

                vec3 TestProjO = u*Triangle.A + v*Triangle.B + w*Triangle.C;

                Assert(VectorsAreEqual(Triangle.A, Triangle.aA - Triangle.bA));
                Assert(VectorsAreEqual(Triangle.B, Triangle.aB - Triangle.bB));
                Assert(VectorsAreEqual(Triangle.C, Triangle.aC - Triangle.bC));

                // NOTE(georgy): 'A' object here is one that is A here "Support(A) - Support(B)""
                vec3 ContactPForA = u*Triangle.aA + v*Triangle.aB + w*Triangle.aC;
                vec3 ContactPForB = u*Triangle.bA + v*Triangle.bB + w*Triangle.bC;

                PenetrationVector = D*ClosestNormal;

                CollisionInfo->ContactP = ContactPForA;
                CollisionInfo->Normal = ClosestNormal;
                CollisionInfo->PenetrationVector = PenetrationVector;

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

        // ObjA->RigidBody.P += PenetrationVector;
    }

    return(Result);
}

internal model
InitSphereMesh(void)
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
        GameState->Hero->RigidBody.P = vec3(1.5f, 9.0f, 0.0f);
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
        GameState->Hero->RigidBody.CoeffOfRestitution = 0.3f;
        GameState->Hero->RigidBody.CoeffOfFriction = 0.1f;

        for(i32 ZOffset = -1;
            ZOffset <= 1;
            ZOffset++)
        {
            for(i32 YOffset = 0;
                YOffset <= 2;
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
                    GameObject->RigidBody.P = vec3(1.0f + 1.5f * XOffset, 2.0f*(YOffset) + 0.001f, 1.2f*ZOffset);
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
                    if(ZOffset == 0)
                        GameObject->RigidBody.Orientation = Rotation3x3(45.0f, vec3(1.0f, 0.3f, 0.0f));
                    else
                        GameObject->RigidBody.Orientation = Rotation3x3(-45.0f, vec3(1.0f, 0.0f, 4.0f));

                    // GameObject->RigidBody.Orientation = Identity3x3();
                    GameObject->RigidBody.InverseOrientation = Transpose3x3(GameObject->RigidBody.Orientation);
                    GameObject->RigidBody.GlobalInverseInertiaTensor = GameObject->RigidBody.Orientation *  
                                                                    GameObject->RigidBody.InverseInertiaTensor * 
                                                                    GameObject->RigidBody.InverseOrientation;
                    GameObject->RigidBody.CoeffOfRestitution = 0.3f;
                    GameObject->RigidBody.CoeffOfFriction = 0.1f;
                }
            }
        }

        game_object *GameObject = &GameState->GameObjects[GameState->GameObjectCount++];
        GameObject->Type = GameObject_Wall;
        GameObject->Model = &GameState->Cube;
        GameObject->Width = 1000.0f;
        GameObject->Height = 10.0f;
        GameObject->Depth = 1000.0f;
        GameObject->RigidBody.P = vec3(0.0f, -5.5f, 0.0f);
        GameObject->RigidBody.Mass = 0.0f;
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
        GameObject->RigidBody.CoeffOfRestitution = 0.3f;
        GameObject->RigidBody.CoeffOfFriction = 0.0f;   

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

    r32 CameraDistanceFromHero = 20.0f;
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
        GameState->Hero->RigidBody.ForceAccumulated += 10000.0f*dt*CameraForward;
    }
    if(Input->MoveBack.EndedDown)
    {
        
        GameState->Hero->RigidBody.ForceAccumulated += 10000.0f*dt*-CameraForward;
    }
    if(Input->MoveRight.EndedDown)
    {
        
        GameState->Hero->RigidBody.ForceAccumulated += 10000.0f*dt*CameraRight;
    }
    if(Input->MoveLeft.EndedDown)
    {
        
        GameState->Hero->RigidBody.ForceAccumulated += 10000.0f*dt*-CameraRight;
    }

    Camera->P = GameState->Hero->RigidBody.P + CameraOffsetFromHero;
    Camera->Dir = CameraForward;

    mat4 Projection = Perspective(45.0f, (r32)WindowWidth/(r32)WindowHeight, 0.1f, 50.0f);
    mat4 View = Camera->GetRotationMatrix();

    Clear(RenderCommandBuffer, vec3(1.0f, 0.0f, 0.0f));

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
            RigidBody->dP *= (1.0f - dt*1.0f);
            RigidBody->P += dt*RigidBody->dP;

            vec3 AngularAcceleration = RigidBody->GlobalInverseInertiaTensor*RigidBody->TorqueAccumulated;
            RigidBody->AngularMomentum += dt*RigidBody->TorqueAccumulated;
            RigidBody->AngularMomentum *= (1.0f - 0.25f*dt*1.0f);
            RigidBody->AngularVelocity = RigidBody->GlobalInverseInertiaTensor*RigidBody->AngularMomentum;
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

    struct collision_data
    {
        game_object *A;
		u32 AIndex;
        game_object *B;
		u32 BIndex;
        collision_info Manifold;
    };

    u32 CollisionCount = 0;
    collision_data Collisions[64];

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

            if(Intersect(&Collisions[CollisionCount].Manifold, GameObject, TestGameObject))
            {
                Collisions[CollisionCount].A = GameObject;
                Collisions[CollisionCount].B = TestGameObject;
                Collisions[CollisionCount].AIndex = GameObjectIndex;
				Collisions[CollisionCount].BIndex = TestObjectIndex;
                CollisionCount++;
            }
        }
    }

    for(u32 CollisionIndex = 0;
        CollisionIndex < CollisionCount;
        CollisionIndex++)
    {
        collision_data *Collision = Collisions + CollisionIndex;
        collision_info *Manifold = &Collision->Manifold;

        rigid_body *A = &Collision->A->RigidBody;
        rigid_body *B = &Collision->B->RigidBody;
        
        // NOTE(georgy): Linear projection
        r32 Slop = 0.01f;
        if(Length(Manifold->PenetrationVector) > Slop)
        {
            r32 Percent = 0.2f;
            vec3 PenetrationVector = Percent*Normalize(Manifold->PenetrationVector)*(Length(Manifold->PenetrationVector) - Slop);
            if(Collision->B->Type == GameObject_Wall)
            {
                A->P += PenetrationVector;
            }
            else
            {
                PenetrationVector *= (1.0f / (A->Mass + B->Mass));
                A->P += PenetrationVector * B->Mass;
                B->P -= PenetrationVector * A->Mass;
            }
        }

        vec3 N = Manifold->Normal;
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

        r32 CoeffOfRestitution = Min(A->CoeffOfRestitution, B->CoeffOfRestitution);
        r32 CoeffOfFriction = Min(A->CoeffOfFriction, B->CoeffOfFriction);

        vec3 aP = Manifold->ContactP - A->P;
        vec3 VelocityWithoutTranslation = Cross(A->AngularVelocity, aP);
        vec3 VelocityA = A->dP + VelocityWithoutTranslation;
        r32 OneOverMassA = 1.0f / A->Mass;
        b32 ImpulseWasUsed = false;

        r32 ImpulseMagnitude = 0.0f;
        r32 ImpulseMagnitudeTangent1 = 0.0f;
        r32 ImpulseMagnitudeTangent2 = 0.0f;
        if(B->Mass == 0.0f)
        {
            if(Dot(N, VelocityA) < 0.0f)
            {
                ImpulseWasUsed = true;

                r32 ImpulseNom = -(1.0f + CoeffOfRestitution)*Dot(VelocityA, N);
                r32 ImpulseDenom = OneOverMassA + Dot(N, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, N), aP)));
                ImpulseMagnitude = ImpulseNom / ImpulseDenom;

                r32 ImpulseNomTangent1 = -(1.0f + CoeffOfRestitution)*Dot(VelocityA, Tangent1);
                r32 ImpulseDenomTangent1 = OneOverMassA + Dot(Tangent1, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, Tangent1), aP)));   
                ImpulseMagnitudeTangent1 = ImpulseNomTangent1 / ImpulseDenomTangent1;
                if(Absolute(ImpulseMagnitudeTangent1) > CoeffOfFriction*ImpulseMagnitude)
                {
                    ImpulseMagnitudeTangent1 = Sign(ImpulseMagnitudeTangent1)*CoeffOfFriction*ImpulseMagnitude;
                }

                r32 ImpulseNomTangent2 = -(1.0f + CoeffOfRestitution)*Dot(VelocityA, Tangent2);
                r32 ImpulseDenomTangent2 = OneOverMassA + Dot(Tangent2, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, Tangent2), aP)));   
                ImpulseMagnitudeTangent2 = ImpulseNomTangent2 / ImpulseDenomTangent2;
                if(Absolute(ImpulseMagnitudeTangent2) > CoeffOfFriction*ImpulseMagnitude)
                {
                    ImpulseMagnitudeTangent2 = Sign(ImpulseMagnitudeTangent2)*CoeffOfFriction*ImpulseMagnitude;
                }
            }
        }
        else
        {
            vec3 bP = Manifold->ContactP - B->P;
            vec3 VelocityB = B->dP + Cross(B->AngularVelocity, bP);
            r32 OneOverMassB = 1.0f / B->Mass;

            vec3 RelativeVelocity = VelocityA - VelocityB;
            if(Dot(N, RelativeVelocity) < 0.0f)
            {
                ImpulseWasUsed = true;

                r32 ImpulseNom = -(1.0f + CoeffOfRestitution)*Dot(RelativeVelocity, N);
                r32 ImpulseDenom = (OneOverMassA + OneOverMassB)+
                                   Dot(N, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, N), aP))) + 
                                   Dot(N, (Cross(B->GlobalInverseInertiaTensor*Cross(bP, N), bP)));
                ImpulseMagnitude = ImpulseNom / ImpulseDenom;

                r32 ImpulseNomTangent1 = -(1.0f + CoeffOfRestitution)*Dot(VelocityA - VelocityB, Tangent1);
                r32 ImpulseDenomTangent1 = (OneOverMassA + OneOverMassB)+
                                           Dot(Tangent1, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, Tangent1), aP))) + 
                                           Dot(Tangent1, (Cross(B->GlobalInverseInertiaTensor*Cross(bP, Tangent1), bP)));
                ImpulseMagnitudeTangent1 = ImpulseNomTangent1 / ImpulseDenomTangent1;
                if(Absolute(ImpulseMagnitudeTangent1) > CoeffOfFriction*ImpulseMagnitude)
                {
                    ImpulseMagnitudeTangent1 = Sign(ImpulseMagnitudeTangent1)*CoeffOfFriction*ImpulseMagnitude;
                }

                r32 ImpulseNomTangent2 = -(1.0f + CoeffOfRestitution)*Dot(VelocityA - VelocityB, Tangent2);
                r32 ImpulseDenomTangent2 = (OneOverMassA + OneOverMassB)+
                                           Dot(Tangent2, (Cross(A->GlobalInverseInertiaTensor*Cross(aP, Tangent2), aP))) + 
                                           Dot(Tangent2, (Cross(B->GlobalInverseInertiaTensor*Cross(bP, Tangent2), bP)));
                ImpulseMagnitudeTangent2 = ImpulseNomTangent2 / ImpulseDenomTangent2;
                if(Absolute(ImpulseMagnitudeTangent2) > CoeffOfFriction*ImpulseMagnitude)
                {
                    ImpulseMagnitudeTangent2 = Sign(ImpulseMagnitudeTangent2)*CoeffOfFriction*ImpulseMagnitude;
                }

                B->dP += (-ImpulseMagnitude * OneOverMassB)*N + 
                         (-ImpulseMagnitudeTangent1 * OneOverMassB)*Tangent1 +
                         (-ImpulseMagnitudeTangent2 * OneOverMassB)*Tangent2;
                B->AngularMomentum += Cross(bP, -ImpulseMagnitude*N) +
                                      Cross(bP, -ImpulseMagnitudeTangent1*Tangent1) +
                                      Cross(bP, -ImpulseMagnitudeTangent2*Tangent2);
                B->AngularVelocity = B->GlobalInverseInertiaTensor*B->AngularMomentum;
            }
        }
        
        A->dP += (ImpulseMagnitude * OneOverMassA)*N + 
                 (ImpulseMagnitudeTangent1 * OneOverMassA)*Tangent1 + 
                 (ImpulseMagnitudeTangent2 * OneOverMassA)*Tangent2;
        A->AngularMomentum += Cross(aP, ImpulseMagnitude*N) +
                              Cross(aP, ImpulseMagnitudeTangent1*Tangent1) + 
                              Cross(aP, ImpulseMagnitudeTangent2*Tangent2);
        A->AngularVelocity = A->GlobalInverseInertiaTensor*A->AngularMomentum;

        // NOTE(georgy): Check normal relative velocity after impulse
        if(B->Mass != 0.0f)
        {
            if(ImpulseWasUsed)
            {
                vec3 ATestAngularVelocity = A->GlobalInverseInertiaTensor*A->AngularMomentum;
                vec3 ATestVelocityWithoutTranslation = Cross(ATestAngularVelocity, aP);
                vec3 TestVelocityA = A->dP + ATestVelocityWithoutTranslation;

                vec3 bP = Manifold->ContactP - B->P;
                vec3 BTestAngularVelocity = B->GlobalInverseInertiaTensor*B->AngularMomentum;
                vec3 BTestVelocityWithoutTranslation = Cross(BTestAngularVelocity, bP);
                vec3 TestVelocityB = B->dP + BTestVelocityWithoutTranslation;

                vec3 TestRelativeVelocity = TestVelocityA - TestVelocityB;
                r32 NormalRelativeVelocity = Dot(N, TestRelativeVelocity);
                Assert(NormalRelativeVelocity >= 0.0f);
            }
        }
        else
        {
            if(ImpulseWasUsed)
            {
                vec3 ATestAngularVelocity = A->GlobalInverseInertiaTensor*A->AngularMomentum;
                vec3 ATestVelocityWithoutTranslation = Cross(ATestAngularVelocity, aP);
                vec3 TestVelocityA = A->dP + ATestVelocityWithoutTranslation;

                r32 NormalRelativeVelocity = Dot(N, TestVelocityA);
                Assert(NormalRelativeVelocity >= 0.0f);
            }
        }
    }

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < GameState->GameObjectCount;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;

        mat4 Model = Translation(GameObject->RigidBody.P) * 
                     Mat4(GameObject->RigidBody.Orientation) * 
                     Scaling(vec3(GameObject->Width, GameObject->Height, GameObject->Depth));
        PushMat4(RenderCommandBuffer, "Model", &Model);
        if(GameObjectIndex == (GameState->GameObjectCount - 1))
        {
            PushVec3(RenderCommandBuffer, "Color", vec3(1.0f, 1.0f, 0.0f));
        }
        else
        {
            PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 0.0f, 1.0f));
        }
        DrawModel(RenderCommandBuffer, GameObject->Model);
    }
}
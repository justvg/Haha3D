#include "haha3d.h"
#include "haha3d_math.cpp"
#include "haha3d_render_command_buffer.cpp"

mat4 
camera::GetRotationMatrix(void)
{ 
    mat4 Result = LookAt(P, P + Dir); 

    return(Result);
} 

internal vec2
GetPointOfContact(rigid_body *RigidBody, r32 t, vec2 DeltaP, 
                  r32 HalfWidth, vec2 OrientationX, r32 HalfHeight, vec2 OrientationY,
                  plane Plane, r32 Radius)
{
    vec2 RigidBodyCorners[4];
    RigidBodyCorners[0] = RigidBody->P + t*DeltaP + 
        (HalfWidth*OrientationX + HalfHeight*OrientationY);
    RigidBodyCorners[1] = RigidBody->P + t*DeltaP + 
        (HalfWidth*OrientationX - HalfHeight*OrientationY);
    RigidBodyCorners[2] = RigidBody->P + t*DeltaP + 
        (-HalfWidth*OrientationX - HalfHeight*OrientationY);
    RigidBodyCorners[3] = RigidBody->P + t*DeltaP + 
        (-HalfWidth*OrientationX + HalfHeight*OrientationY);

    u32 ClosestCount = 0;
    u32 ClosestIndecies[4];

    u32 TheClosestToPlaneIndex = 0;
    r32 ClosestDistance = FLT_MAX;

    for(u32 Corner = 0;
        Corner < 4;
        Corner++)
    {
        r32 Distance = Absolute(Dot(RigidBodyCorners[Corner], Plane.N) - Plane.D);
        if(Distance > 0.0f)
        {
            r32 Diff = ClosestDistance - Distance;
            if(Absolute(Diff) <= Radius*0.1f)
            {
                ClosestIndecies[ClosestCount++] = Corner;
            }
            else if(Diff > 0.0f)
            {
                TheClosestToPlaneIndex = Corner;
                ClosestDistance = Distance;

                ClosestCount = 0;
                ClosestIndecies[ClosestCount++] = Corner;
            }
        }
    }

    // Assert(ClosestCount <= 2);

    vec2 PointOfContact;
    if(ClosestCount > 2)
    {
        PointOfContact = Lerp(RigidBodyCorners[ClosestIndecies[0]], RigidBodyCorners[ClosestIndecies[1]], 0.5f);
    }
    else
    {
        PointOfContact = RigidBodyCorners[TheClosestToPlaneIndex];
    }

    return(PointOfContact);
}

internal b32
Contains(vec2 *Vertices, u32 Count, vec2 Check)
{
    b32 Result = false;

    for(u32 Index = 0;
        Index < Count;
        Index++)
    {
        vec2 Vertex = Vertices[Index];
        b32 XIsEqual = Absolute(Vertex.x - Check.x) <= Epsilon;
        b32 YIsEqual = Absolute(Vertex.y - Check.y) <= Epsilon;
        if(XIsEqual && YIsEqual)
        {
            Result = true;
            break;
        }
    }

    return(Result);
}

internal u32
CalculateMinkowskiDiff(vec2 *Diff, u32 ACount, vec2 *AVertices, u32 BCount, vec2 *BVertices)
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
            vec2 NewVertex = AVertices[AVertex] - BVertices[BVertex];
            if(!Contains(Diff, VertexCount, NewVertex))
            {
                Diff[VertexCount++] = NewVertex;
            }
        }
    }

    return(VertexCount);
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
        b32 ZIsEqual = Absolute(Vertex.y - Check.y) <= Epsilon;
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

internal b32
DoSimplex(u32 *SimplexCount, vec3 *Simplex, vec3 *D)
{
    b32 Result = false;

    if(*SimplexCount == 2)
    {
        vec3 B = Simplex[0];
        vec3 A = Simplex[1];

        vec3 AB = B - A;
        vec3 AO = -A;

        if(Dot(AB, AO) > 0.0f)
        {
            *D = Cross(Cross(AB, AO), AB);
        }
        else
        {
            *SimplexCount = 1;
            Simplex[0] = A;
            *D = AO;
        }
    }
    else if(*SimplexCount == 3)
    {
        vec3 C = Simplex[0];
        vec3 B = Simplex[1]; 
        vec3 A = Simplex[2];

        vec3 AB = B - A;
        vec3 AC = C - A;
        vec3 AO = -A;

        vec3 ABCNormal = Cross(AB, AC);
        vec3 ACNormal = Cross(ABCNormal, AC);
        vec3 ABNormal = Cross(AB, ABCNormal);

        if(Dot(ACNormal, AO) > 0.0f)
        {
            if(Dot(AC, AO) > 0.0f)
            {
                *SimplexCount = 2;
                Simplex[0] = C;
                Simplex[1] = A;
                *D = Cross(Cross(AC, AO), AC);
            }
            else
            {
                if(Dot(AB, AO) > 0.0f)
                {
                    *SimplexCount = 2;
                    Simplex[0] = B;
                    Simplex[1] = A;
                    *D = Cross(Cross(AB, AO), AB);
                }
                else
                {
                    *SimplexCount = 1;
                    Simplex[0] = A;
                    *D = AO;
                }
            }
        }
        else
        {
            if(Dot(ABNormal, AO) > 0.0f)
            {
                if(Dot(AB, AO) > 0.0f)
                {
                    *SimplexCount = 2;
                    Simplex[0] = B;
                    Simplex[1] = A;
                    *D = Cross(Cross(AB, AO), AB);
                }
                else
                {
                    *SimplexCount = 1;
                    Simplex[0] = A;
                    *D = AO;
                }
            }
            else
            {
                Result = true;
            }
        }
    }

    return(Result);
}

internal b32
Intersect(game_object *A, game_object *B)
{
    vec2 AAxisX = vec2(Cos(Radians(A->RigidBody.Orientation)), Sin(Radians(A->RigidBody.Orientation)));
    vec2 AAxisY = Perp(AAxisX);
    // TODO(georgy): This function works for any convex polygon!
    // NOTE(georgy): CCW order
#if 0
    vec3 APolygonVertices[4] = 
    {
        vec3(A->RigidBody.P + 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P + 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY, 0.0f),
    };
#else
    vec3 APolygonVertices[8] = 
    {
        vec3(A->RigidBody.P + 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P + 0.25f*A->Width*AAxisX + 0.75f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.25f*A->Width*AAxisX + 0.75f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.25f*A->Width*AAxisX - 0.75f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P + 0.25f*A->Width*AAxisX - 0.75f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P + 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY, 0.0f),
    };
#endif
    vec2 BAxisX = vec2(Cos(Radians(B->RigidBody.Orientation)), Sin(Radians(B->RigidBody.Orientation)));
    vec2 BAxisY = Perp(BAxisX);
    // TODO(georgy): This function works for any convex polygon!
    // NOTE(georgy): CCW order
    vec3 BPolygonVertices[4] = 
    {
        vec3(B->RigidBody.P + 0.5f*B->Width*BAxisX + 0.5f*B->Height*BAxisY, 0.0f),
        vec3(B->RigidBody.P - 0.5f*B->Width*BAxisX + 0.5f*B->Height*BAxisY, 0.0f),
        vec3(B->RigidBody.P - 0.5f*B->Width*BAxisX - 0.5f*B->Height*BAxisY, 0.0f),
        vec3(B->RigidBody.P + 0.5f*B->Width*BAxisX - 0.5f*B->Height*BAxisY, 0.0f),
    };

    vec3 MD[32];
    u32 MDCount = CalculateMinkowskiDiff(MD, ArrayCount(APolygonVertices), APolygonVertices, ArrayCount(BPolygonVertices), BPolygonVertices);

    vec3 S = Support(vec3(1.0f, 0.0f, 0.0f), MDCount, MD);
    u32 SimplexCount = 1;
    vec3 Simplex[3] = {S, vec3(0, 0, 0), vec3(0, 0, 0)};
    vec3 D = -S;

    while(true)
    {
        vec3 A = Support(D, MDCount, MD);
        if(Dot(D, A) < 0.0f) return(false); // NOTE(georgy): No intersection
        Simplex[SimplexCount++] = A;
        if(DoSimplex(&SimplexCount, Simplex, &D)) return(true); // NOTE(georgy): Intersection
    }
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
            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 0.0f, 1.0f, 0.0f,
            -0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f
            // Front face
            -0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            -0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            -0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 0.0f,
            // Left face
            -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f, -0.5f, -0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f, -0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
            -0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 0.0f,
            // Right face
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
            0.5f,  0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 1.0f,
            0.5f,  0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
            0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 1.0f,
            // Bottom face
            -0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f, -0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            0.5f, -0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f, -0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            // Top face
            -0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f,  0.5f , 0.5f, 0.5f, 0.5f, 0.5f, 
            0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            0.5f,  0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f,  0.5f, -0.5f, 0.5f, 0.5f, 0.5f,
            -0.5f,  0.5f,  0.5f, 0.5f, 0.5f, 0.5f,
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

        GameState->CubeOrientation = Identity();

        GameState->GameObjectCount = 1;
        GameState->Hero = GameState->GameObjects;
        GameState->Hero->Type = GameObject_Rectangle;
        GameState->Hero->Model = &GameState->Quad;
        GameState->Hero->Width = 0.1f;
        GameState->Hero->Height = 0.05f;
        GameState->Hero->RigidBody.P = vec2(0.0f, -0.33f);
        GameState->Hero->RigidBody.Mass = 3.0f;
        GameState->Hero->RigidBody.MomentOfInertia = (1.0f / 12.0f) * GameState->Hero->RigidBody.Mass * 
                                                      (Square(GameState->Hero->Width) + Square(GameState->Hero->Height));

        for(i32 YOffset = -1;
            YOffset <= 1;
            YOffset++)
        {
            for(i32 XOffset = -3;
                XOffset <= 3;
                XOffset++)
            {
                u32 GameObjectIndex = GameState->GameObjectCount;
                GameState->GameObjects[GameObjectIndex].Type = GameObject_Rectangle;
                GameState->GameObjects[GameObjectIndex].Model = &GameState->Quad;
                GameState->GameObjects[GameObjectIndex].RigidBody.P = vec2(XOffset*0.15f, YOffset*0.1f);
                GameState->GameObjects[GameObjectIndex].Width = 0.1f;
                GameState->GameObjects[GameObjectIndex].Height = 0.05f;
                GameState->GameObjects[GameObjectIndex].RigidBody.Mass = 3.0f;
                GameState->GameObjects[GameObjectIndex].RigidBody.MomentOfInertia = 
                        (1.0f / 12.0f) * GameState->GameObjects[GameObjectIndex].RigidBody.Mass * 
                        (Square(GameState->GameObjects[GameObjectIndex].Width) + Square(GameState->GameObjects[GameObjectIndex].Height));

                GameState->GameObjectCount++;
            }
        }

        GameState->IsInitialized = true;
    }

    r32 dt = Input->dtForFrame;

    camera *Camera = &GameState->Camera;
#if 0
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
    vec3 CameraRight = Normalize(Cross(CameraForward, vec3(0.0f, 1.0f, 1.0f)));
    r32 Theta = Degrees(ATan2(CameraForward.z, CameraForward.x)) - 90.0f;
    if(Input->MoveForward.EndedDown)
    {
        GameState->HeroRotation = Theta;
        GameState->HeroP += 20.0f*dt*CameraForward;
    }
    if(Input->MoveBack.EndedDown)
    {
        GameState->HeroRotation = Theta + 180.0f;
        GameState->HeroP += 20.0f*dt*-CameraForward;
    }
    if(Input->MoveRight.EndedDown)
    {
        GameState->HeroRotation = Theta - 90.0f;
        GameState->HeroP += 20.0f*dt*CameraRight;
    }
    if(Input->MoveLeft.EndedDown)
    {
        GameState->HeroRotation = Theta + 90.0f;
        GameState->HeroP += 20.0f*dt*-CameraRight;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        GameState->HeroRotation = Theta - 45.0f;
    }
    if(Input->MoveForward.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        GameState->HeroRotation = Theta + 45.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveRight.EndedDown)
    {
        GameState->HeroRotation = Theta - 135.0f;
    }
    if(Input->MoveBack.EndedDown && 
        Input->MoveLeft.EndedDown)
    {
        GameState->HeroRotation = Theta + 135.0f;
    }

    Camera->P = GameState->HeroP + CameraOffsetFromHero;
    Camera->Dir = CameraForward;
#else

    GameState->Hero->RigidBody.ForceAccumulated = vec2(0.0f, 0.0f);
    GameState->Hero->RigidBody.TorqueAccumulated = 0.0f;

    static b32 ForceActsOnCM = true;
    if(WasDown(&Input->MouseLeft))
    {
        ForceActsOnCM = !ForceActsOnCM;
    }

    static b32 GravityEnabled = false;
    if(WasDown(&Input->MouseRight))
    {
        GravityEnabled = !GravityEnabled;
    }

    rigid_body *HeroRigidBody = &GameState->Hero->RigidBody;
    r32 HeroHalfWidth = 0.5f*GameState->Hero->Width;
    r32 HeroHalfHeight = 0.5f*GameState->Hero->Height;
    vec2 HeroOrientationAxisX = vec2(Cos(Radians(HeroRigidBody->Orientation)), Sin(Radians(HeroRigidBody->Orientation)));
    vec2 HeroOrientationAxisY = Perp(HeroOrientationAxisX);

    // NOTE(georgy): At the moment these forces are acting only on the CM of the hero
    if(Input->MoveForward.EndedDown)
    {
        vec2 Force = 200.0f*vec2(0.0f, 1.0f);
        HeroRigidBody->ForceAccumulated += Force;
        if(!ForceActsOnCM)
        {
            HeroRigidBody->TorqueAccumulated += Cross2D(HeroHalfHeight*HeroOrientationAxisY + 0.2f*HeroHalfWidth*HeroOrientationAxisX, Force);
        }
    }
    if(Input->MoveBack.EndedDown)
    {
        vec2 Force = -200.0f*vec2(0.0f, 1.0f);
        HeroRigidBody->ForceAccumulated += Force;
        if(!ForceActsOnCM)
        {
            HeroRigidBody->TorqueAccumulated += Cross2D(-HeroHalfHeight*HeroOrientationAxisY - 0.2f*HeroHalfWidth*HeroOrientationAxisX, Force);
        }
    }
    if(Input->MoveRight.EndedDown)
    {
        vec2 Force = 200.0f*vec2(1.0f, 0.0f);
        HeroRigidBody->ForceAccumulated += Force;
        if(!ForceActsOnCM)
        {
            HeroRigidBody->TorqueAccumulated += Cross2D(-HeroHalfWidth*HeroOrientationAxisX + 0.2f*HeroHalfHeight*HeroOrientationAxisY, Force);
        }
    }
    if(Input->MoveLeft.EndedDown)
    {
        vec2 Force = -200.0f*vec2(1.0f, 0.0f);
        HeroRigidBody->ForceAccumulated += Force;
        if(!ForceActsOnCM)
        {
            HeroRigidBody->TorqueAccumulated += Cross2D(HeroHalfWidth*HeroOrientationAxisX - 0.2f*HeroHalfHeight*HeroOrientationAxisY, Force);
        }
    }

    Camera->P = vec3(0.0f, 0.0f, 1.0f);
    Camera->Dir = vec3(0.0f, 0.0f, -1.0f);

#endif
    mat4 Projection = Perspective(45.0f, (r32)WindowWidth/(r32)WindowHeight, 0.1f, 50.0f);
    mat4 View = Camera->GetRotationMatrix();

    Clear(RenderCommandBuffer, vec3(1.0f, 0.0f, 0.0f));

    PushShader(RenderCommandBuffer, GameState->Shader);
    PushMat4(RenderCommandBuffer, "Projection", &Projection);
    PushMat4(RenderCommandBuffer, "View", &View);
#if 0
    mat4 Model = Translation(GameState->HeroP) * Rotation(GameState->HeroRotation, vec3(0.0f, 1.0f, 0.0f));
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &GameState->Cube);

    Model = Translation(vec3(0.0f, -2.0f, 0.0f)) * Rotation(-90.0f, vec3(1.0f, 0.0f, 0.0f)) * Scaling(6.0f);
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &GameState->Plane);

    GameState->CubeOrientation = Rotation(dt*180.0f, vec3(0.0f, 1.0f, 0.0f)) * GameState->CubeOrientation;
    Model = Translation(vec3(0.0f, -1.5f, 0.0f)) * GameState->CubeOrientation;
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawModel(RenderCommandBuffer, &GameState->Cube);
#else
#if 0

    plane Planes[4] = 
    {
        {vec2(0.0f, 1.0f), -0.4f},
        {vec2(-1.0f, 0.0f), -0.8f},
        {vec2(0.0f, -1.0f), -0.4f},
        {vec2(1.0f, 0.0f), -0.8f},
    };

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < GameState->GameObjectCount;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;
        rigid_body *RigidBody = &GameObject->RigidBody;

        RigidBody->dP += dt*((RigidBody->ForceAccumulated*(1.0f / RigidBody->Mass)) + (GravityEnabled ? vec2(0.0f, -9.8f) : vec2(0.0f, 0.0f)));
        RigidBody->AngularSpeed += dt*(RigidBody->TorqueAccumulated / RigidBody->MomentOfInertia);
        RigidBody->Orientation += dt*Degrees(RigidBody->AngularSpeed);
        vec2 BodyOrientationAxisX = vec2(Cos(Radians(RigidBody->Orientation)), Sin(Radians(RigidBody->Orientation)));
        vec2 BodyOrientationAxisY = Perp(BodyOrientationAxisX);
        r32 BodyHalfWidth = 0.5f*GameObject->Width;
        r32 BodyHalfHeight = 0.5f*GameObject->Height;

        vec2 GameObjectDeltaP = dt*RigidBody->dP;
        r32 dtRemaining = dt;
        for(u32 Iteration = 0;
            Iteration < 100;
            Iteration++)
        {
            r32 GameObjectDeltaPLength = Length(GameObjectDeltaP);
            if((GameObjectDeltaPLength > 0.0f) && (dtRemaining > 0.0f))
            {
                r32 t = 1.0f;
                vec2 PointOfContact;
                vec2 CollisionNormal;
                b32 Collision = false;
                game_object *CollisionGameObject = 0;

                vec2 DesiredP = RigidBody->P + GameObjectDeltaP;

                for(u32 PlaneIndex = 0;
                    PlaneIndex < ArrayCount(Planes);
                    PlaneIndex++)
                {
                    plane *Plane = Planes + PlaneIndex;
                    b32 CollisionWithPlane = false;

                    r32 Radius = (BodyHalfWidth*Absolute(Dot(BodyOrientationAxisX, Plane->N)) + BodyHalfHeight*Absolute(Dot(BodyOrientationAxisY, Plane->N)));
                    r32 DistFromHeroCenterToPlane = Dot(Plane->N, RigidBody->P) - Plane->D;

                    u32 MaxDepth = 100;
                    while((Absolute(DistFromHeroCenterToPlane) < Radius) && (MaxDepth > 0))
                    {
                        // NOTE(georgy): Penetration
                        RigidBody->P += (1.1f*Radius - Absolute(DistFromHeroCenterToPlane))*Plane->N;
                        DistFromHeroCenterToPlane = Dot(Plane->N, RigidBody->P) - Plane->D;

                        MaxDepth--;
                    }   

                    Assert(Absolute(DistFromHeroCenterToPlane) >= Radius);

                    {
                        r32 Denom = Dot(Plane->N, GameObjectDeltaP);
                        if((Denom * DistFromHeroCenterToPlane) >= 0.0f)
                        {
                            // NOTE(georgy): Moving parallel to or away from the plane
                            CollisionWithPlane = false;
                        }
                        else
                        {
                            r32 PlaneDisplace = (DistFromHeroCenterToPlane > 0.0f) ? Radius : -Radius;
                            r32 NewT = (PlaneDisplace - DistFromHeroCenterToPlane) / Denom;
                            if((NewT <= t) && (NewT >= 0.0f))
                            {
                                t = NewT - 0.1f;
                                if(t < 0.0f) t = NewT;
                                CollisionWithPlane = true;

                                PointOfContact = GetPointOfContact(RigidBody, t, GameObjectDeltaP, 
                                                                   BodyHalfWidth, BodyOrientationAxisX, BodyHalfHeight, BodyOrientationAxisY,
                                                                   *Plane, Radius);
                                CollisionNormal = Plane->N;
                            }
                        }
                    }

                    if(CollisionWithPlane) 
                    {
                        Collision = true;
                    }
                }

                for(u32 TestGameObjectIndex = 0;
                    TestGameObjectIndex < GameState->GameObjectCount;
                    TestGameObjectIndex++)
                {
                    if(TestGameObjectIndex != GameObjectIndex)
                    {
                        game_object *TestGameObject = GameState->GameObjects + TestGameObjectIndex;
                        
                        r32 NewT = FLT_MAX;
                        vec2 NewPointOfContact = {};
                        vec2 NewCollisionNormal = {};
                        if(TestIntersection(GameObject, GameObjectDeltaP, TestGameObject, 
                                            &NewT, NewPointOfContact, NewCollisionNormal))
                        {
                            if((NewT <= t) && (NewT >= 0.0f))
                            {
                                t = NewT - 0.1f;
                                PointOfContact = NewPointOfContact;
                                CollisionNormal = NewCollisionNormal;

                                Assert((Absolute(PointOfContact.x) > Epsilon) || (Absolute(PointOfContact.y) > Epsilon));
                                Assert((Absolute(CollisionNormal.x) > Epsilon) || (Absolute(CollisionNormal.y) > Epsilon));

                                Collision = true;
                                CollisionGameObject = TestGameObject;
                            }
                        }
                    }
                }

                dtRemaining -= t*dtRemaining;
                RigidBody->P += t*GameObjectDeltaP;
                GameObjectDeltaP = DesiredP - RigidBody->P;
                if(Collision)
                {
                    r32 CoeffOfRestitution = 0.25f;
                    
                    r32 OneOverMassA = 1.0f / RigidBody->Mass;
                    r32 OneOverMomentOfInertiaA = 1.0f / RigidBody->MomentOfInertia;
                    vec2 PointRotationDirA = Perp(PointOfContact - RigidBody->P);
                    vec2 VelocityA = RigidBody->dP + RigidBody->AngularSpeed*PointRotationDirA;

                    r32 ImpulseMagnitude = 0.0f;
                    if(!CollisionGameObject)
                    {
                        r32 ImpulseNom = -(1.0f + CoeffOfRestitution)*Dot(VelocityA, CollisionNormal);
                        r32 ImpulseDenom = Dot(CollisionNormal, CollisionNormal)*OneOverMassA + 
                                                OneOverMomentOfInertiaA*Square(Dot(PointRotationDirA, CollisionNormal));

                        ImpulseMagnitude = ImpulseNom / ImpulseDenom;
                    }
                    else
                    {
                        rigid_body *RigidBodyB = &CollisionGameObject->RigidBody;
                        r32 OneOverMassB = 1.0f / RigidBodyB->Mass;
                        r32 OneOverMomentOfInertiaB = 1.0f / RigidBodyB->MomentOfInertia;
                        vec2 PointRotationDirB = Perp(PointOfContact - RigidBodyB->P);
                        vec2 VelocityB = RigidBodyB->dP + RigidBodyB->AngularSpeed*PointRotationDirB;

                        r32 ImpulseNom = -(1.0f + CoeffOfRestitution)*Dot(VelocityA - VelocityB, CollisionNormal);
                        r32 ImpulseDenom = Dot(CollisionNormal, CollisionNormal)*(OneOverMassA + OneOverMassB) + 
                                            OneOverMomentOfInertiaA*Square(Dot(PointRotationDirA, CollisionNormal)) + 
                                            OneOverMomentOfInertiaB*Square(Dot(PointRotationDirB, CollisionNormal));

                        ImpulseMagnitude = ImpulseNom / ImpulseDenom;

                        RigidBodyB->dP += (-ImpulseMagnitude * OneOverMassB)*CollisionNormal;
                        RigidBodyB->AngularSpeed += Dot(PointRotationDirB, -ImpulseMagnitude*CollisionNormal) * OneOverMomentOfInertiaB;
                    }

                    RigidBody->dP += (ImpulseMagnitude * OneOverMassA)*CollisionNormal;
                    RigidBody->AngularSpeed += Dot(PointRotationDirA, ImpulseMagnitude*CollisionNormal) * OneOverMomentOfInertiaA;

                    GameObjectDeltaP = dtRemaining*RigidBody->dP;
                }
            }
            else
            {
                break;
            }
        }
    }

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < GameState->GameObjectCount;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;
        mat4 Model = Translation(vec3(GameObject->RigidBody.P, 0.0f)) * 
                     Rotation(GameObject->RigidBody.Orientation, vec3(0.0f, 0.0f, 1.0f)) * 
                     Scaling(vec3(GameObject->Width, GameObject->Height, 1.0f));
        PushMat4(RenderCommandBuffer, "Model", &Model);
        DrawModel(RenderCommandBuffer, GameObject->Model);
    }

    mat4 Model = Identity();
    PushMat4(RenderCommandBuffer, "Model", &Model);
    DrawLine(RenderCommandBuffer, vec3(-1.0f, Planes[0].D*Planes[0].N.y, 0.0f), vec3(1.0f, Planes[0].D*Planes[0].N.y, 0.0f));
    DrawLine(RenderCommandBuffer, vec3(Planes[1].D*Planes[1].N.x, -1.0f, 0.0f), vec3(Planes[1].D*Planes[1].N.x, 1.0f, 0.0f));
    DrawLine(RenderCommandBuffer, vec3(-1.0f, Planes[2].D*Planes[2].N.y, 0.0f), vec3(1.0f, Planes[2].D*Planes[2].N.y, 0.0f));
    DrawLine(RenderCommandBuffer, vec3(Planes[3].D*Planes[3].N.x, -1.0f, 0.0f), vec3(Planes[3].D*Planes[3].N.x, 1.0f, 0.0f));
#else

    if(Intersect(&GameState->GameObjects[0], &GameState->GameObjects[1]))
    {
        PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 1.0f, 0.0f));
    }
    else
    {
        PushVec3(RenderCommandBuffer, "Color", vec3(0.0f, 0.0f, 1.0f));
    }

    game_object *A = &GameState->GameObjects[0];
    vec2 AAxisX = vec2(Cos(Radians(A->RigidBody.Orientation)), Sin(Radians(A->RigidBody.Orientation)));
    vec2 AAxisY = Perp(AAxisX);
    // TODO(georgy): This function works for any convex polygon!
    // NOTE(georgy): CCW order
    vec3 APolygonVertices[8] = 
    {
        vec3(A->RigidBody.P + 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P + 0.25f*A->Width*AAxisX + 0.75f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.25f*A->Width*AAxisX + 0.75f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.5f*A->Width*AAxisX + 0.5f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P - 0.25f*A->Width*AAxisX - 0.75f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P + 0.25f*A->Width*AAxisX - 0.75f*A->Height*AAxisY, 0.0f),
        vec3(A->RigidBody.P + 0.5f*A->Width*AAxisX - 0.5f*A->Height*AAxisY, 0.0f),
    };

    mat4 Model = Identity();
    PushMat4(RenderCommandBuffer, "Model", &Model);
    for(u32 I0 = 0, I1 = ArrayCount(APolygonVertices) - 1;
        I0 < ArrayCount(APolygonVertices);
        I1 = I0, I0++)
    {
        vec3 A = APolygonVertices[I1];
        vec3 B = APolygonVertices[I0];
        DrawLine(RenderCommandBuffer, A, B);
    }

    for(u32 GameObjectIndex = 0;
        GameObjectIndex < 2;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;
        rigid_body *RigidBody = &GameObject->RigidBody;

        RigidBody->dP += dt*((RigidBody->ForceAccumulated*(1.0f / RigidBody->Mass)) + (GravityEnabled ? vec2(0.0f, -9.8f) : vec2(0.0f, 0.0f)));
        RigidBody->P += dt*RigidBody->dP;
    }

    for(u32 GameObjectIndex = 1;
        GameObjectIndex < 2;
        GameObjectIndex++)
    {
        game_object *GameObject = GameState->GameObjects + GameObjectIndex;
        rigid_body *RigidBody = &GameObject->RigidBody;

        mat4 Model = Translation(vec3(RigidBody->P, 0.0f)) * 
                     Rotation(RigidBody->Orientation, vec3(0.0f, 0.0f, 1.0f)) * 
                     Scaling(vec3(GameObject->Width, GameObject->Height, 1.0f));
        PushMat4(RenderCommandBuffer, "Model", &Model);
        DrawModel(RenderCommandBuffer, GameObject->Model);
    }

#endif
#endif
}
